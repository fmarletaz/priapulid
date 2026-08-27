#!/usr/bin/env python3
"""
Quantify homology-detection failure (HDF) among inferred gene family losses.

Profiles are built for BOTH families the target lineage retains and families it
is inferred to have lost, with the target lineage always excluded from the
alignment. Searching them against the full target proteome then gives:

  retained families -> leave-one-out recovery = detection sensitivity for that
                       lineage, measured directly and independently of any loss
                       inference
  lost families     -> candidate rescues, scored against the sensitivity
                       calibration rather than an arbitrary E-value

Multithreaded for a single large machine: hmmbuild runs in a process pool,
and the search is split into chunks of profiles run concurrently (HMMER scales
poorly past ~8 threads within one process).

Inputs
  --alignments DIR     untrimmed alignments, <family>.fa, headers <Species>_<gene>
  --families   TSV     two columns: family <TAB> retained|lost
  --proteome   FASTA   ALL predicted proteins of the target
  --exclude    CODE    species code to drop from the profiles (repeatable)

  python3 hdf_rescue.py --alignments aln/ --families families.tsv \\
      --proteome proteomes/Caeele_all.faa --exclude Caeele --exclude Trispi \\
      --jobs 100 --prefix nematoda
"""
from __future__ import annotations
import argparse, csv, math, os, shutil, subprocess, sys, tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed


# ------------------------------------------------------------------ helpers
def read_fasta(path):
    name, seq, out = None, [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: out.append((name, "".join(seq)))
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name: out.append((name, "".join(seq)))
    return out


def wilson(k, n, z=1.96):
    if n == 0: return float("nan"), float("nan")
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def quantile(v, q):
    if not v: return float("nan")
    v = sorted(v); k = (len(v) - 1) * q; f = int(k)
    return v[f] if f + 1 >= len(v) else v[f] + (v[f + 1] - v[f]) * (k - f)


def load_membership(path):
    """gene -> family. Accepts a 2-column TSV, or an OrthoFinder Orthogroups.tsv
    (wide: family in column 1, then one column per species with comma-separated
    gene ids). Genes absent from the map are treated as unassigned."""
    if not path:
        return {}
    m = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        wide = len(header) > 2
        if not wide and header and header[0].lower() not in ("gene", "sequence"):
            fh.seek(0)                                   # no header, 2 columns
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            if wide:
                fam = f[0]
                for cell in f[1:]:
                    for g in cell.split(","):
                        g = g.strip()
                        if g:
                            m[g] = fam
            else:
                m[f[0].strip()] = f[1].strip()
    return m


# ------------------------------------------------------------------ stage 1
def build_one(args):
    fam, aln_dir, hmm_dir, excl, min_seqs = args
    src = os.path.join(aln_dir, fam + ".al.fa")
    if not os.path.exists(src):
        return fam, "no_alignment"
    recs = [(h, s) for h, s in read_fasta(src) if h.split("_")[0] not in excl]
    if len(recs) < min_seqs:
        return fam, "too_few_seqs"
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for h, s in recs:
            fh.write(f">{h}\n{s}\n")
        sub = fh.name
    out = os.path.join(hmm_dir, fam + ".hmm")
    try:
        subprocess.run(["hmmbuild", "--amino", "-n", fam, out, sub],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return fam, "ok"
    except subprocess.CalledProcessError:
        return fam, "hmmbuild_failed"
    finally:
        os.unlink(sub)


# ------------------------------------------------------------------ stage 2
def search_chunk(args):
    chunk, proteome, evalue = args
    tbl = chunk + ".tbl"
    subprocess.run(["hmmsearch", "--cpu", "1", "-E", str(evalue), "--noali",
                    "--tblout", tbl, chunk, proteome],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return tbl


def chunk_profiles(hmm_dir, work, n_chunks):
    files = sorted(f for f in os.listdir(hmm_dir) if f.endswith(".hmm"))
    cdir = os.path.join(work, "chunks")
    os.makedirs(cdir, exist_ok=True)
    handles = [open(os.path.join(cdir, f"c{i:04d}.hmm"), "w") for i in range(n_chunks)]
    for i, f in enumerate(files):
        with open(os.path.join(hmm_dir, f)) as fh:
            handles[i % n_chunks].write(fh.read())
    for h in handles: h.close()
    return [h.name for h in handles if os.path.getsize(h.name) > 0], len(files)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignments", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--exclude", action="append", default=[])
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--evalue", type=float, default=1e-3,
                    help="permissive; the bitscore calibration does the real filtering")
    ap.add_argument("--min-seqs", type=int, default=4)
    ap.add_argument("--membership", default=None,
                    help="gene->family map (2-col TSV or OrthoFinder Orthogroups.tsv); "
                         "genes absent from it count as unassigned")
    ap.add_argument("--sensitivity", type=float, default=0.95,
                    help="fraction of retained families the threshold must recover")
    ap.add_argument("--prefix", default="hdf")
    ap.add_argument("--keep-tmp", action="store_true")
    a = ap.parse_args()

    for tool in ("hmmbuild", "hmmsearch"):
        if not shutil.which(tool):
            sys.exit(f"{tool} not found in PATH")

    status = {}
    with open(a.families) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if len(r) >= 2 and r[0].strip() and r[0] != "family":
                status[r[0].strip()] = r[1].strip()
    n_ret = sum(1 for v in status.values() if v == "retained")
    n_lost = sum(1 for v in status.values() if v == "lost")
    print(f"[in] {len(status)} families ({n_ret} retained, {n_lost} lost); "
          f"{a.jobs} workers")

    work = f"{a.prefix}_work"; hmm_dir = os.path.join(work, "hmm")
    os.makedirs(hmm_dir, exist_ok=True)
    excl = set(a.exclude)

    # ---- build ----------------------------------------------------------
    tally = defaultdict(int)
    jobs = [(f, a.alignments, hmm_dir, excl, a.min_seqs) for f in status]
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(build_one, j) for j in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            _, st = fut.result(); tally[st] += 1
            if i % 1000 == 0:
                print(f"  built {i}/{len(jobs)}", file=sys.stderr)
    print(f"[1/3] profiles: {dict(tally)}")

    # ---- search ---------------------------------------------------------
    chunks, nprof = chunk_profiles(hmm_dir, work, a.jobs)
    print(f"[2/3] searching {nprof} profiles vs {a.proteome} in {len(chunks)} chunks")
    tbls = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(search_chunk, (c, a.proteome, a.evalue)) for c in chunks]
        for i, fut in enumerate(as_completed(futs), 1):
            tbls.append(fut.result())
            if i % 10 == 0:
                print(f"  {i}/{len(chunks)} chunks", file=sys.stderr)

    best = {}
    for t in tbls:
        with open(t) as fh:
            for line in fh:
                if line.startswith("#"): continue
                f = line.split()
                fam, target, score, ev = f[2], f[0], float(f[5]), float(f[4])
                if fam not in best or score > best[fam][1]:
                    best[fam] = (target, score, ev)

    member = load_membership(a.membership)
    if member:
        print(f"[in] membership map: {len(member)} genes assigned to families")

    def classify(fam, hit):
        """Is the recovered gene already accounted for?"""
        if not member:
            return "unknown"
        fam_of = member.get(hit)
        if fam_of is None:
            return "unassigned"          # missed entirely -> genuine HDF
        if fam_of == fam:
            return "same_family"         # inconsistent with the inferred loss
        return "other_family"            # likely a paralogue, not HDF

    with open(f"{a.prefix}_besthits.tsv", "w") as fh:
        fh.write("family\tstatus\tbest_hit\tbitscore\tevalue\thit_class\thit_family\n")
        for fam, st in sorted(status.items()):
            if fam in best:
                t, s, e = best[fam]
                fh.write(f"{fam}\t{st}\t{t}\t{s:.1f}\t{e:.2g}\t"
                         f"{classify(fam, t)}\t{member.get(t, '-')}\n")
            else:
                fh.write(f"{fam}\t{st}\t-\t-\t-\t-\t-\n")

    # ---- score ----------------------------------------------------------
    ret_scores = [best[f][1] for f, s in status.items() if s == "retained" and f in best]
    thr = quantile(ret_scores, 1 - a.sensitivity)
    r_hit = len(ret_scores)
    r_thr = sum(1 for s in ret_scores if s >= thr)
    lost_scores = [best[f][1] for f, s in status.items() if s == "lost" and f in best]
    l_thr = sum(1 for s in lost_scores if s >= thr)

    cls = defaultdict(int)
    for f, s in status.items():
        if s == "lost" and f in best and best[f][1] >= thr:
            cls[classify(f, best[f][0])] += 1
    hdf_n = cls["unassigned"] if member else l_thr

    lo1, hi1 = wilson(r_hit, n_ret)
    lo2, hi2 = wilson(hdf_n, n_lost)
    out = (
        f"\n{'':<36}{'n':>7}{'recovered':>12}{'%':>8}   95% CI\n"
        f"{'retained (leave-one-out)':<36}{n_ret:>7}{r_hit:>12}"
        f"{100*r_hit/max(n_ret,1):>8.1f}   {100*lo1:.1f}-{100*hi1:.1f}\n"
        f"{'lost, any hit':<36}{n_lost:>7}{len(lost_scores):>12}"
        f"{100*len(lost_scores)/max(n_lost,1):>8.1f}\n"
        f"{'lost, above calibrated threshold':<36}{n_lost:>7}{l_thr:>12}"
        f"{100*l_thr/max(n_lost,1):>8.1f}\n"
    )
    if member:
        out += (
            f"\nof those above threshold, the recovered gene was:\n"
            f"  {'unassigned to any family (= HDF)':<40}{cls['unassigned']:>6}\n"
            f"  {'already in another family (paralogue)':<40}{cls['other_family']:>6}\n"
            f"  {'already in the same family (inconsistent)':<40}{cls['same_family']:>6}\n"
        )
    out += (
        f"\nthreshold = {thr:.1f} bits "
        f"({100*a.sensitivity:.0f}% recovery of retained families)\n"
        f"detection sensitivity in this lineage = {100*r_hit/max(n_ret,1):.1f}%\n"
        f"inferred losses attributable to HDF = {hdf_n} of {n_lost} "
        f"({100*hdf_n/max(n_lost,1):.1f}%, 95% CI {100*lo2:.1f}-{100*hi2:.1f})\n"
    )
    print(out)
    open(f"{a.prefix}_summary.txt", "w").write(out)
    if not a.keep_tmp:
        shutil.rmtree(os.path.join(work, "chunks"), ignore_errors=True)
    print(f"[3/3] {a.prefix}_besthits.tsv  {a.prefix}_summary.txt")


if __name__ == "__main__":
    main()
