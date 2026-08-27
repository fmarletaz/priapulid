#!/usr/bin/env python3
"""
Re-score an existing <prefix>_besthits.tsv against a (corrected) membership map.
No hmmbuild/hmmsearch: everything is recomputed from the stored best hits.

  python3 rescore.py nematoda_besthits.tsv orthogroups.tsv
  python3 rescore.py nematoda_besthits.tsv orthogroups.tsv --sensitivity 0.90
"""
from __future__ import annotations
import argparse, csv, math, re, sys
from collections import defaultdict


def load_membership(path):
    """gene -> family.

    The orientation is decided once for the whole file, not line by line: if
    most rows list several ids after the first column, the file is a
    family->genes listing (singleton families are then read correctly too).
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 2:
                continue
            ids = [g for cell in f[1:] for g in re.split(r"[,\s]+", cell.strip()) if g]
            rows.append((f[0].strip(), ids))
    if not rows:
        return {}
    multi = sum(1 for _, ids in rows if len(ids) > 1)
    family_to_genes = multi > 0.1 * len(rows)      # a listing, not a 2-column map
    m = {}
    for first, ids in rows:
        if family_to_genes:
            for g in ids:
                m[g] = first
        else:
            m[first] = ids[0]
    return m


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


ap = argparse.ArgumentParser()
ap.add_argument("besthits")
ap.add_argument("membership")
ap.add_argument("--sensitivity", type=float, default=0.95)
ap.add_argument("--gainloss", default=None,
                help="gain/loss table (ID NbSpecies NbDesc Gain Losses); lets the "
                     "script ask WHERE the family holding a rescued gene originated")
ap.add_argument("--target", default=None, help="target node name, e.g. Nematoda")
ap.add_argument("--members", default=None,
                help="comma-separated species codes of the target clade")
ap.add_argument("--out", default=None)
a = ap.parse_args()

gain_of = {}
if a.gainloss:
    with open(a.gainloss) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gain_of[r["ID"]] = r["Gain"].strip()
    print(f"[in] gain node known for {len(gain_of)} families")
within = set()
if a.target:  within.add(a.target)
if a.members: within.update(s.strip() for s in a.members.split(","))

member = load_membership(a.membership)
print(f"[in] membership map: {len(member)} genes")

rows = []
with open(a.besthits) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.append(r)

# sanity: do the ids in the map look like the ids in the hits?
hits = [r["best_hit"] for r in rows if r["best_hit"] != "-"]
overlap = sum(1 for h in hits[:2000] if h in member)
print(f"[check] {overlap}/{min(len(hits),2000)} sampled hit ids found in the map"
      + ("" if overlap else "   <-- ids do not match, check formats"))

n_ret = sum(1 for r in rows if r["status"] == "retained")
n_lost = sum(1 for r in rows if r["status"] == "lost")
ret_scores = [float(r["bitscore"]) for r in rows
              if r["status"] == "retained" and r["bitscore"] != "-"]
thr = quantile(ret_scores, 1 - a.sensitivity)


def classify(fam, hit):
    f = member.get(hit)
    if f is None: return "unassigned"
    if f == fam:  return "same_family"
    if gain_of and within:
        g = gain_of.get(f)
        if g is None:
            return "other_family_unknown_origin"
        # a family that arose inside the target lineage is not an old paralogue:
        # it is what a diverged orthologue looks like after being reclustered
        return "other_family_lineage_specific" if g in within else "other_family_ancient"
    return "other_family"


cls_ret, cls_lost = defaultdict(int), defaultdict(int)
lost_any = lost_thr = 0
for r in rows:
    if r["best_hit"] == "-":
        continue
    c = classify(r["family"], r["best_hit"])
    if r["status"] == "retained":
        cls_ret[c] += 1
    else:
        lost_any += 1
        if float(r["bitscore"]) >= thr:
            lost_thr += 1
            cls_lost[c] += 1

# a rescued gene counts as detection failure if it was in no family at all, or in
# a family that only came into existence inside the target lineage
hdf = cls_lost["unassigned"] + cls_lost["other_family_lineage_specific"]
lo1, hi1 = wilson(len(ret_scores), n_ret)
lo2, hi2 = wilson(hdf, n_lost)

out = (
    f"\n{'':<36}{'n':>7}{'recovered':>12}{'%':>8}   95% CI\n"
    f"{'retained (leave-one-out)':<36}{n_ret:>7}{len(ret_scores):>12}"
    f"{100*len(ret_scores)/max(n_ret,1):>8.1f}   {100*lo1:.1f}-{100*hi1:.1f}\n"
    f"{'lost, any hit':<36}{n_lost:>7}{lost_any:>12}"
    f"{100*lost_any/max(n_lost,1):>8.1f}\n"
    f"{'lost, above calibrated threshold':<36}{n_lost:>7}{lost_thr:>12}"
    f"{100*lost_thr/max(n_lost,1):>8.1f}\n"
    f"\nretained hits fall in: "
    + ", ".join(f"{k} {v}" for k, v in sorted(cls_ret.items()))
    + "   <- should be dominated by same_family\n"
    f"lost hits above threshold fall in:\n"
    + "".join(f"  {k:<42}{v:>6}\n" for k, v in sorted(cls_lost.items()))
    + "  (HDF = unassigned + other_family_lineage_specific)\n"
    f"\nthreshold = {thr:.1f} bits ({100*a.sensitivity:.0f}% recovery of retained)\n"
    f"detection sensitivity = {100*len(ret_scores)/max(n_ret,1):.1f}%\n"
    f"losses attributable to HDF = {hdf} of {n_lost} "
    f"({100*hdf/max(n_lost,1):.1f}%, 95% CI {100*lo2:.1f}-{100*hi2:.1f})\n"
)
print(out)
if a.out:
    open(a.out, "w").write(out)

# rewrite the table with corrected classes
newf = a.besthits.replace(".tsv", "_reclassified.tsv")
with open(newf, "w") as fh:
    fh.write("family\tstatus\tbest_hit\tbitscore\tevalue\thit_class\thit_family\n")
    for r in rows:
        c = classify(r["family"], r["best_hit"]) if r["best_hit"] != "-" else "-"
        fh.write(f"{r['family']}\t{r['status']}\t{r['best_hit']}\t{r['bitscore']}\t"
                 f"{r['evalue']}\t{c}\t{member.get(r['best_hit'], '-')}\n")
print(f"-> {newf}")
