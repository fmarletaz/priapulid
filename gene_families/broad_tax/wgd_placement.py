#!/usr/bin/env python3
"""
Where on the tree does the duplication sit? Fraction of P. caudatus paralogue
pairs lying on ohnologous chromosome pairs, computed separately for each
duplication node (Figure 3a, made quantitative).

The chromosome pairing is held fixed and the node of origin is the variable. A
whole-genome duplication puts both copies on chromosomes that are themselves an
ohnologous pair; small-scale duplication does not. So paralogues assigned to the
node carrying the WGD should be WGD-consistent far more often than paralogues
assigned to any other node, and the profile across nodes localises the event.

Expected rate under chance is not flat: chromosomes hold very different numbers
of genes. For each pair the null probability is computed conditional on where
one member actually is -- if gene A sits on chromosome X, the chance its partner
lands on a chromosome ohnologous to X is the gene share of those chromosomes.
The total then follows a Poisson-binomial, giving an exact expectation and tail.

Inputs (TSV, header required)
  pairs.tsv   gene_a  gene_b  node   [family]
  chrom.tsv   gene    chrom
  ohno.tsv    chrom_a chrom_b        one line per ohnologous chromosome pair;
                                     a chromosome may appear in several pairs,
                                     which is handled as a graph rather than a
                                     strict pairing

  python3 wgd_placement.py pairs.tsv chrom.tsv ohno.tsv --out wgd_placement.tsv


"""
from __future__ import annotations
import argparse, csv, math, sys
from collections import defaultdict


# ------------------------------------------------------------------ stats
def poisson_binomial(probs):
    """pmf over the number of successes for independent, unequal trials."""
    pmf = [1.0]
    for p in probs:
        new = [0.0] * (len(pmf) + 1)
        for k, v in enumerate(pmf):
            new[k] += v * (1 - p)
            new[k + 1] += v * p
        pmf = new
    return pmf


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def fisher_greater(a, b, c, d):
    """One-sided P(X >= a) for the 2x2 [[a,b],[c,d]], hypergeometric."""
    def lc(n, k):
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    n1, n2, t = a + b, c + d, a + c
    tot = lc(n1 + n2, n1)
    p = 0.0
    for i in range(a, min(n1, t) + 1):
        v = lc(t, i) + lc(n1 + n2 - t, n1 - i) - tot
        if v > -745:
            p += math.exp(v)
    return min(1.0, p)


# ------------------------------------------------------------------ input
def read_tsv(path, need):
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = set(rd.fieldnames or [])
        missing = [c for c in need if c not in cols]
        if missing:
            sys.exit(f"{path}: missing column(s) {missing}; found {sorted(cols)}")
        return list(rd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs")
    ap.add_argument("chrom")
    ap.add_argument("ohno")
    ap.add_argument("--focal", help="node to contrast against all others; "
                                    "default is the one with most pairs")
    ap.add_argument("--min-pairs", type=int, default=10)
    ap.add_argument("--sep", default="_",
                    help="species-prefix separator, used to reconcile gene ids "
                         "between pairs.tsv and chrom.tsv when one carries the "
                         "prefix and the other does not")
    ap.add_argument("--keep-same-chrom", action="store_true",
                    help="count pairs on one chromosome in the denominator. By "
                         "default they are reported but excluded, since tandem "
                         "duplicates cannot be WGD-consistent by construction")
    ap.add_argument("--only-paired-chroms", action="store_true",
                    help="restrict to genes on chromosomes that appear in the "
                         "ohnology table, dropping unplaced scaffolds and any "
                         "sequence with no assigned partner. This is the "
                         "conservative choice: it raises the chance expectation, "
                         "because a gene confined to the main chromosomes is "
                         "much more likely to land on an ohnologous one")
    ap.add_argument("--leave-one-out", action="store_true",
                    help="see the note on circularity in the header")
    ap.add_argument("--out", default="wgd_placement.tsv")
    a = ap.parse_args()

    gene_chrom = {r["gene"]: r["chrom"] for r in read_tsv(a.chrom, ["gene", "chrom"])}

    # The two files rarely agree on whether gene ids carry the species prefix:
    # orthogroups give 'Pricau_XP_001', a GFF gives 'XP_001'. Index both forms
    # so either file can use either convention.
    def bare(g):
        return g.split(a.sep, 1)[1] if a.sep in g else g

    bare_chrom = {}
    for g, c in gene_chrom.items():
        bare_chrom.setdefault(bare(g), c)

    def chrom_of(gene):
        if gene in gene_chrom:
            return gene_chrom[gene]
        b = bare(gene)
        return gene_chrom.get(b) or bare_chrom.get(b) or bare_chrom.get(gene)
    n_genes = defaultdict(int)
    for c in gene_chrom.values():
        n_genes[c] += 1
    total_genes = sum(n_genes.values())
    print(f"[in] {len(gene_chrom)} genes on {len(n_genes)} chromosomes",
          file=sys.stderr)

    # ohnology as a graph: chromosome -> set of ohnologous partners
    partners = defaultdict(set)
    for r in read_tsv(a.ohno, ["chrom_a", "chrom_b"]):
        x, y = r["chrom_a"].strip(), r["chrom_b"].strip()
        partners[x].add(y)
        partners[y].add(x)
    unknown = [c for c in partners if c not in n_genes]
    if unknown:
        print(f"[warn] chromosome(s) in the ohnology table with no genes: "
              f"{unknown[:5]}", file=sys.stderr)
    print(f"[in] {len(partners)} chromosomes in the ohnology graph, "
          f"{sum(len(v) for v in partners.values()) // 2} pairs", file=sys.stderr)

    # Optionally confine everything to the chromosomes that carry an ohnologous
    # partner. Both the pairs AND the gene counts behind the null must be
    # restricted together: filtering the pairs alone would leave the expectation
    # computed over a gene pool the genes can no longer occupy.
    if a.only_paired_chroms:
        keep = set(partners)
        dropped_chroms = sorted(c for c in n_genes if c not in keep)
        n_dropped_genes = sum(n_genes[c] for c in dropped_chroms)
        n_genes = {c: v for c, v in n_genes.items() if c in keep}
        total_genes = sum(n_genes.values())
        print(f"[filter] restricted to {len(n_genes)} paired chromosome(s); "
              f"dropped {len(dropped_chroms)} sequence(s) carrying "
              f"{n_dropped_genes} gene(s)", file=sys.stderr)
        if total_genes == 0:
            sys.exit("[error] no genes left; do the chromosome names in the "
                     "ohnology table match those in the gene table?")

    rows = read_tsv(a.pairs, ["gene_a", "gene_b", "node"])
    print(f"[in] {len(rows)} paralogue pairs", file=sys.stderr)

    def null_p(x, y):
        """chance a partner of a gene on x lands on a chromosome ohnologous to
        x, averaged over which member is conditioned on."""
        def one(src):
            pool = total_genes - (n_genes.get(src, 0) if not a.keep_same_chrom else 0)
            if pool <= 0:
                return 0.0
            return sum(n_genes.get(z, 0) for z in partners.get(src, ())) / pool
        return 0.5 * (one(x) + one(y))

    by_node = defaultdict(list)
    n_nochrom = n_unpaired = 0
    unmatched = []
    for r in rows:
        ca, cb = chrom_of(r["gene_a"]), chrom_of(r["gene_b"])
        if ca is None or cb is None:
            n_nochrom += 1
            if len(unmatched) < 5:
                unmatched.append(r["gene_a"] if ca is None else r["gene_b"])
            continue
        if a.only_paired_chroms and (ca not in n_genes or cb not in n_genes):
            n_unpaired += 1
            continue
        same = ca == cb
        consistent = (not same) and (cb in partners.get(ca, ()))
        by_node[r["node"].strip()].append((ca, cb, same, consistent))

    if n_unpaired:
        print(f"[filter] {n_unpaired} pair(s) had a gene outside the paired "
              f"chromosomes, dropped", file=sys.stderr)
    if n_nochrom:
        print(f"[warn] {n_nochrom} of {len(rows)} pair(s) had a gene with no "
              f"chromosome, dropped", file=sys.stderr)
        if n_nochrom > 0.5 * len(rows):
            # a wholesale failure is an identifier mismatch, not missing data
            print("\n[error] most genes did not match. The two files probably "
                  "use different identifier conventions.", file=sys.stderr)
            print(f"  unmatched from pairs.tsv : {unmatched}", file=sys.stderr)
            print(f"  examples from chrom.tsv  : "
                  f"{list(gene_chrom)[:5]}", file=sys.stderr)
            print(f"  (--sep is currently {a.sep!r}; both the full name and the "
                  f"part after the first separator were tried)", file=sys.stderr)
            sys.exit(1)

    out = []
    for node, items in sorted(by_node.items(), key=lambda kv: -len(kv[1])):
        n_same = sum(1 for *_x, s, _c in items if s)
        use = items if a.keep_same_chrom else [i for i in items if not i[2]]
        n = len(use)
        if n < a.min_pairs:
            continue
        k = sum(1 for *_x, _s, c in use if c)
        ps = [null_p(ca, cb) for ca, cb, _s, _c in use]
        pmf = poisson_binomial(ps)
        exp = sum(ps)
        p_up = sum(pmf[k:]) if k <= n else 0.0
        lo, hi = wilson(k, n)
        out.append(dict(
            node=node, n_pairs=len(items), n_same_chrom=n_same, n_tested=n,
            n_consistent=k, pct=round(100 * k / n, 2),
            ci_lo=round(100 * lo, 2), ci_hi=round(100 * hi, 2),
            expected=round(exp, 1), pct_expected=round(100 * exp / n, 2),
            ratio=round(k / exp, 2) if exp > 0 else "",
            p=f"{p_up:.3g}"))

    if not out:
        sys.exit("no node reached --min-pairs")

    print(f"\n{'node':>18} {'pairs':>6} {'same':>5} {'tested':>7} {'consist':>8} "
          f"{'%':>6} {'95% CI':>14} {'%exp':>6} {'ratio':>6} {'p':>9}")
    for r in out:
        ci = f"{r['ci_lo']:.1f}-{r['ci_hi']:.1f}"
        print(f"{r['node']:>18} {r['n_pairs']:>6} {r['n_same_chrom']:>5} "
              f"{r['n_tested']:>7} {r['n_consistent']:>8} {r['pct']:>6.1f} "
              f"{ci:>14} {r['pct_expected']:>6.1f} {str(r['ratio']):>6} "
              f"{r['p']:>9}")

    # focal node against the rest, which is the comparison that matters: it asks
    # whether this node is special, not whether it beats a null model
    focal = a.focal or out[0]["node"]
    f = next((r for r in out if r["node"] == focal), None)
    if f is None:
        print(f"\n[warn] focal node '{focal}' not among the tested nodes",
              file=sys.stderr)
    else:
        rk = sum(r["n_consistent"] for r in out if r["node"] != focal)
        rn = sum(r["n_tested"] for r in out if r["node"] != focal)
        p = fisher_greater(f["n_consistent"], f["n_tested"] - f["n_consistent"],
                           rk, rn - rk)
        odds = ((f["n_consistent"] / max(f["n_tested"] - f["n_consistent"], 1)) /
                max(rk / max(rn - rk, 1), 1e-9))
        print(f"\n{focal}: {f['n_consistent']}/{f['n_tested']} "
              f"({f['pct']:.1f}%) WGD-consistent")
        print(f"all other nodes: {rk}/{rn} ({100*rk/max(rn,1):.1f}%)")
        print(f"Fisher one-sided p = {p:.3g}, odds ratio = {odds:.2f}")

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        for r in out:
            w.writerow(r)
    print(f"\n-> {a.out}")
    print("Plot: pct with ci_lo/ci_hi as error bars, pct_expected as a dashed "
          "line, nodes ordered from the root down.")


if __name__ == "__main__":
    main()