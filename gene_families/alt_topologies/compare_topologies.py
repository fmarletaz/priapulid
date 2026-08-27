#!/usr/bin/env python3
"""
Compare the gene family content inferred under two topologies.

Same orthogroups, same species, two species trees -- so any difference in the
reconstructed content of a node is caused by the topology alone. This writes a
long-format matrix of Jaccard similarities between every node of topology A and
every node of topology B, ready for geom_tile.

Two things can be compared, chosen with --what:

  present  the repertoire reconstructed AT each node: which families the
           ancestor is inferred to have had. This is the biologically
           meaningful comparison and the one to plot.
  lost     the families lost ON the branch leading to each node. More volatile,
           because a family lost once under one topology may be lost twice under
           the other, but it is the quantity Figure 4a shows.

Output (TSV, long)
  node_a node_b n_a n_b n_shared jaccard overlap expected_shared ratio same_name

  'overlap' is |A n B| / min(|A|,|B|), which unlike Jaccard is not dragged down
  when the two nodes hold very different numbers of families. 'expected_shared'
  is |A|*|B|/universe, so 'ratio' says how much more the two nodes share than
  two unrelated sets of the same sizes would.

  python3 compare_topologies.py \\
      crypto_families.tsv crypto_tree.nw tri_families.tsv tri_tree.nw \\
      --what present --name-a Cryptovermes --name-b Triradiata \\
      --out topology_content.tsv

Plot:
  library(ggplot2); library(dplyr)
  d <- read.delim("topology_content.tsv")
  ggplot(d, aes(node_b, node_a, fill = jaccard)) +
    geom_tile(colour = "white", linewidth = .3) +
    geom_text(data = subset(d, same_name == "yes"),
              aes(label = sprintf("%.2f", jaccard)), size = 2.2) +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1)) +
    labs(x = "Triradiata", y = "Cryptovermes", fill = "Jaccard") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
"""
from __future__ import annotations
import argparse, csv, sys
from ete3 import Tree


def read_families(path):
    """family_id -> (gains, losses)."""
    out = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = set(rd.fieldnames or [])
        stats = "lost_clades" in cols
        if not stats and "Losses" not in cols:
            sys.exit(f"{path}: unrecognised columns {sorted(cols)}")
        for r in rd:
            fid = r.get("family_id") or r.get("ID")
            if stats:
                g = r.get("gain_nodes") or r.get("origin") or ""
                gains = [s.strip() for s in g.split(",") if s.strip()]
                losses = {s.strip() for s in r["lost_clades"].split(",") if s.strip()}
            else:
                gains = [r["Gain"].strip()] if r["Gain"].strip() else []
                losses = {s.strip() for s in r["Losses"].split(",") if s.strip()}
            if gains:
                out[fid] = (gains, losses)
    return out


def content(tree_f, table, what):
    """node name -> set of families present at it, or lost on its branch."""
    t = Tree(tree_f, format=1)
    for i, n in enumerate(t.traverse()):
        if not n.name:
            n.name = f"node{i}"
    node = {n.name: n for n in t.traverse()}
    lineage = {nm: {x.name for x in n.get_ancestors()} | {nm}
               for nm, n in node.items()}

    fams = read_families(table)
    out = {nm: set() for nm in node}
    for fid, (gs, L) in fams.items():
        if what == "lost":
            for nm in L:
                if nm in out:
                    out[nm].add(fid)
            continue
        for nm in node:
            for g in gs:
                if g not in lineage.get(nm, ()):
                    continue
                path = lineage[nm] - lineage.get(g, set()) | {nm}
                if L & (path - {nm}):
                    continue
                out[nm].add(fid)
                break
    return out, set(fams), t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("families_a")
    ap.add_argument("tree_a")
    ap.add_argument("families_b")
    ap.add_argument("tree_b")
    ap.add_argument("--what", choices=["present", "lost"], default="present")
    ap.add_argument("--name-a", default="A")
    ap.add_argument("--name-b", default="B")
    ap.add_argument("--internal-only", action="store_true",
                    help="drop terminal branches, which carry annotation "
                         "artefacts and are identical between topologies")
    ap.add_argument("--min-size", type=int, default=1,
                    help="skip nodes holding fewer families than this")
    ap.add_argument("--out", default="topology_content.tsv")
    a = ap.parse_args()

    ca, fams_a, ta = content(a.tree_a, a.families_a, a.what)
    cb, fams_b, tb = content(a.tree_b, a.families_b, a.what)

    universe = fams_a | fams_b
    only_a, only_b = fams_a - fams_b, fams_b - fams_a
    if only_a or only_b:
        print(f"[warn] {len(only_a)} families only in {a.name_a}, "
              f"{len(only_b)} only in {a.name_b}. The two runs should use the "
              f"same orthogroups; differences here will look like topological "
              f"effects.", file=sys.stderr)

    tips_a = {l.name for l in ta.get_leaves()}
    tips_b = {l.name for l in tb.get_leaves()}
    if tips_a != tips_b:
        print(f"[warn] the two trees do not have the same species: "
              f"{sorted(tips_a ^ tips_b)[:6]}", file=sys.stderr)

    def usable(cont, tips):
        return {nm: s for nm, s in cont.items()
                if len(s) >= a.min_size and not (a.internal_only and nm in tips)}

    ua, ub = usable(ca, tips_a), usable(cb, tips_b)
    print(f"[in] {len(ua)} nodes in {a.name_a}, {len(ub)} in {a.name_b}; "
          f"universe {len(universe)} families", file=sys.stderr)

    U = len(universe)
    rows = []
    for na, sa in sorted(ua.items()):
        for nb, sb in sorted(ub.items()):
            inter = len(sa & sb)
            union = len(sa | sb)
            exp = len(sa) * len(sb) / U if U else 0
            rows.append(dict(
                node_a=na, node_b=nb, n_a=len(sa), n_b=len(sb),
                n_shared=inter,
                jaccard=round(inter / union, 4) if union else 0.0,
                overlap=round(inter / min(len(sa), len(sb)), 4)
                        if min(len(sa), len(sb)) else 0.0,
                expected_shared=round(exp, 1),
                ratio=round(inter / exp, 3) if exp > 0 else "",
                same_name="yes" if na == nb else "no"))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # the diagonal is what matters: how stable is a node's content when the
    # topology changes underneath it
    diag = [r for r in rows if r["same_name"] == "yes"]
    if diag:
        diag.sort(key=lambda r: r["jaccard"])
        print(f"\nnodes present in both topologies, least to most stable "
              f"({a.what})")
        print(f"{'node':>18} {'n_' + a.name_a[:6]:>9} {'n_' + a.name_b[:6]:>9} "
              f"{'shared':>7} {'Jaccard':>8}")
        for r in diag:
            print(f"{r['node_a']:>18} {r['n_a']:>9} {r['n_b']:>9} "
                  f"{r['n_shared']:>7} {r['jaccard']:>8.3f}")
        med = sorted(r["jaccard"] for r in diag)[len(diag) // 2]
        print(f"\nmedian Jaccard across shared nodes: {med:.3f}")
        print("A node whose content barely changes is insensitive to the "
              "topology; a low value marks a conclusion that depends on which "
              "tree is assumed.")

    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()