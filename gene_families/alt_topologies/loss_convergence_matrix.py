#!/usr/bin/env python3
"""
All-against-all overlap of gene family losses within ONE reconstruction, as a
long-format matrix for a heatmap.

Output columns match compare_topologies.py, so the same plotting code works.

Every cell carries a 'relationship' label, because two of the three kinds are
constrained by parsimony and must not be read as biology:

  sister    the two branches descend from the same parent. A family absent from
            both is explained by ONE loss at that parent, so shared losses are
            structurally impossible and the cell is ~0 by construction. Grey
            these out rather than interpreting them.
  nested    one node is an ancestor of the other. A family lost at the ancestor
            is gone before the descendant exists, so again the overlap is
            forced, this time in the other direction.
  disjoint  neither of the above: the only cells where a high overlap can mean
            convergence.

'ratio' is the observed overlap over what two independent branches with the
same loss counts would share. It is the number to plot, not the raw count --
Nematoda and Tardigrada lose 38% and 47% of the repertoire, so they share
thousands of families whatever happens. Even so, ratio is inflated for distant
pairs by the same parsimony structure described above, so compare cells with
each other rather than against 1.

  python3 loss_convergence_matrix.py families_gainloss.tsv species_tree.nw \\
      --nodes Priapulida,Nematoda,Tardigrada,Chelicerata,Myriapoda,Pancrustacea \\
      --universe Ecdysozoa --out loss_convergence.tsv

Plot (identical to the topology heatmap):
  library(ggplot2); library(dplyr)
  d <- read.delim("loss_convergence.tsv") %>%
       mutate(ratio = ifelse(relationship == "disjoint", ratio, NA))
  ggplot(d, aes(node_b, node_a, fill = ratio)) +
    geom_tile(colour = "white", linewidth = .3) +
    geom_text(aes(label = ifelse(is.na(ratio), "", sprintf("%.2f", ratio)),
                  colour = ratio > 2.5), size = 2, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
    scale_fill_distiller(palette = "PuBu", direction = 1, na.value = "grey93") +
    scale_y_discrete(limits = rev) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
"""
from __future__ import annotations
import argparse, csv, sys
from collections import defaultdict
from itertools import combinations
from ete3 import Tree


def read_families(path):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("families")
    ap.add_argument("tree")
    ap.add_argument("--nodes", help="comma-separated; default = all internal "
                                    "nodes with >= --min-losses")
    ap.add_argument("--min-losses", type=int, default=100)
    ap.add_argument("--universe", default="intersection",
                    help="'intersection' (per pair), 'all', or a node name")
    ap.add_argument("--out", default="loss_convergence.tsv")
    a = ap.parse_args()

    t = Tree(a.tree, format=1)
    for i, n in enumerate(t.traverse()):
        if not n.name:
            n.name = f"node{i}"
    node = {n.name: n for n in t.traverse()}
    lineage = {nm: {x.name for x in n.get_ancestors()} | {nm}
               for nm, n in node.items()}

    fams = read_families(a.families)
    present, lost = defaultdict(set), defaultdict(set)
    for fid, (gs, L) in fams.items():
        for nm in L:
            if nm in node:
                lost[nm].add(fid)
        for nm in node:
            for g in gs:
                if g not in lineage.get(nm, ()):
                    continue
                path = lineage[nm] - lineage.get(g, set()) | {nm}
                if L & (path - {nm}):
                    continue
                present[nm].add(fid)
                break

    if a.nodes:
        nodes = [c.strip() for c in a.nodes.split(",") if c.strip()]
        missing = [c for c in nodes if c not in node]
        if missing:
            sys.exit(f"not in the tree: {missing}")
    else:
        nodes = sorted(nm for nm, n in node.items()
                       if not n.is_root() and not n.is_leaf()
                       and len(lost[nm]) >= a.min_losses)
    print(f"[in] {len(nodes)} nodes, {len(fams)} families", file=sys.stderr)

    def relationship(x, y):
        if x == y:
            return "self"
        if x in lineage[y] or y in lineage[x]:
            return "nested"
        px, py = node[x].up, node[y].up
        if px is not None and py is not None and px is py:
            return "sister"
        return "disjoint"

    def universe_for(x, y):
        if a.universe == "all":
            return set(fams)
        if a.universe != "intersection":
            return present[a.universe]
        ux = present[node[x].up.name] if node[x].up else set(fams)
        uy = present[node[y].up.name] if node[y].up else set(fams)
        return ux & uy

    if a.universe not in ("intersection", "all") and a.universe not in node:
        sys.exit(f"--universe: '{a.universe}' is not a node in the tree")

    rows = []
    for x, y in [(x, y) for x in nodes for y in nodes]:
        U = universe_for(x, y)
        sa, sb = lost[x] & U, lost[y] & U
        inter, union = len(sa & sb), len(sa | sb)
        exp = len(sa) * len(sb) / len(U) if U else 0
        rows.append(dict(
            node_a=x, node_b=y, n_a=len(sa), n_b=len(sb), n_shared=inter,
            n_universe=len(U),
            jaccard=round(inter / union, 4) if union else 0.0,
            overlap=round(inter / min(len(sa), len(sb)), 4)
                    if min(len(sa), len(sb)) else 0.0,
            expected_shared=round(exp, 1),
            ratio=round(inter / exp, 3) if exp > 0 else "",
            relationship=relationship(x, y)))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    dis = [r for r in rows if r["relationship"] == "disjoint"
           and r["node_a"] < r["node_b"]]
    dis.sort(key=lambda r: -(r["ratio"] if r["ratio"] != "" else 0))
    print(f"\ninterpretable pairs only (disjoint), most to least shared")
    print(f"{'pair':>36} {'n_a':>6} {'n_b':>6} {'shared':>7} {'exp':>7} "
          f"{'ratio':>6} {'Jacc':>6}")
    for r in dis:
        print(f"{r['node_a'] + ' & ' + r['node_b']:>36} {r['n_a']:>6} "
              f"{r['n_b']:>6} {r['n_shared']:>7} {r['expected_shared']:>7.0f} "
              f"{str(r['ratio']):>6} {r['jaccard']:>6.3f}")

    n_sis = sum(1 for r in rows if r["relationship"] == "sister")
    n_nest = sum(1 for r in rows if r["relationship"] == "nested")
    print(f"\n{len(dis)} interpretable pairs; {n_sis//2} sister and "
          f"{n_nest//2} nested pairs excluded as parsimony-constrained.")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()