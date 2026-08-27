#!/usr/bin/env python3
"""
Loss rate per branch, normalised by the repertoire available to lose, with an
explicit outlier rule for calling a "gene loss episode" (Referee #3, point 10).

Raw loss counts are not comparable across branches: an ancient node inherits a
large ancestral repertoire and so accumulates more losses than a shallow one at
the same intensity. Here

    loss rate(b) = losses(b) / families present at the ancestral node of b

and a branch is called an episode if its rate exceeds Tukey's upper fence
(Q3 + 1.5 x IQR) across all branches. No molecular dating required.

Inputs
  species_tree.nw    labelled internal nodes
  family table       either layout, auto-detected:
                       ID  NbSpecies  NbDesc  Gain  Losses
                       gene_family_stats.py output (family_id, n_species,
                       n_gains, origin, gain_nodes, lost_clades, ...)

  python3 loss_rate.py species_tree.nw families_gainloss.tsv [prefix]

All outputs are named <prefix>_*, so the Cryptovermes and Triradiata runs can
sit side by side:

  python3 loss_rate.py crypto.nw crypto_families.tsv crypto
  python3 loss_rate.py triradiata.nw tri_families.tsv triradiata
"""
from __future__ import annotations
import csv, sys
from collections import defaultdict
from ete3 import Tree


def read_families(path):
    """-> list of (gains, losses, n_species). Handles both table layouts.

    gene_family_stats.py may report several independent gains per family under
    the Wagner model, so gains is always a list.
    """
    out = []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = set(rd.fieldnames or [])
        stats = "lost_clades" in cols                   # gene_family_stats.py
        if not stats and "Losses" not in cols:
            sys.exit(f"{path}: unrecognised columns {sorted(cols)}")
        print(f"[in] {'gene_family_stats.py' if stats else 'gain/loss'} table")
        for r in rd:
            if stats:
                g = r.get("gain_nodes") or r.get("origin") or ""
                gains = [s.strip() for s in g.split(",") if s.strip()]
                losses = {s.strip() for s in r["lost_clades"].split(",") if s.strip()}
                nsp = int(r["n_species"])
            else:
                gains = [r["Gain"].strip()] if r["Gain"].strip() else []
                losses = {s.strip() for s in r["Losses"].split(",") if s.strip()}
                nsp = int(r["NbSpecies"])
            if gains:
                out.append((gains, losses, nsp))
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    tree_f, table = sys.argv[1], sys.argv[2]
    prefix = sys.argv[3] if len(sys.argv) > 3 else "loss"
    out = f"{prefix}_rates.tsv"

    t = Tree(tree_f, format=1)
    for i, n in enumerate(t.traverse()):                 # make sure every node has a label
        if not n.name:
            n.name = f"node{i}"
    node = {n.name: n for n in t.traverse()}

    # ancestors of each node, by name (inclusive of itself)
    lineage = {nm: {a.name for a in n.get_ancestors()} | {nm}
               for nm, n in node.items()}

    losses = defaultdict(int)          # losses called on each branch
    present = defaultdict(int)         # families present at each node
    gains = defaultdict(int)
    unknown = set()

    for gs, L, _nsp in read_families(table):
        for g in gs:
            gains[g] += 1
        for nm in L:
            losses[nm] += 1
        unknown |= ({*gs} | L) - set(node)
        # a family is present at node N if it descends from one of its gains
        # and nothing on the path between that gain and N removed it
        for nm in node:
            for g in gs:
                if g not in lineage.get(nm, ()):
                    continue
                path = lineage[nm] - lineage.get(g, set()) | {nm}
                if L & (path - {nm}):        # lost on the way down, not at N
                    continue
                present[nm] += 1
                break                        # count the family once

    if unknown:
        print(f"[warn] {len(unknown)} node name(s) in the table are absent from "
              f"the tree, e.g. {sorted(unknown)[:5]}", file=sys.stderr)

    rows = []
    for nm, n in node.items():
        if n.is_root():
            continue
        parent = n.up.name
        avail = present.get(parent, 0)       # repertoire available on this branch
        if avail == 0:
            continue
        rows.append(dict(branch=nm, parent=parent, losses=losses.get(nm, 0),
                         available=avail, rate=losses.get(nm, 0) / avail,
                         gains=gains.get(nm, 0), is_tip=n.is_leaf()))

    # Tukey upper fence over internal branches (tips carry annotation artefacts)
    r = sorted(x["rate"] for x in rows if not x["is_tip"])
    def q(v, p):
        k = (len(v) - 1) * p; f = int(k)
        return v[f] if f + 1 >= len(v) else v[f] + (v[f+1] - v[f]) * (k - f)
    q1, q3 = q(r, .25), q(r, .75)
    fence = q3 + 1.5 * (q3 - q1)

    for x in rows:
        x["episode"] = "yes" if (x["rate"] > fence and not x["is_tip"]) else "no"

    rows.sort(key=lambda x: -x["rate"])
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for x in rows:
            w.writerow({**x, "rate": round(x["rate"], 4)})

    print(f"Q1 {q1:.3f}  Q3 {q3:.3f}  upper fence {fence:.3f} "
          f"(internal branches, n={len(r)})\n")
    print(f"{'branch':>18} {'losses':>7} {'avail':>7} {'rate':>7}  episode")
    for x in rows:
        if x["is_tip"] and x["rate"] <= fence:
            continue
        print(f"{x['branch']:>18} {x['losses']:>7} {x['available']:>7} "
              f"{x['rate']:>7.3f}  {x['episode']}")
    # ---- trees with the rate as branch length --------------------------------
    by_branch = {x["branch"]: x for x in rows}
    orig_name = {id(n): n.name for n in t.traverse()}

    def write_tree(key, path, label=None, scale=1.0):
        """label: None keeps the clade name; a column name puts that value in the
        node label instead; 'both' gives Clade|value."""
        for n in t.traverse():
            nm = orig_name[id(n)]
            b = by_branch.get(nm, {})
            n.dist = 0.0 if n.is_root() else b.get(key, 0.0) * scale
            if label is None or n.is_root():
                n.name = nm
            elif label == "both":
                v = b.get("losses", 0)
                n.name = f"{nm}|{v}" if nm else str(v)
            else:
                n.name = str(b.get(label, 0))
        res = t.write(format=1, outfile=path)
        for n in t.traverse():                    # restore
            n.name = orig_name[id(n)]
        return res

    write_tree("rate", f"{prefix}_rate.tre")
    write_tree("rate", f"{prefix}_rate_nlosses.tre", label="losses")
    write_tree("rate", f"{prefix}_rate_labelled.tre", label="both")
    write_tree("losses", f"{prefix}_counts.tre")

    with open(f"{prefix}_itol.txt", "w") as fh:   # branch labels for iTOL/FigTree
        fh.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for x in rows:
            fh.write(f"{x['branch']}\t{x['losses']}\n")

    print(f"\n-> {out}")
    print(f"-> {prefix}_rate.tre            length = losses / repertoire available")
    print(f"-> {prefix}_rate_nlosses.tre    same lengths, labels = number of losses")
    print(f"-> {prefix}_rate_labelled.tre   same lengths, labels = Clade|losses")
    print(f"-> {prefix}_counts.tre          length = raw loss count (for comparison)")
    print(f"-> {prefix}_itol.txt            branch -> loss count, for iTOL/FigTree")


if __name__ == "__main__":
    main()