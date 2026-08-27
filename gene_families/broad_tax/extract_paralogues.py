#!/usr/bin/env python3
"""
Paralogue pairs and the node they arose at, from GeneRax reconciled gene trees.

Feeds wgd_placement.py:  family  gene_a  gene_b  node  ...

GeneRax writes NHX tags on internal nodes: D=Y/N for duplication, S=<name> for
the species-tree node the event maps to. Those are used when present. With
--species-tree the same is derived from topology instead -- a node is a
duplication if its two child subtrees share a species, and it maps to the MRCA
in the species tree of every species below it -- so trees without tags still
work.

Input is either
  a file of 'OG_id<TAB or :>newick' lines  (the generax_*.nhx layout), or
  --dir DIR with one newick per file, the family taken from the filename.

Pair selection matters. A duplication with several genes on each side yields
many cross-pairs, and counting them all lets one large family dominate. The
default takes ONE representative pair per duplication node -- the closest gene
to the node on each side -- so every duplication contributes equally. Use
--pairs all for the full cross-product.

  python3 extract_paralogues.py generax_0707_crvm.nhx --species Pricau \\
      --out pairs.tsv

  python3 extract_paralogues.py trees.nhx --species Pricau \\
      --species-tree species.nw --out pairs.tsv
"""
from __future__ import annotations
import argparse, csv, glob, os, sys
from collections import defaultdict
from ete3 import Tree


def iter_trees(path, dirmode, pattern):
    """yield (family_id, newick string)."""
    if dirmode:
        files = sorted(glob.glob(os.path.join(path, pattern)))
        if not files:
            sys.exit(f"no files matching {pattern} in {path}")
        for f in files:
            fam = os.path.basename(f).split(".")[0]
            txt = open(f).read().strip()
            if txt:
                yield fam, txt
        return
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # 'OG_1<TAB>newick' or 'OG_1:newick'; the newick itself contains
            # colons, so split once and only on the first separator
            if "\t" in line.split("(")[0]:
                fam, _, nwk = line.partition("\t")
            else:
                fam, _, nwk = line.partition(":")
            if not nwk:
                print(f"[warn] line {ln}: no tree found, skipped", file=sys.stderr)
                continue
            yield fam.strip(), nwk.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trees")
    ap.add_argument("--dir", action="store_true",
                    help="'trees' is a directory of newick files")
    ap.add_argument("--pattern", default="*.nhx")
    ap.add_argument("--species", required=True,
                    help="the diagnostic species, e.g. Pricau")
    ap.add_argument("--sep", default="_",
                    help="species prefix separator in gene names")
    ap.add_argument("--species-tree",
                    help="labelled species tree. Duplications are mapped to the "
                         "MRCA of the species below them, so the node names in "
                         "the output are yours and match the rest of the "
                         "pipeline. Needed when the gene trees carry no S tags. "
                         "Duplication status still comes from the D tag.")
    ap.add_argument("--ignore-nhx", action="store_true",
                    help="ignore the D and S tags entirely and decide both from "
                         "topology: a node is a duplication if two of its child "
                         "subtrees share a species")
    ap.add_argument("--pairs", choices=["one", "all"], default="one")
    ap.add_argument("--allow-polytomy", action="store_true",
                    help="by default a polytomy at a duplication node is a "
                         "fatal error, since reconciled trees should be binary")
    ap.add_argument("--nodes",
                    help="keep only duplications mapping to these species-tree "
                         "nodes: comma-separated, or @file for one per line")
    ap.add_argument("--out", default="pairs.tsv")
    a = ap.parse_args()

    sp_tree = Tree(a.species_tree, format=1) if a.species_tree else None
    if sp_tree:
        for i, n in enumerate(sp_tree.traverse()):
            if not n.name:
                n.name = f"node{i}"
    if not a.nodes:
        want = None
    elif a.nodes.startswith("@"):
        want = {l.strip() for l in open(a.nodes[1:])
                if l.strip() and not l.startswith("#")}
    else:
        want = {c.strip() for c in a.nodes.split(",") if c.strip()}
    if want:
        print(f"[filter] keeping {len(want)} node(s)", file=sys.stderr)
        if sp_tree:
            absent = want - {n.name for n in sp_tree.traverse() if n.name}
            if absent:
                print(f"[warn] not in the species tree: {sorted(absent)}",
                      file=sys.stderr)

    def species_of(gene):
        return gene.split(a.sep, 1)[0] if a.sep in gene else gene

    def map_node(species):
        """species-tree node holding exactly this set of species."""
        species = {s for s in species if s in sp_leaves}
        if not species:
            return None
        if len(species) == 1:
            return next(iter(species))
        return sp_tree.get_common_ancestor(list(species)).name

    sp_leaves = set(sp_tree.get_leaf_names()) if sp_tree else set()

    rows, n_tree, n_dup, n_notag, n_poly, unmapped = [], 0, 0, 0, 0, set()
    for fam, nwk in iter_trees(a.trees, a.dir, a.pattern):
        try:
            t = Tree(nwk, format=1)
        except Exception as e:
            print(f"[warn] {fam}: {type(e).__name__}: {e}", file=sys.stderr)
            continue
        n_tree += 1

        for node in t.traverse():
            if node.is_leaf() or len(node.children) < 2:
                continue

            # every child subtree, so a polytomy is not silently truncated to
            # its first two branches
            subs = [(c, {species_of(x) for x in c.get_leaf_names()})
                    for c in node.children]
            all_sp = set().union(*[s for _c, s in subs])

            # Duplication status and node mapping are decided independently.
            # GeneRax's D tag is better than species overlap, which misses a
            # duplication whenever one copy has lost every species the other
            # kept, so the tag is used even when a species tree is supplied for
            # the mapping.
            d = None if a.ignore_nhx else getattr(node, "D", None)
            if d is not None:
                is_dup = str(d).upper().startswith("Y")
            else:
                if not a.ignore_nhx:
                    n_notag += 1
                is_dup = any(sa & sb
                             for i, (_a, sa) in enumerate(subs)
                             for _b, sb in subs[i + 1:])
            if not is_dup:
                continue
            n_dup += 1
            if len(subs) > 2:
                # a reconciled gene tree should be binary; a polytomy at a
                # duplication node means something went wrong upstream, and
                # silently splitting it across subtrees would hide that
                n_poly += 1
                msg = (f"{fam}: duplication node with {len(subs)} children "
                       f"(polytomy). Reconciled trees should be binary -- check "
                       f"the GeneRax output. Pass --allow-polytomy to form pairs "
                       f"across every pair of child subtrees instead.")
                if not a.allow_polytomy:
                    sys.exit("[error] " + msg)
                print("[warn] " + msg, file=sys.stderr)

            # where does it map on the species tree?
            # mapping: a supplied species tree wins, since its node names are
            # the ones every other script in the pipeline uses
            if sp_tree:
                where = map_node(all_sp)
            else:
                s = None if a.ignore_nhx else getattr(node, "S", None)
                where = str(s) if s is not None else None
            if where is None:
                unmapped.add(fam)
                continue
            if want and where not in want:
                continue

            # genes of the diagnostic species in each child subtree. Pairs are
            # formed ACROSS subtrees only: two genes within the same subtree
            # coalesce at a deeper node, which is scored separately with its own
            # label. (A,B),(C,D) therefore gives A-C, A-D, B-C, B-D, not A-B.
            side = []
            for c, _s in subs:
                g = [x for x in c.get_leaf_names() if species_of(x) == a.species]
                if g:
                    side.append((c, g))
            if len(side) < 2:
                continue          # not informative for this species

            chosen = []
            for i in range(len(side)):
                for j in range(i + 1, len(side)):
                    (ca, ga), (cb, gb) = side[i], side[j]
                    if a.pairs == "all":
                        chosen += [(x, y, len(ga), len(gb)) for x in ga for y in gb]
                    else:
                        # closest gene to the node on each side, so a large
                        # family does not outweigh a small one
                        def pick(sub, cand):
                            try:
                                return min(cand, key=lambda g: sub.get_distance(g))
                            except Exception:
                                return cand[0]
                        chosen.append((pick(ca, ga), pick(cb, gb),
                                       len(ga), len(gb)))

            for x, y, na, nb in chosen:
                rows.append(dict(family=fam, gene_a=x, gene_b=y, node=where,
                                 n_a=na, n_b=nb, n_species=len(all_sp)))

    if not rows:
        sys.exit("no paralogue pairs found; check --species and --sep")

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    per_node = defaultdict(int)
    for r in rows:
        per_node[r["node"]] += 1
    print(f"[in] {n_tree} trees, {n_dup} duplication nodes", file=sys.stderr)
    if n_notag:
        print(f"[warn] {n_notag} node(s) had no D tag; duplication status "
              f"inferred from species overlap", file=sys.stderr)
    if n_poly:
        print(f"[note] {n_poly} duplication node(s) were polytomies; pairs were "
              f"formed across every pair of child subtrees", file=sys.stderr)
    if unmapped:
        print(f"[warn] {len(unmapped)} family(ies) had a duplication that could "
              f"not be mapped to a species-tree node", file=sys.stderr)
    print(f"\n{len(rows)} pairs for {a.species}, by node:")
    for nm, c in sorted(per_node.items(), key=lambda kv: -kv[1]):
        print(f"{nm:>20} {c:>7}")
    if want:
        never = sorted(want - set(per_node))
        if never:
            print(f"\n[warn] requested but yielded no pairs: {never}. Check the "
                  f"spelling against the node names in the reconciliation.",
                  file=sys.stderr)
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()