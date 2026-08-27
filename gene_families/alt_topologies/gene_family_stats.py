#!/usr/bin/env python3
"""
gene_family_stats.py -- Dollo gain/loss analysis and summary statistics for
gene families, given an orthogroup clustering and a species tree.

For every gene family the script reports:

  * basic size statistics (genes, species, copy number distribution)
  * the inferred family origin: the MRCA of the species in which the family
    is present (Dollo parsimony -- a single gain, losses thereafter)
  * secondary losses: the maximal clades within the origin subtree in which
    the family is entirely absent
  * optionally, species in which the family is significantly expanded
    relative to the other members of the family (one-sided Fisher exact test,
    Benjamini-Hochberg corrected across families within each species)

Optional secondary outputs: a per-node summary of gains/losses, the species
tree annotated with those counts in Newick, and a family x species count
matrix.

Input formats are auto-detected:

  broccoli     two tab-separated columns, "family_id <TAB> gene gene gene",
               with species encoded in the gene name (Species_geneid)
  orthomcl     OrthoFinder Orthogroups.txt, "OG0000000: gene gene gene"
  orthofinder  OrthoFinder Orthogroups.tsv, one column per species,
               comma-separated gene lists

Requires: ete3, numpy, scipy, statsmodels, pandas.

All ete3 calls are confined to the "Tree operations" section below, so a
future port to ete4 only needs to touch that block. The ete4 equivalents are
noted in comments there.

Examples
--------
    ./gene_family_stats.py -g orthologous_groups_0623.txt -t metazoa-5.tre \
        -o families_gainloss.tsv

    ./gene_family_stats.py -g Orthogroups.tsv -t metazoa-5.tre \
        -o families_gainloss.tsv --node-summary nodes.tsv \
        --annotate-tree metazoa --counts-matrix counts.tsv \
        --exclude Hofmia Epibro

    ./gene_family_stats.py -g clusters.txt -t tree.tre -o out.tsv --no-expansion
"""

from __future__ import annotations

import argparse
import copy
import csv
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from ete3 import PhyloTree
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# ==========================================================================
# Tree operations (the only ete3-dependent code)
# ==========================================================================

def load_tree(path: str):
    """Load the species tree, label unnamed internal nodes, and return
    (tree, {leaf_name: node}, n_relabelled).

    ete4: PhyloTree(open(path), parser=1); node.is_leaf is a property.
    """
    with open(path) as fh:
        tree = PhyloTree(fh.read(), format=1)

    relabelled = 0
    seen = set()
    for i, node in enumerate(tree.traverse("preorder")):
        if not node.name:
            if node.is_leaf():
                raise ValueError("species tree contains an unnamed tip")
            node.name = "n{}".format(i)
            relabelled += 1
        elif node.name in seen and not node.is_leaf():
            node.name = "{}.{}".format(node.name, i)
        seen.add(node.name)

    leaves = {l.name: l for l in tree.get_leaves()}
    if len(leaves) != len(tree.get_leaves()):
        raise ValueError("species tree contains duplicated tip labels")
    return tree, leaves, relabelled


def restrict_tree(tree, keep: set[str]):
    """Prune the tree to `keep`, preserving internal node labels.

    ete3's prune() collapses single-child nodes for us. Returns a copy so the
    caller's tree is untouched.
    """
    pruned = tree.copy()
    pruned.prune(sorted(keep), preserve_branch_length=True)
    leaves = {l.name: l for l in pruned.get_leaves()}
    missing = keep - set(leaves)
    if missing:
        raise ValueError("pruning lost species: {}".format(", ".join(sorted(missing))))
    return pruned, leaves


def set_presence(leaves: dict, present: set[str]) -> None:
    """Mark every tip as present ('1') or absent ('0').

    The notebook deep-copied the whole tree once per family, which dominates
    the runtime on tens of thousands of families. Re-stamping the property in
    place gives identical results, since the previous family's values are
    fully overwritten each time.

    ete4: node.add_props(presence=...) and node.props['presence'] to read.
    """
    for name, leaf in leaves.items():
        leaf.add_features(presence=("1" if name in present else "0"))


def find_origin(tree, leaves: dict, present: set[str]):
    """MRCA of the species in which the family is present.

    ete4: tree.common_ancestor(...)
    """
    names = sorted(present)
    if len(names) == 1:
        return leaves[names[0]]
    return tree.get_common_ancestor(names)


def absent_clades(origin) -> list:
    """Maximal clades below `origin` in which the family is wholly absent.

    ete4: origin.get_monophyletic(values=['0'], prop='presence')
    """
    if origin.is_leaf():
        return []
    return list(origin.get_monophyletic(values=["0"], target_attr="presence"))


# -- Asymmetric Wagner parsimony -------------------------------------------
#
# Dollo forbids re-gain: a family arises once, at the MRCA of the species that
# have it, and every absence below that point must be a loss. Losses are free,
# so the origin is always as deep as the data allow, and a single spurious gene
# in a distant lineage buys an arbitrarily deep origin plus a cascade of
# invented losses.
#
# Wagner parsimony (Farris 1970) instead minimises the total number of state
# changes, so a patchy distribution can be explained by several independent
# gains rather than one deep gain and many losses. The asymmetric version
# (Csuros 2008) charges `gain_cost` per gain against 1 per loss: at gain_cost=1
# this is plain Wagner, and as gain_cost grows it converges on Dollo. Sweeping
# the ratio therefore measures how much of an origin call is data and how much
# is Dollo's assumption.
#
# Implemented as two-state Sankoff: a postorder pass accumulating the cost of
# each state at each node, then a preorder traceback fixing states.

_INF = float("inf")


def wagner_reconstruct(tree, present: set[str], gain_cost: float,
                       loss_cost: float = 1.0, eps: float = 1e-9):
    """Return (gain_nodes, loss_nodes) under asymmetric Wagner parsimony.

    A gain node is one whose state is 'present' while its parent's is 'absent'
    (the root counts as a gain if present, since the family has to arise
    somewhere). A loss node is the converse: absent below a present parent,
    i.e. the stem of the lost clade, matching the Dollo convention above.

    `eps` is added to each gain so that ties are broken toward *fewer*
    independent gains. Without it, ratios sitting exactly on a threshold
    resolve arbitrarily, and the sweep looks noisy for no real reason.
    """
    cost: dict[int, tuple[float, float]] = {}

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            here = node.name in present
            cost[id(node)] = ((_INF, 0.0) if here else (0.0, _INF))
            continue
        c_absent = c_present = 0.0
        for child in node.children:
            ch_absent, ch_present = cost[id(child)]
            # parent absent: child absent (no event) or child present (a gain)
            c_absent += min(ch_absent, ch_present + gain_cost + eps)
            # parent present: child present (no event) or child absent (a loss)
            c_present += min(ch_present, ch_absent + loss_cost)
        cost[id(node)] = (c_absent, c_present)

    root_absent, root_present = cost[id(tree)]
    root_state = (root_present + gain_cost + eps) < root_absent

    gains, losses = ([tree] if root_state else []), []

    def descend(node, state: bool) -> None:
        for child in node.children:
            ch_absent, ch_present = cost[id(child)]
            if state:
                ch_state = ch_present <= ch_absent + loss_cost
                if not ch_state:
                    losses.append(child)
            else:
                ch_state = (ch_present + gain_cost + eps) < ch_absent
                if ch_state:
                    gains.append(child)
            descend(child, ch_state)

    descend(tree, root_state)
    return gains, losses


def reconstruct(tree, leaves, present: set[str], model: str, gain_cost: float):
    """Uniform interface over both models: (gain_nodes, loss_nodes)."""
    if model == "dollo":
        set_presence(leaves, present)
        origin = find_origin(tree, leaves, present)
        return [origin], absent_clades(origin)
    return wagner_reconstruct(tree, present, gain_cost)


def descendant_names(node) -> list[str]:
    """ete4: list(node.leaf_names())"""
    return node.get_leaf_names()


def origin_support(origin, present: set[str]) -> int:
    """Smallest number of present species on any side of the origin split.

    A value of 1 means the origin rests on a single gene in one lineage: drop
    that gene and the origin moves tipwards, taking its whole cascade of
    inferred losses with it. This is the main way an orthogroup false positive
    turns into a deep origin plus spurious losses, so it is worth carrying as
    a column rather than discovering it later.
    """
    if origin.is_leaf():
        return len(present)
    return min(len(present.intersection(descendant_names(child)))
               for child in origin.children
               if present.intersection(descendant_names(child)))


def annotated_newick(tree, counts: dict) -> str:
    """Copy of the tree with each node label suffixed by its count.

    ete4: tree.write(parser=1)
    """
    out = copy.deepcopy(tree)
    for node in out.traverse("preorder"):
        node.name = "{}_{}".format(node.name, counts.get(node.name, 0))
    return out.write(format=1)


# ==========================================================================
# Orthogroup input
# ==========================================================================

def sniff_format(path: str) -> str:
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                first = line.rstrip("\n")
                break
        else:
            raise ValueError("{} appears to be empty".format(path))

    head = [h.strip().lower() for h in first.split("\t")]
    if head[0] in {"orthogroup", "og", "hog", "family"} and len(head) >= 2:
        return "orthofinder"
    # A table round-tripped through pandas or R loses the "Orthogroup" header
    # and starts with an empty index field instead. Broccoli output never has
    # an empty first field, so this is unambiguous.
    if head[0] == "" and len(head) >= 3:
        return "orthofinder"
    if ":" in head[0] and len(head) == 1:
        return "orthomcl"
    return "broccoli" if len(head) >= 2 else "orthomcl"


def read_species_map(path: str) -> dict[str, str]:
    """Two-column TSV mapping names as they appear in the orthogroup file to
    names as they appear in the tree. Blank lines and #comments are skipped."""
    mapping = {}
    with open(path, encoding="utf-8-sig") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 2 or not parts[0] or not parts[1]:
                raise ValueError(
                    "{} line {}: expected two tab-separated names, got {!r}".format(
                        path, lineno, line))
            mapping[parts[0]] = parts[1]
    return mapping


def apply_species_map(fams: dict, mapping: dict) -> dict:
    """Rename species in place using `mapping`, leaving unlisted names alone.

    OrthoFinder names its columns after the input filenames, so a table can
    say Aiptasia_pallida where the tree says Apal. Renaming here rather than
    editing the tree keeps the tree as the canonical reference.
    """
    if not mapping:
        return fams
    seen, renamed = set(), set()
    out = {}
    for fid, gbs in fams.items():
        new = {}
        for sp, genes in gbs.items():
            seen.add(sp)
            target = mapping.get(sp, sp)
            if target != sp:
                renamed.add(sp)
            # Two source names can legitimately map onto one tree tip.
            new.setdefault(target, []).extend(genes)
        out[fid] = new
    unused = set(mapping) - seen
    if unused:
        print("[species] warning: {} entries in the species map never matched: "
              "{}".format(len(unused), ", ".join(sorted(unused))), file=sys.stderr)
    print("[species] renamed {} of {} species via the map".format(
        len(renamed), len(seen)), file=sys.stderr)
    return out


def species_of(gene: str, sep: str, where: str) -> tuple[str, str]:
    if where == "prefix":
        sp, _, gid = gene.partition(sep)
        return sp, (gid or gene)
    idx = gene.rfind(sep)
    return (gene, gene) if idx == -1 else (gene[idx + 1:], gene[:idx])


# Columns in OrthoFinder's hierarchical orthogroup files (N0.tsv and friends)
# that are metadata rather than species.
_HOG_METADATA_COLUMNS = {"hog", "og", "gene tree parent clade"}


def read_unassigned(path: str, prefix: str = "UNASSIGNED") -> dict:
    """Read OrthoFinder's Orthogroups_UnassignedGenes.tsv as singleton families.

    OrthoFinder excludes genes with no homolog anywhere from Orthogroups.tsv
    and writes them here instead. Those genes are exactly the species-specific
    families that would be gained at tips, so leaving them out systematically
    undercounts recent gains -- while putting them in means every mispredicted
    ORF, unannotated fragment and contaminant becomes a "family gained at this
    tip". Neither choice is neutral, so families from this file are given ids
    prefixed with `prefix` and flagged in the output table, letting the
    downstream analysis include or exclude them explicitly.

    The layout is the same as Orthogroups.tsv: one row per orthogroup, one
    column per species. Each row here holds a single gene, but the parser does
    not rely on that.
    """
    fams: dict[str, dict[str, list[str]]] = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rdr = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = [h.strip() for h in next(rdr)]
        first_species = 1
        while (first_species < len(header)
               and header[first_species].lower() in _HOG_METADATA_COLUMNS):
            first_species += 1
        spnames = header[first_species:]
        for row in rdr:
            if not row or not row[0].strip():
                continue
            gbs = {}
            for sp, cell in zip(spnames, row[first_species:]):
                genes = [g.strip() for g in cell.split(",") if g.strip()]
                if genes:
                    gbs[sp] = genes
            if gbs:
                fams["{}_{}".format(prefix, row[0].strip())] = gbs

    n_multi = sum(1 for g in fams.values() if len(g) > 1)
    print("[unassigned] {} singleton families from {}{}".format(
        len(fams), path,
        "" if not n_multi else
        " (warning: {} span more than one species, which is unexpected "
        "here)".format(n_multi)), file=sys.stderr)
    return fams


def read_orthogroups(path: str, fmt: str, sep: str, where: str,
                     keep_full_names: bool = False) -> dict:
    """Return {family_id: {species: [gene, ...]}}.

    By default the species tag is stripped from gene identifiers, since the
    species is already the key. Pass keep_full_names=True to retain them as
    written, which is what downstream joins against external annotation files
    (eggNOG, InterProScan) need in order to match.

    Species with no genes are never stored. This matters: in the notebook
    Fams[fid] was a defaultdict(list), so a bare lookup such as gbs['Pricau']
    inserted an empty list as a side effect, and checkFam's
    len(gbs.keys()) then counted that phantom species. Plain dicts here.
    """
    fams: dict[str, dict[str, list[str]]] = {}

    if fmt == "orthofinder":
        # utf-8-sig strips a BOM if one is present; QUOTE_NONE because gene
        # identifiers derived from FASTA headers may contain quote characters,
        # which csv would otherwise treat as field delimiters and silently
        # merge columns across.
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rdr = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
            header = [h.strip() for h in next(rdr)]
            # N0.tsv and other hierarchical outputs carry several metadata
            # columns before the species; Orthogroups.tsv carries exactly one.
            first_species = 1
            while (first_species < len(header)
                   and header[first_species].lower() in _HOG_METADATA_COLUMNS):
                first_species += 1
            spnames = header[first_species:]
            if not spnames:
                raise ValueError(
                    "{}: no species columns found after the header fields {}".format(
                        path, ", ".join(header[:first_species]) or "(none)"))
            print("[input] {} species columns: {}".format(
                len(spnames), ", ".join(spnames)), file=sys.stderr)

            for lineno, row in enumerate(rdr, first_species + 1):
                if not row or not row[0].strip():
                    continue
                cells = row[first_species:]
                if len(cells) > len(spnames):
                    raise ValueError(
                        "{} line {}: {} fields but the header declares {}. A gene "
                        "identifier containing a tab would do this.".format(
                            path, lineno, len(row), len(header)))
                # Short rows are normal: some writers trim trailing empty fields.
                gbs = {}
                for sp, cell in zip(spnames, cells):
                    genes = [g.strip() for g in cell.split(",") if g.strip()]
                    if genes:
                        gbs[sp] = genes
                if gbs:
                    fams[row[0].strip()] = gbs
        if not fams:
            raise ValueError("{}: parsed 0 families -- is this really an "
                             "OrthoFinder table?".format(path))
        return fams

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if fmt == "orthomcl":
                fid, _, rest = line.partition(":")
            else:
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                fid, rest = parts[0], "\t".join(parts[1:])
            gbs = defaultdict(list)
            for gene in rest.split():
                sp, gid = species_of(gene, sep, where)
                gbs[sp].append(gene if keep_full_names else gid)
            if gbs:
                fams[fid.strip()] = dict(gbs)
    return fams


# ==========================================================================
# Expansion test
# ==========================================================================

def expansion_pvalues(counts: dict, totals: dict, universe: str) -> dict:
    """One-sided Fisher exact test per species: is this species' copy number
    in this family high relative to the *median* member of the family, given
    its genome-wide gene count?

    universe='present' compares against the other species with at least one
    gene in the family (the notebook's behaviour). universe='all' counts
    absent species as zero, which is stricter for patchily-distributed
    families.

    The notebook also computed a second test against the summed rest of the
    family and then discarded its p-value; only the median-based test is
    kept here. As in the notebook the test is one-sided, so contractions are
    not detected.
    """
    cbs = ({sp: counts.get(sp, 0) for sp in totals} if universe == "all"
           else {sp: c for sp, c in counts.items() if c > 0})

    pvals = {}
    for sp, obs in cbs.items():
        if obs == 0:
            continue
        others = [c for s, c in cbs.items() if s != sp]
        oth_tot = [t for s, t in totals.items() if s != sp]
        if not others or not oth_tot:
            continue
        ref = int(round(np.median(others)))
        ref_tot = int(round(np.median(oth_tot)))
        table = [[obs, max(totals[sp] - obs, 0)], [ref, max(ref_tot - ref, 0)]]
        pvals[sp] = fisher_exact(table, alternative="greater")[1]
    return pvals


# ==========================================================================
# Main analysis
# ==========================================================================

def analyse(fams, tree, leaves, args):
    totals = defaultdict(int)
    for gbs in fams.values():
        for sp, genes in gbs.items():
            totals[sp] += len(genes)
    totals = dict(totals)

    # The notebook wrapped checkFam in a bare `except`, which on failure left
    # the previous family's nspec/ndesc/oritax/lost still bound -- and those
    # stale values were written into idGainLoss[fid] on the very next line.
    # Rather than swap that for a narrower catch, assert the invariant that
    # made it necessary: every species in every family is a tip of the tree.
    # Species sets are reconciled in main(), so a violation here is a bug in
    # this script, not bad input, and should be loud.
    stray = {sp for gbs in fams.values() for sp in gbs} - set(leaves)
    if stray:
        raise AssertionError(
            "species absent from the tree survived reconciliation: {}".format(
                ", ".join(sorted(stray))))

    origins, losses = defaultdict(int), defaultdict(int)
    rows = []

    for fid, gbs in fams.items():
        present = {sp for sp, genes in gbs.items() if genes}
        if not present:
            continue

        gains, lost = reconstruct(tree, leaves, present, args.model, args.gain_cost)
        # Under Wagner a family can arise more than once, so the clade it is
        # "expected" in is the union of the subtrees below every gain.
        expected = {n for g in gains for n in descendant_names(g)}

        counts = np.array([len(g) for g in gbs.values() if g])
        for g in gains:
            origins[g.name] += 1
        for nd in lost:
            losses[nd.name] += 1

        single = gains[0] if len(gains) == 1 else None
        rows.append({
            "family_id": fid,
            "from_unassigned": bool(args.unassigned) and fid.startswith(
                args.unassigned_prefix + "_"),
            "n_genes": int(counts.sum()),
            "n_species": len(present),
            "n_species_in_clade": len(expected),
            "n_species_missing": len(expected) - len(present),
            "occupancy": round(len(present) / len(expected), 4),
            "mean_copies": round(float(counts.mean()), 3),
            "median_copies": float(np.median(counts)),
            "max_copies": int(counts.max()),
            "single_copy": bool(counts.max() == 1),
            "n_gains": len(gains),
            "origin": single.name if single else "multiple",
            "gain_nodes": ",".join(n.name for n in gains),
            "origin_is_tip": bool(single.is_leaf()) if single else False,
            "origin_min_side": origin_support(single, present) if single else np.nan,
            "n_loss_events": len(lost),
            "n_species_lost": sum(len(descendant_names(nd)) for nd in lost),
            "lost_clades": ",".join(nd.name for nd in lost),
        })

    df = pd.DataFrame(rows).set_index("family_id")

    if args.no_expansion:
        df["expanded_species"] = ""
        df["expansion_min_q"] = np.nan
        return df, dict(origins), dict(losses), totals

    per_sp_fam, per_sp_p = defaultdict(list), defaultdict(list)
    tested = 0
    for fid in df.index:
        counts = {sp: len(g) for sp, g in fams[fid].items() if g}
        if len(counts) < args.min_species or sum(counts.values()) < args.min_genes:
            continue
        tested += 1
        for sp, p in expansion_pvalues(counts, totals, args.expansion_universe).items():
            per_sp_fam[sp].append(fid)
            per_sp_p[sp].append(p)

    expanded, min_q = defaultdict(list), {}
    for sp, plist in per_sp_p.items():
        reject, qvals = multipletests(pvals=np.array(plist), alpha=args.alpha,
                                      method="fdr_bh")[:2]
        for fid, sig, qv in zip(per_sp_fam[sp], reject, qvals):
            if sig:
                expanded[fid].append("{}:{:.2e}".format(sp, qv))
                if fid not in min_q or qv < min_q[fid]:
                    min_q[fid] = qv

    df["expanded_species"] = pd.Series(
        {f: ",".join(v) for f, v in expanded.items()}).reindex(df.index).fillna("")
    df["expansion_min_q"] = pd.Series(min_q).reindex(df.index)
    print("[expansion] tested {} families (>= {} species, >= {} genes); {} with at "
          "least one expanded species at FDR < {}".format(
              tested, args.min_species, args.min_genes, len(expanded), args.alpha),
          file=sys.stderr)

    return df, dict(origins), dict(losses), totals


def report_size_distribution(fams, stream=sys.stderr) -> None:
    """Family size quantiles, plus a warning on outsized clusters.

    Single-linkage clustering (Broccoli in particular) chains unrelated gene
    families into giant components. A giant family is present in nearly every
    species, so its origin is the root and it contributes no losses, but it
    inflates the gene counts and quietly removes real families from the
    analysis by absorbing them.
    """
    sizes = np.array([sum(len(g) for g in gbs.values()) for gbs in fams.values()])
    total = int(sizes.sum())
    qs = np.percentile(sizes, [50, 90, 99, 99.9])
    print("[sizes] genes per family: median {:.0f}, p90 {:.0f}, p99 {:.0f}, "
          "p99.9 {:.0f}, max {}".format(*qs, sizes.max()), file=stream)

    if len(fams) < 50:
        return                      # the 1% rule below is meaningless here

    big = sorted(((sum(len(g) for g in gbs.values()), fid)
                  for fid, gbs in fams.items()), reverse=True)[:5]
    # Both tests must fire: >1% of all genes, and far above the typical family.
    # Either alone gives false positives -- the first on small datasets, the
    # second whenever the median family is tiny.
    flagged = [(n, fid) for n, fid in big
               if n > 0.01 * total and n > 50 * max(qs[0], 1)]
    if flagged:
        print("[sizes] warning: {} family/families hold >1% of all genes each, "
              "which usually means single-linkage chaining:".format(len(flagged)),
              file=stream)
        for n, fid in flagged:
            print("          {}\t{} genes ({:.1f}% of total)".format(
                fid, n, 100 * n / total), file=stream)


def run_sweep(path, fams, tree, leaves, gain_costs, stream=sys.stderr):
    """Recompute origins across a range of gain:loss cost ratios.

    Dollo is the limit of asymmetric Wagner as the gain cost grows, so the
    sweep separates origins that the data support from origins that are an
    artefact of Dollo charging nothing for a loss. A family whose origin is
    stable across the whole range can be claimed; one that collapses tipwards
    at a modest gain penalty is resting on Dollo's assumption.
    """
    per_family: dict[str, list[str]] = {}
    for fid, gbs in fams.items():
        present = {sp for sp, genes in gbs.items() if genes}
        if not present:
            continue
        row = []
        for g in gain_costs:
            gains, _ = wagner_reconstruct(tree, present, g)
            row.append(gains[0].name if len(gains) == 1 else
                       "multiple:" + ",".join(sorted(n.name for n in gains)))
        # Dollo reference, i.e. the gain cost -> infinity limit
        set_presence(leaves, present)
        row.append(find_origin(tree, leaves, present).name)
        per_family[fid] = row

    cols = ["gain_cost_{:g}".format(g) for g in gain_costs] + ["dollo"]
    df = pd.DataFrame.from_dict(per_family, orient="index", columns=cols)
    df.index.name = "family_id"
    df.to_csv(path, sep="\t")

    stable = (df.nunique(axis=1) == 1).sum()
    print("[sweep] {} families; {} ({:.1f}%) keep the same origin at every gain "
          "cost tested".format(len(df), stable, 100 * stable / max(len(df), 1)),
          file=stream)
    for col in cols:
        agree = int((df[col] == df["dollo"]).sum())
        print("          {:<18} agrees with Dollo for {} families ({:.1f}%)".format(
            col, agree, 100 * agree / max(len(df), 1)), file=stream)
    return df


def events_newick(tree, origins, losses, fmt="G{gains}_L{losses}",
                  keep_tip_names=True) -> str:
    """Single tree whose node labels carry the event counts.

    Internal labels are replaced outright, since the clade name is what you
    already know from the input tree. Tip labels keep the species name with the
    counts appended, because a tree of bare counts with no taxon names is
    unreadable -- pass keep_tip_names=False if you want them replaced too.

    `fmt` is a format string over {gains}, {losses} and {name}. Avoid Newick's
    reserved characters -- , : ; ( ) [ ] -- or the file will not parse
    downstream; "/" is safe in most viewers but not all.
    """
    out = copy.deepcopy(tree)
    # Only prepend the species name to tips when the template does not already
    # interpolate it, otherwise it appears twice.
    prepend = keep_tip_names and "{name}" not in fmt
    for node in out.traverse("preorder"):
        label = fmt.format(gains=origins.get(node.name, 0),
                           losses=losses.get(node.name, 0),
                           name=node.name)
        node.name = ("{}_{}".format(node.name, label)
                     if node.is_leaf() and prepend else label)
    return out.write(format=1)


def write_node_summary(path, tree, origins, losses):
    rows = [{
        "node": node.name,
        "is_tip": node.is_leaf(),
        "n_descendant_species": len(descendant_names(node)),
        "gains": origins.get(node.name, 0),
        "losses": losses.get(node.name, 0),
    } for node in tree.traverse("preorder")]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def write_counts_matrix(path, fams, species):
    species = sorted(species)
    with open(path, "w") as out:
        out.write("family_id\t" + "\t".join(species) + "\n")
        for fid, gbs in fams.items():
            out.write(fid + "\t" + "\t".join(
                str(len(gbs.get(sp, ()))) for sp in species) + "\n")


def reconcile_species(fams, leaves, exclude, prune, allow_mismatch):
    """Require the tree tips and the orthogroup species to be the same set.

    A mismatch is an error rather than a warning, because neither direction is
    harmless and both are silent:

      * a species in the tree but absent from every orthogroup is scored as
        absent in every family, so it contributes one loss to each -- inflating
        the loss count on its own branch by the size of the whole dataset;
      * a species in the orthogroups but not in the tree has all of its genes
        dropped, deflating gene counts and potentially moving origins tipwards.

    Resolve deliberately with --exclude (drop from both sides), --prune-tree
    (cut the tree down to the species that have data), or --allow-mismatch
    (proceed anyway, having understood the above).

    Returns the set of species to analyse.
    """
    excl = set(exclude)
    all_data = {sp for gbs in fams.values() for sp in gbs}
    all_tree = set(leaves)
    in_data, in_tree = all_data - excl, all_tree - excl

    unknown = excl - (all_data | all_tree)
    if unknown:
        print("[species] warning: --exclude names {} not found in either the tree "
              "or the orthogroups (typo?): {}".format(
                  len(unknown), ", ".join(sorted(unknown))), file=sys.stderr)

    only_data = sorted(in_data - in_tree)
    only_tree = sorted(in_tree - in_data)

    if only_data or only_tree:
        lines = ["species sets do not match between the tree and the "
                 "orthogroup file:"]
        if only_tree:
            lines.append("  {} in the tree with no genes in any orthogroup: {}".format(
                len(only_tree), ", ".join(only_tree)))
            lines.append("      -> each would be scored as a loss in every family")
        if only_data:
            lines.append("  {} in the orthogroups but not in the tree: {}".format(
                len(only_data), ", ".join(only_data)))
            lines.append("      -> all of their genes would be discarded")
        if not allow_mismatch and not (prune and not only_data):
            lines.append("")
            lines.append("Resolve with --exclude <sp ...>, --prune-tree, or "
                         "--allow-mismatch. Refusing to guess.")
            sys.exit("error: " + "\n".join(lines))
        print("[species] " + "\n           ".join(lines), file=sys.stderr)

    keep = in_data & in_tree
    if not keep:
        sys.exit("error: no species shared between the tree and the orthogroup file")
    return keep


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Gene family gain/loss and summary statistics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-g", "--orthogroups", required=True,
                    help="Broccoli / OrthoFinder orthogroup file")
    ap.add_argument("-t", "--tree", required=True,
                    help="species tree in Newick, with named internal nodes")
    ap.add_argument("-o", "--out", required=True, help="output family table (TSV)")
    ap.add_argument("--format", choices=["auto", "broccoli", "orthomcl", "orthofinder"],
                    default="auto")
    ap.add_argument("--sep", default="_",
                    help="separator between species tag and gene id")
    ap.add_argument("--species-position", choices=["prefix", "suffix"], default="prefix")
    ap.add_argument("--unassigned", metavar="TSV",
                    help="OrthoFinder Orthogroups_UnassignedGenes.tsv; adds these "
                         "genes as singleton families, flagged in the output as "
                         "from_unassigned so they can be filtered")
    ap.add_argument("--unassigned-prefix", default="UNASSIGNED",
                    help="family id prefix for --unassigned entries")
    ap.add_argument("--species-map", metavar="TSV",
                    help="two-column TSV renaming species from the orthogroup "
                         "file to the names used in the tree, e.g. "
                         "'Aiptasia_pallida<TAB>Apal'")
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="species to drop from both the tree and the orthogroups")
    ap.add_argument("--prune-tree", action="store_true",
                    help="cut the tree down to the species that have genes, "
                         "instead of erroring on a species-set mismatch")
    ap.add_argument("--allow-mismatch", action="store_true",
                    help="proceed despite a species-set mismatch (see --prune-tree "
                         "first; this option lets tree-only species be scored as "
                         "a loss in every family)")
    ap.add_argument("--model", choices=["dollo", "wagner"], default="dollo",
                    help="ancestral reconstruction model. dollo forces a single "
                         "origin at the MRCA; wagner allows independent gains, "
                         "priced by --gain-cost")
    ap.add_argument("--gain-cost", type=float, default=2.0,
                    help="cost of a gain relative to a loss (--model wagner). "
                         "1 is symmetric Wagner; large values converge on Dollo")
    ap.add_argument("--sweep", metavar="COSTS",
                    help="comma-separated gain costs, e.g. 1,1.5,2,3,5,10; writes "
                         "a per-family origin-stability table to --sweep-out")
    ap.add_argument("--sweep-out", default="origin_sweep.tsv",
                    help="output path for --sweep")
    ap.add_argument("--min-species", type=int, default=3,
                    help="minimum species for the expansion test")
    ap.add_argument("--min-genes", type=int, default=5,
                    help="minimum genes for the expansion test")
    ap.add_argument("--alpha", type=float, default=0.05, help="FDR threshold")
    ap.add_argument("--expansion-universe", choices=["present", "all"], default="present")
    ap.add_argument("--no-expansion", action="store_true")
    ap.add_argument("--node-summary", help="write per-node gains/losses TSV")
    ap.add_argument("--annotate-tree", metavar="PREFIX",
                    help="write PREFIX.origins.tre and PREFIX.losses.tre")
    ap.add_argument("--events-tree", metavar="PATH",
                    help="write a single tree whose node labels are the gain and "
                         "loss counts, for viewing in FigTree/iTOL")
    ap.add_argument("--label-format", default="G{gains}_L{losses}",
                    help="label template for --events-tree, over {gains}, "
                         "{losses} and {name}")
    ap.add_argument("--replace-tip-names", action="store_true",
                    help="use the label template for tips too, instead of "
                         "appending it to the species name")
    ap.add_argument("--counts-matrix", help="write a family x species count matrix")
    args = ap.parse_args(argv)

    tree, leaves, relabelled = load_tree(args.tree)
    if relabelled:
        print("[tree] warning: {} internal node(s) had no label and were named nN by "
              "preorder index; gains/losses there will be hard to interpret".format(
                  relabelled), file=sys.stderr)

    fmt = sniff_format(args.orthogroups) if args.format == "auto" else args.format
    fams = read_orthogroups(args.orthogroups, fmt, args.sep, args.species_position)
    if args.unassigned:
        extra = read_unassigned(args.unassigned, args.unassigned_prefix)
        clash = set(extra) & set(fams)
        if clash:
            sys.exit("error: {} family ids appear in both files, e.g. {}. Change "
                     "--unassigned-prefix.".format(
                         len(clash), ", ".join(sorted(clash)[:3])))
        fams.update(extra)
    if args.species_map:
        fams = apply_species_map(fams, read_species_map(args.species_map))
    print("[input] {} read as {}: {} families, {} tips in tree".format(
        args.orthogroups, fmt, len(fams), len(leaves)), file=sys.stderr)

    keep = reconcile_species(fams, leaves, args.exclude, args.prune_tree,
                             args.allow_mismatch)
    if set(leaves) != keep:
        tree, leaves = restrict_tree(tree, keep)

    fams = {fid: kept for fid, kept in
            ((f, {sp: g for sp, g in gbs.items() if sp in keep})
             for f, gbs in fams.items()) if kept}
    print("[species] analysing {} species, {} families".format(len(keep), len(fams)),
          file=sys.stderr)
    report_size_distribution(fams)

    df, origins, losses, totals = analyse(fams, tree, leaves, args)

    df.to_csv(args.out, sep="\t")
    print("[output] {} families -> {}".format(len(df), args.out), file=sys.stderr)

    # Dollo assigns exactly one origin per family; Wagner may assign several.
    total_gains, expected_gains = sum(origins.values()), int(df["n_gains"].sum())
    if total_gains != expected_gains:
        print("[output] warning: {} origins tallied per node vs {} summed over "
              "families -- these should be equal".format(
                  total_gains, expected_gains), file=sys.stderr)
    tips = {l.name for l in tree.get_leaves()}
    tip_gains = sum(v for k, v in origins.items() if k in tips)
    print("[output] gains: {} at tips, {} at internal nodes".format(
        tip_gains, sum(origins.values()) - tip_gains), file=sys.stderr)
    if args.unassigned:
        by_src = df.groupby("from_unassigned")["origin_is_tip"].agg(["sum", "count"])
        for src, row in by_src.iterrows():
            print("           {:<14} {} families, {} with a tip origin".format(
                "unassigned:" if src else "orthogroups:",
                int(row["count"]), int(row["sum"])), file=sys.stderr)

    multi = int((df["n_gains"] > 1).sum())
    if multi:
        print("[output] {} families ({:.1f}%) have more than one independent "
              "gain".format(multi, 100 * multi / max(len(df), 1)), file=sys.stderr)
    thin = int((df["origin_min_side"] == 1).sum())
    print("[output] {} families ({:.1f}%) have an origin resting on a single "
          "species on one side (origin_min_side == 1)".format(
              thin, 100 * thin / max(len(df), 1)), file=sys.stderr)

    if args.sweep:
        costs = [float(c) for c in args.sweep.split(",") if c.strip()]
        run_sweep(args.sweep_out, fams, tree, leaves, costs)
        print("[output] origin sweep -> {}".format(args.sweep_out), file=sys.stderr)

    if args.node_summary:
        write_node_summary(args.node_summary, tree, origins, losses)
        print("[output] node summary -> {}".format(args.node_summary), file=sys.stderr)
    if args.annotate_tree:
        for tag, table in (("origins", origins), ("losses", losses)):
            with open("{}.{}.tre".format(args.annotate_tree, tag), "w") as out:
                out.write(annotated_newick(tree, table) + "\n")
        print("[output] annotated trees -> {}.{{origins,losses}}.tre".format(
            args.annotate_tree), file=sys.stderr)
    if args.events_tree:
        with open(args.events_tree, "w") as out:
            out.write(events_newick(tree, origins, losses, args.label_format,
                                    not args.replace_tip_names) + "\n")
        print("[output] events tree -> {}".format(args.events_tree), file=sys.stderr)
    if args.counts_matrix:
        write_counts_matrix(args.counts_matrix, fams, keep)
        print("[output] count matrix -> {}".format(args.counts_matrix), file=sys.stderr)

    print("\n[summary] genes per species:", file=sys.stderr)
    for sp in sorted(totals, key=lambda s: -totals[s]):
        print("           {}\t{}".format(sp, totals[sp]), file=sys.stderr)


if __name__ == "__main__":
    main()