# Datasets and code for the priapulid genome manuscript

[![DOI](https://zenodo.org/badge/1186296803.svg)](https://doi.org/10.5281/zenodo.22126220)

This repository contains the datasets and code used to analyse the priapulid genomes and to
generate the figures of the manuscript.

## Genomes

`genomes` holds the primary datasets.

*Priapulus caudatus*

- `PriCau_v1.fa.gz` — genome assembly after Hi-C scaffolding, with scaffolds renamed as chromosomes
- `PriCau_v1_ahfp.gtf.gz` — genome annotation
- `PriCau_v1.rep.bed` — location of annotated repeats
- `PriCau_v1_ahfp.cds` — CDS sequences (FASTA)
- `PriCau_v1_ahfp.pep` — peptides

*Halicryptus spinulosus*

- `Hali_chhp_rn.fasta.gz` — genome assembly
- `Hspi.aug.AGAT.rn.gtf.gz` — genome annotation, protein-coding genes (GTF)
- `Hspi_aug_agat.pep` — peptides
- `Hali_chhp_rn.rep.bed` — location of annotated repeats

## Synteny and ancestral linkage groups

`synteny` contains the files used for the dotplots and the ancestral linkage group
reconstruction. `prepSynt.py` generates the pairwise files (e.g. `Braflo-Caeele_synt.txt`)
from the reference BED files; the R scripts plot them and compute the enrichments
(e.g. `synt_ecdy.R`).

## Gene families

`gene_families` contains the main analysis notebooks. `gene_families.ipynb` examines the
phylogenetic distribution of gene loss and gene duplication, taking as input the families
built with Broccoli (`orthologous_groups_0623.txt`) and the reconciled trees from GeneRax
(`Priap_0124_nhx.lab.tre`). Its main outputs are `ortho0623_GainLoss.txt`, listing the gene
families with their inferred gains and losses, and `Priap_dup_alg.txt`, listing the ohnologue
candidates. `molecular_dating.ipynb` prepares the gene families used for dating.
`loss_rate2.py` converts a gain/loss table into loss rates along the branches of a species
tree, and `panther` holds the PANTHER family assignments per species used as an independent
check on the loss patterns.

Four subdirectories hold analyses added during revision.

- `broad_tax` — gene family reconstruction under an extended taxonomic sampling that adds
  three kinorhynch transcriptomes (`orthologous_groups_jun26.txt`, `generax_0707_crvm.nhx`,
  `orthoJun26_GainLoss_Crypto.txt`), used to test whether the duplication is restricted to
  priapulids (Supplementary Fig. 1). `extract_paralogues.py` pulls paralogue pairs and their
  duplication node out of the reconciliations, and `wgd_placement.py` scores each pair against
  the ohnologous chromosome pairs (`chrom_pairs.txt`, `Pricau.bed`) to test where on the tree
  the duplication sits.

- `alt_topologies` — gene gain and loss under both the Cryptovermes (`metazoa-panarthro.tre`)
  and Triradiata (`metazoa-triradiate.tre`) topologies, giving `ortho_gainloss_panartho.txt`
  and `ortho_gainloss_triradiata.txt`. `gene_family_stats.py` performs the Dollo and Wagner
  reconstructions; `loss_rate2.py` expresses losses as a rate over the repertoire available at
  the ancestral node of each branch and calls loss episodes by Tukey's upper fence (Fig. 4a,
  Extended Data Fig. 8). `compare_topologies.py` and `loss_convergence_matrix.py` compare the
  reconstructed gene content between topologies and between clades.

- `cafe_GO` — `cafe_GO.ipynb` runs the gene family expansion tests and the Gene Ontology
  enrichment of the families lost at each node, using the families available at the ancestral
  node as background (Fig. 4b). GO annotations come from eggNOG-mapper
  (`PriCau_v1_ahfp.eggnog.emapper.annotations`) formatted with `format_go_annotations.py`;
  `curated_GO_subset.txt` lists the terms shown in the figure.

- `hdf_rescue` — test for homology-detection failure. `hdf_rescue.py` builds a profile HMM for
  each family excluding the clade of interest, searches it against the target proteomes, and
  applies a bitscore threshold calibrated by leave-one-out on retained families. `rescore.py`
  reclassifies the resulting hits without repeating the search, according to whether a hit
  already belongs to another family and, if so, whether that family predates the node. Results
  are given for Nematoda and Tardigrada (`*_besthits_reclassified.tsv`, `nematoda_summary.txt`).

## Phylogenomics

`phylogenomics` contains the phylogenetic reconstruction based on the genomes analysed here
(Supplementary Fig. 2). `markers_filt_fa.tgz` and `markers_filt_tree.tgz` hold the filtered
single-marker alignments and their trees; `Concat-P0207.fa`/`.phy` is the concatenated
supermatrix with its partition file (`Concat-P0207.part.txt`), occupancy
(`Concat-P0207.compl.tsv`) and summary statistics (`Concat-P0207.stats.txt`).
`Concat-P0207.dh6.*` is the Dayhoff-6 recoded matrix and `Concat-P0207S.*` the reduced taxon
set. `P0207_C60.*` are the IQ-TREE outputs under the C60 profile mixture model, including the
ML tree (`.treefile`), the consensus tree (`.contree`) and the full log (`.iqtree`);
`P0207_LG4X.tre` is the tree obtained under LG4X. `phylobayes_runs.tgz` contains the
corresponding Bayesian analyses.

## Molecular dating

`dating_graham` contains the molecular dating of the priapulid duplication. Each
`*_ln_sy_sample.chronogram` file holds the posterior sample of node ages for one gene family,
obtained with PhyloBayes under a log-normal relaxed clock with fossil calibrations. The R
scripts summarise these across families to produce the age distributions and credibility
intervals shown in Fig. 2e and Supplementary Fig. 3.

## Homeobox genes

`homeodomain` contains the alignment of homeodomain sequences across protostomes
(`NewAllProtosHboxAln.fa`) and the resulting phylogeny, together with the reference homeobox
complements used for comparison (Extended Data Fig. 5).

## Expression clustering

`mfuzz` contains the code used to run fuzzy clustering (Mfuzz) on the staged expression data.
