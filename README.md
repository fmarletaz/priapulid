
# Datasets and code for the priapulid genome manuscript 
This repository contains datasets and code used to analyse the priapulid genomes and generate figures from the manuscript. 

* `gene_families` contains severak ipython notebooks, notably `gene_family.ipynb` used to examine the phylogenetic distribution of gene loss and gene duplication taking families generated using broccoli : `orthologous_groups_0623.txt` and reconciliated trees obtained with generax: `Priap_0124_nhx.lab.tre` as input. Relevant output files are : `ortho0623_GainLoss.txt` listing gene families and `Priap_dup_alg.txt` listing ohnologue candidates. 
* `synteny` contains files used for dotplot and ancestral linkage group reconstruction, particularly the `prepSynt.py` script used to generate pairwise files (e.g. `Braflo-Caeele_synt.txt`) together with reference bed files. Some R files contain code to plot and computate enrichments (e.g. `synt_ecdy.R`). 
* `mfuzz` contains code to run the fuzzy clustering mfuzz on staged data. 
* `genomes` contains important primary datasets including in particular:
	For priapulus
	- `PriCau_v1.fa.gz`: the P. caudatus genome assembly after HiC scaffolding and With scaffolds renamed into chromsomes
	- `PriCau_v1_ahfp.gtf.gz`: the annotatin of P. caudatus genome 
	- `PriCau_v1.rep.bed`: location of annotated repeats
	- `PriCau_v1_ahfp.cds`: cds sequence (fasta)
	- `PriCau_v1_ahfp.pep`: peptives
	For Halicryptus: 
	- `Hali_chhp_rn.fasta`: H. spinulosus genome assembly
	- `Hspi.aug.AGAT.rn.gtf`:  H. spinulosus genome annotation in GTF format (protein coding genes)
	- `Hspi_aug_agat.pep`: peptides 
	- `Hali_chhp_rn.rep.bed`: location of annotated repeats


/Live-Manuscripts/Priapulid/zenodo/genomes