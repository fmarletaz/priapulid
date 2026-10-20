## Mfuzz cluster enrichment: WGD ohnologues and phylostrata
## Priapulus caudatus developmental expression (Nature revision)

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)

#setwd("~/Dropbox/Genomes/Priapulid/Analyses/mfuzz")

## ---- Parameters -------------------------------------------------------------
prefix   <- "Pcau_dev_m13"                              # swap for the dev-only run
cl_order <- c(4,12,1,7,13,9,6,2,10,3,8,11) # display order (update for new run)

pal_expr <- colorRampPalette(brewer.pal(8, "YlOrRd"))(100)
pal_p    <- colorRampPalette(brewer.pal(8, "PuBu"))(100)

## ---- Inputs -----------------------------------------------------------------
gcl     <- read.table(paste0(prefix, "_clusters.tsv"), header = TRUE)
centers <- read.table(paste0(prefix, "_centers.tsv"),  header = TRUE)  # rows = clusters
ohno    <- read.table("../Gene_families/Pc_dup_raw.txt")$V1
strata  <- read.table("Pcau_phylostrata.txt", col.names = c("gene", "stratum", "family"))

clusters <- tibble(gene = rownames(gcl), cluster = as.character(gcl$cl.cluster))
universe <- clusters$gene

## ---- Hypergeometric enrichment, one gene set vs every cluster ---------------
## Returns cluster size, observed, expected and P for every cluster (incl. zeros)
cluster_enrichment <- function(set, set_name) {
  set <- intersect(set, universe)
  N <- length(universe)
  n <- length(set)
  clusters %>%
    group_by(cluster) %>%
    summarise(cluster_size = n(), observed = sum(gene %in% set), .groups = "drop") %>%
    mutate(set      = set_name,
           set_size = n,
           expected = n * cluster_size / N,
           log2_OE  = log2(observed / expected),
           pval     = phyper(observed - 1, cluster_size, N - cluster_size, n,
                             lower.tail = FALSE))
}

gene_sets <- c(list(ohnologues = ohno),
               split(strata$gene, strata$stratum),
               list(Priapulus = setdiff(universe, strata$gene)))

enrich <- imap_dfr(gene_sets, cluster_enrichment) %>%
  #mutate(padj = p.adjust(pval, method = "BH")) %>%   # BH across all tests
  group_by(set) %>% mutate(padj = p.adjust(pval, "BH")) %>% ungroup() %>%
  relocate(set, set_size, cluster, cluster_size, observed, expected)

cp <- imap_dfr(gene_sets, function(g, nm) {
  enricher(g, TERM2GENE = select(clusters, cluster, gene), universe = universe,
           minGSSize = 1, maxGSSize = 1e5, pvalueCutoff = 1, qvalueCutoff = 1)@result %>%
    transmute(set = nm, cluster = ID, cp_pval = pvalue, cp_padj = p.adjust)
})

cmp <- enrich %>% left_join(cp, by = c("set", "cluster"))
all.equal(cmp$pval[!is.na(cmp$cp_pval)], cmp$cp_pval[!is.na(cmp$cp_pval)])  # should be TRUE
plot(-log10(cmp$padj), -log10(cmp$cp_padj), xlim = c(0, 5), ylim = c(0, 5))
abline(0, 1); abline(h = -log10(0.05), v = -log10(0.05), lty = 2)

## Supplementary table requested by Referee 5
write_tsv(enrich, paste0(prefix, "_cluster_enrichment.tsv"))

## Numbers to report in the text / legend
cat("Genes clustered:", length(universe), "\n")
print(count(clusters, cluster))

## ---- Figures ----------------------------------------------------------------
## Cluster centres
ctr <- as.matrix(centers)
rownames(ctr) <- seq_len(nrow(ctr))
pheatmap(ctr[as.character(cl_order), ],
         cluster_rows = FALSE, cluster_cols = FALSE, color = pal_expr)

## Ohnologue enrichment per cluster (adjusted P)
ohno_p <- enrich %>%
  filter(set == "ohnologues") %>%
  select(cluster, padj) %>%
  column_to_rownames("cluster")
pheatmap(-log10(ohno_p[as.character(cl_order), , drop = FALSE]),
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = seq(0, 10, length.out = 101), color = pal_p)

## Phylostrata x cluster enrichment (adjusted P)
str_order=c('Metazoa','Bilateria','Protostomia','Ecdysozoa','Priapulida','Priapulus','ohnologues')
strata_p <- enrich %>%
  filter(set != "ohnologues") %>%
  select(set, cluster, padj) %>%
  pivot_wider(names_from = cluster, values_from = padj) %>%
  column_to_rownames("set")
pheatmap(-log10(t(strata_p[str_order, as.character(cl_order)])),   # reorder rows by age if needed
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = seq(0, 10, length.out = 101), color = pal_p)

strata_p <- enrich %>%
  select(set, cluster, padj) %>%
  pivot_wider(names_from = cluster, values_from = padj) %>%
  column_to_rownames("set")
pheatmap(-log10(t(strata_p[str_order, as.character(cl_order)])),
         cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_col = length(str_order) - 1,              # separates ohnologue row
         breaks = seq(0, 10, length.out = 101), color = pal_p)

## Cluster profiles coloured by ohnologue enrichment
centers_long <- as_tibble(ctr, rownames = "cluster") %>%
  pivot_longer(-cluster, names_to = "stage", values_to = "center") %>%
  mutate(stage = factor(stage, levels = colnames(ctr))) %>%
  left_join(rownames_to_column(ohno_p, "cluster"), by = "cluster")

ggplot(centers_long, aes(stage, center, group = cluster, colour = -log10(padj))) +
  geom_line() +
  theme_classic()

## ---- Expression divergence of ohnologues ------------------------------------
## Duplication-event table: one row per event, P. caudatus genes in `pcg` as
## "sideA_genes/sideB_genes"; ALG = NA when not on an ancestral chromosome pair
events <- read_tsv("../Gene_families/Pca_dup_alg2.txt") %>%   # adapt file name
  mutate(event = paste(FID, Age, row_number(), sep = "_"))

print(count(events, Age, on_ALG_pair = !is.na(ALG)))         # compare with 418 / 112

dup <- events %>%
  filter(!is.na(ALG)) %>%                                    # ancestral chromosome-pair filter
  select(event, FID, Age, pcg) %>%
  separate_wider_delim(pcg, "/", names = c("A", "B")) %>%
  pivot_longer(c(A, B), names_to = "side", values_to = "gene") %>%
  separate_longer_delim(gene, ",") %>%
  filter(!is.na(gene), gene != "") %>%
  group_by(event) %>% filter(n_distinct(side) == 2) %>% ungroup()

cat("Events:", n_distinct(dup$event), " Genes:", n_distinct(dup$gene), "(expect ~1679)\n")
write_tsv(dup, "../Gene_families/Pc_wgd_events.tsv")

## Keep events with clustered genes on both sides
dup_cl <- dup %>%
  inner_join(clusters, by = "gene") %>%
  group_by(event) %>% filter(n_distinct(side) == 2) %>% ungroup()

## Asymmetric = the two daughter lineages share no expression cluster
## (later lineage-specific duplicates within one side are not counted as divergence)
frac_asym <- function(d) {
  d %>% group_by(event) %>%
    summarise(asym = length(intersect(cluster[side == "A"], cluster[side == "B"])) == 0,
              .groups = "drop") %>%
    pull(asym) %>% mean()
}

## 1. Cluster-based: permute labels across all clustered genes
##    (keeps cluster sizes and event sizes fixed)
obs_asym <- frac_asym(dup_cl)

set.seed(1)
perm_asym <- replicate(1000, {
  shuffled <- setNames(sample(clusters$cluster), clusters$gene)
  dup_cl %>% mutate(cluster = shuffled[gene]) %>% frac_asym()
})

cat("Events tested:", n_distinct(dup_cl$event), "\n",
    "Observed asymmetric:", round(obs_asym, 3), "\n",
    "Null mean [95%]:", round(mean(perm_asym), 3),
    round(quantile(perm_asym, c(0.025, 0.975)), 3), "\n",
    "P(null <= obs):", (sum(perm_asym <= obs_asym) + 1) / 1001, "\n",
    "P(null >= obs):", (sum(perm_asym >= obs_asym) + 1) / 1001, "\n")

## 2. Cluster-free: correlation of expression profiles, paralogues vs random pairs
expr_z <- exprs(eset)                                   # standardised dev-only matrix

## Ohnologue pairs = every side-A gene x every side-B gene of the same event
para_pairs <- inner_join(filter(dup_cl, side == "A"), filter(dup_cl, side == "B"),
                         by = "event", suffix = c("1", "2"),
                         relationship = "many-to-many") %>%
  filter(gene1 %in% rownames(expr_z), gene2 %in% rownames(expr_z)) %>%
  transmute(event, g1 = gene1, g2 = gene2)

pair_cor <- function(g1, g2) map2_dbl(g1, g2, ~ cor(expr_z[.x, ], expr_z[.y, ]))

para_pairs$r <- pair_cor(para_pairs$g1, para_pairs$g2)

set.seed(2)
rand <- tibble(g1 = sample(rownames(expr_z), 10000, replace = TRUE),
               g2 = sample(rownames(expr_z), 10000, replace = TRUE)) %>%
  filter(g1 != g2)
rand$r <- pair_cor(rand$g1, rand$g2)

wilcox.test(para_pairs$r, rand$r)

## ---- Expression specificity (tau) -------------------------------------------
## tau (Yanai et al. 2005) on log2(TPM + 1) across the developmental stages
## (not the standardised matrix: tau needs absolute levels)
tpms <- read.csv(file = 'Pcau_summed_averaged_tpms.csv', row.names = 1, header = TRUE, sep = "\t")
tpms <- tpms[,c(0:7)]
#dev_stages <- colnames(exprs(eset))                       # same 7 stages as the clustering
tpm <- tpms %>% as.matrix()

tau_of <- function(x) {
  x <- log2(x + 1)
  if (max(x) == 0) return(NA_real_)
  sum(1 - x / max(x)) / (length(x) - 1)
}
gene_tau <- tibble(gene = rownames(tpm), tau = apply(tpm, 1, tau_of))

## (a) Are ohnologues more stage-specific than other genes?
gene_tau %>%
  filter(gene %in% universe) %>%
  mutate(class = if_else(gene %in% dup$gene, "ohnologue", "other")) %>%
  { print(group_by(., class) %>% summarise(n = n(), median_tau = median(tau, na.rm = TRUE)));
    print(wilcox.test(tau ~ class, data = .)) }

## corrected by checking only genes older than duplication as contrast
ohno_strata=c('Metazoa','Bilateria','Protostomia','Ecdysozoa','Priapulida')
old <- strata$gene[strata$stratum %in% ohno_strata]   # strata older than the WGD
gene_tau %>% filter(gene %in% universe, gene %in% c(old, dup$gene)) %>%
  mutate(class = if_else(gene %in% dup$gene, "ohnologue", "other")) %>%
  { print(group_by(., class) %>% summarise(median(tau, na.rm = TRUE))); wilcox.test(tau ~ class, data = .) }


## (b) tau-matched null for the asymmetry test: shuffle cluster labels only
##     among genes in the same tau decile (controls for broad vs specific expression)
clusters_tau <- clusters %>%
  left_join(gene_tau, by = "gene") %>%
  mutate(tau_bin = ntile(tau, 10))

set.seed(3)
perm_asym_tau <- replicate(1000, {
  shuffled <- clusters_tau %>%
    group_by(tau_bin) %>% mutate(cluster = sample(cluster)) %>% ungroup() %>%
    { setNames(.$cluster, .$gene) }
  dup_cl %>% mutate(cluster = shuffled[gene]) %>% frac_asym()
})
cat("tau-matched null mean [95%]:", round(mean(perm_asym_tau), 3),
    round(quantile(perm_asym_tau, c(0.025, 0.975)), 3),
    " P(null <= obs):", (sum(perm_asym_tau <= obs_asym) + 1) / 1001, "\n")

## (c) Complementarity: does combining the two sides broaden expression?
##     drop = mean(tau_A, tau_B) - tau(A + B); larger = more complementary
side_sum <- function(genes) colSums(tpm[intersect(genes, rownames(tpm)), , drop = FALSE])

compl <- dup %>%
  filter(gene %in% rownames(tpm)) %>%
  group_by(event) %>% filter(n_distinct(side) == 2) %>%
  summarise(tau_A  = tau_of(side_sum(gene[side == "A"])),
            tau_B  = tau_of(side_sum(gene[side == "B"])),
            tau_AB = tau_of(side_sum(gene)), .groups = "drop") %>%
  mutate(drop = (tau_A + tau_B) / 2 - tau_AB)

## Null: random gene pairs, same calculation
set.seed(4)
rand_drop <- replicate(5000, {
  g <- sample(rownames(tpm), 2)
  (tau_of(tpm[g[1], ]) + tau_of(tpm[g[2], ])) / 2 - tau_of(colSums(tpm[g, ]))
})
cat("Complementarity (tau drop): ohnologues median", round(median(compl$drop, na.rm = TRUE), 3),
    " random median", round(median(rand_drop, na.rm = TRUE), 3), "\n")
print(wilcox.test(compl$drop, rand_drop))
bind_rows(mutate(para_pairs, set = "ohnologue pairs"),
          mutate(rand, set = "random pairs")) %>%
  ggplot(aes(r, colour = set)) + geom_density() + theme_classic() +
  labs(x = "Pearson r between expression profiles")

asym_tau <- dup %>%
  left_join(gene_tau, by = "gene") %>%
  group_by(event, side) %>% summarise(tau = mean(tau, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = side, values_from = tau) %>%
  mutate(d = abs(A - B))

pool <- intersect(old, gene_tau$gene[!is.na(gene_tau$tau)])
set.seed(5)
rand_d <- replicate(5000, { g <- sample(pool, 2); abs(diff(gene_tau$tau[match(g, gene_tau$gene)])) })

median(asym_tau$d, na.rm = TRUE); median(rand_d)
wilcox.test(asym_tau$d, rand_d)



## all genes from any WGD-age duplication event, filtered or not
all_dup_genes <- events %>%
  pull(pcg) %>% str_split("[/,]") %>% unlist() %>% discard(~ .x == "")

## single-copy: P. caudatus has exactly one gene in the family, same age range as ohnologues
single_copy <- strata %>%
  filter(stratum %in% ohno_strata) %>%
  add_count(family, name = "n_pcau") %>%
  filter(n_pcau == 1, !gene %in% all_dup_genes) %>%
  pull(gene)
length(single_copy)

gene_tau %>%
  filter(gene %in% universe, gene %in% c(single_copy, dup$gene)) %>%
  mutate(class = if_else(gene %in% dup$gene, "ohnologue", "single-copy")) %>%
  { print(group_by(., class) %>% summarise(n = n(), median_tau = median(tau, na.rm = TRUE)));
    print(wilcox.test(tau ~ class, data = .)) }



lvl <- tibble(gene = rownames(tpm), mean_expr = log2(rowMeans(tpm) + 1))

d_sc <- gene_tau %>%
  filter(gene %in% universe, gene %in% c(single_copy, dup$gene)) %>%
  mutate(class = if_else(gene %in% dup$gene, "ohnologue", "single-copy")) %>%
  left_join(lvl, by = "gene")

summary(lm(tau ~ class + mean_expr, data = d_sc))      # class effect after controlling for level
d_sc %>% mutate(bin = ntile(mean_expr, 5)) %>%
  group_by(bin, class) %>% summarise(median_tau = median(tau, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = class, values_from = median_tau)

library(ggplot2)

p1 <- ggplot(d_sc, aes(class, tau, fill = class)) +
  geom_violin(colour = NA, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(x = NULL, y = "Expression specificity (τ)") +
  theme_classic() + theme(legend.position = "none")

p2 <- d_sc %>% mutate(quintile = ntile(mean_expr, 5)) %>%
  group_by(quintile, class) %>% summarise(tau = median(tau, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(quintile, tau, colour = class)) +
  geom_line() + geom_point() +
  labs(x = "Expression level (quintile)", y = "Median τ", colour = NULL) +
  theme_classic()

ggsave("Fig2d_tau.pdf", p1, width = 2.5, height = 3)
ggsave("Fig2d_tau_quintiles.pdf", p2, width = 3, height = 3)
