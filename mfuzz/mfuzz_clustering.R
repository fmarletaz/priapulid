## Mfuzz clustering of P. caudatus developmental expression, with TMM normalisation


library(tidyverse)
library(edgeR)
library(Mfuzz)

#setwd("~/Dropbox/Genomes/Priapulid/Analyses/mfuzz")

prefix     <- "Pcau_dev_m13"
dev_stages <- c("0d", "1d", "3d", "5d", "7d", "9d", "Larva")

## Replicate libraries -> stage (adult erythrocyte / nephridia samples left out)
stage_of <- c(Pc1_0d = "0d", Pc2_0d = "0d",
              Pc1_1d = "1d", Pc2_1d = "1d",
              Pc1_3d = "3d", Pc2_3d = "3d",
              Pc1_5d = "5d", Pc2_5d = "5d",
              Pc2_7d = "7d",                       # single replicate
              Pc1_9d = "9d", Pc2_9d = "9d",
              Larva1 = "Larva", PcauLarva2 = "Larva")

## ---- TMM normalisation across replicate libraries ---------------------------
tpm <- read_tsv("Pcau_stgtis_tpms.txt") %>%         
  mutate(Name = sub("\\.\\d+$", "", Name)) %>%         
  column_to_rownames("Name") %>% as.matrix()

tpm   <- tpm[, names(stage_of)]
stage <- unname(stage_of[colnames(tpm)])
tpm   <- tpm[rowSums(tpm) > 0, ]

nf      <- calcNormFactors(tpm, method = "TMM")
print(round(nf, 3))
tpm_tmm <- sweep(tpm, 2, nf, "/")

## Average replicates per stage; save for the tau analysis
expr <- sapply(dev_stages, function(s) rowMeans(tpm_tmm[, stage == s, drop = FALSE]))
write.csv(expr, "Pcau_stgtis_tpms_summ_dev_tmm_avg.txt")

## ---- Mfuzz ------------------------------------------------------------------
eset <- new("ExpressionSet", exprs = log2(expr + 1))
eset <- filter.std(eset, min.std = 0, visu = FALSE)     # drop genes with no variance
eset <- standardise(eset)                                # z-score per gene
cat("Genes clustered:", nrow(eset), "\n")

m <- mestimate(eset)
cat("m =", round(m, 3), "\n")

pdf(paste0(prefix, "_Dmin.pdf"), 4, 4)
Dmin(eset, m = m, crange = seq(4, 30, 2), repeats = 3, visu = TRUE)
dev.off()

set.seed(1)
cl <- mfuzz(eset, c = 13, m = m)
print(table(cl$cluster))

## ---- Outputs (read by mfuzz_enrichment.R) -----------------------------------
write.table(data.frame(cl.cluster = cl$cluster), paste0(prefix, "_clusters.tsv"),
            sep = "\t", quote = FALSE)
write.table(cl$centers, paste0(prefix, "_centers.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
saveRDS(list(eset = eset, cl = cl, m = m), paste0(prefix, ".rds"))

pdf(paste0(prefix, "_profiles.pdf"), 12, 9)
mfuzz.plot2(eset, cl, mfrow = c(4, 4), centre = TRUE, x11 = FALSE, time.labels = dev_stages)
dev.off()
par(mfrow = c(1, 1))
