library('optparse')
library('Mfuzz')
library('RColorBrewer')
library('cluster')
library('factoextra')


setwd('/Users/matar/Desktop/ATAC_and_RNA_Analysis/Lpict/RNAseq/amorannot/')
#option_list = list(
#  make_option(c("-i", "--infile"), type="character", , 
#             help="Input expression matrix"),
# make_option(c("-k", "--nbClusters"), type="numeric", , 
#              help="Number of clusters"),
#  make_option(c("-o", "--outDir"), type="character", default="mfuzz_clusters/", 
#              help="Output directory [default= %default]")
#);

#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);

df <- read.csv(file = 'Lpic_tmm_rnaseg_aamroannot.csv', row.names = 1, header = TRUE, sep = "\t")
df <- df[apply(df, 1, var) != 0, ]
colnames(df) <- gsub("^x", "", colnames(df))


mat <- as.matrix(df)
eset <- new("ExpressionSet", exprs=mat)
eset <- filter.NA(eset, thres=0.01)
eset <- standardise(eset)
colnames(eset) <- gsub("^x", "", colnames(eset))

m1 <- mestimate(eset)


dmin = Dmin(eset, m1, crange=seq(4, 25, 1), repeats=5, visu=TRUE)



cl <- mfuzz(eset, c=16, m=m1)

#plots
pdf(paste( file = "Lpic_amroannot_run1_16.pdf", sep=''))
mfuzz.plot2(eset, cl=cl, mfrow=c(3, 3), x11 = FALSE)
dev.off()

#clusters
outfile <- paste(file = 'Lpic_amroannot_clusters_run1_16.tsv', sep='')
write.table(as.data.frame(cl$cluster), file=outfile, quote=FALSE, sep='\t', col.names = NA)


#centers
outfile <- paste(file = 'Lpic_amroannot_centres_run1_16.pdf.tsv', sep='')
write.table(as.data.frame(cl$centers), file=outfile, quote=FALSE, sep='\t', col.names = NA)

#r.object
save(cl, file = paste(file = 'Lpic_amroannot_run1_16.pdf.RData', sep=''))