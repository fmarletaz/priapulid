require(tidyverse)
require(pheatmap)
require(clusterProfiler)
#install.packages("clusterProfiler")
setwd("~/Dropbox/Genomes/Priapulid/Analyses/mfuzz")


gcl=read.table("priap_clusters_run1_14.tsv",h=T)
head(gcl)
pdr=read.table("../Gene_families/Pc_dup_raw.txt")
dl=as.vector(pdr$V1)


?enricher
gtc=data.frame(term=gcl$cl.cluster,gene=row.names(gcl))


enpd<-enricher(gene = dl,
  TERM2GENE = gtc,  # 2-column df: cluster_id, gene_id
  universe = as.vector(gtc$gene),
  maxGSSize=5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
gtc %>% filter(term=='5') %>% nrow()
enpdd %>% arrange(as.numeric(ID)) %>% select(ID,p.adjust) -> pvhde

as.data.frame(enpd@result) -> enpdd
clct=read.table("priap_centres_run1_14.pdf.tsv",h=T)
clct %>% rowid_to_column(var='cluster') %>% pivot_longer(X0d:NephBR,names_to="stage",values_to="center") -> clctl
ggplot(clctl,aes(x=stage,y=center,group=cluster,color=-log10(pval)))+geom_line()

length(c(9,13,6,10,11,3,4,12,14,1,2,4,7,8))
pheatmap(clct[c(9,13,8,6,10,7,11,3,1,4,12,14,2),c(1:9)],cluster_rows = F,cluster_cols = F,col=coos)
pheatmap(-log10(clct[c(9,13,8,6,10,7,11,3,1,4,12,14,2),10]),cluster_rows = F,cluster_cols = F,breaks=seq(0,10,length.out=100),col=cols)
clctl$pval=as.vector(pvhde[,"p.adjust"])
require(RColorBrewer)
cols=colorRampPalette(brewer.pal(8,'PuBu'))(100)
coos=colorRampPalette(brewer.pal(8,'YlOrRd'))(100)

gag=read.table("Pcau_phylostrata.txt")
names(gag)=c('gid','strata','fid')

for l in levels(as.factor(gag$strata)) {
  print(l)
}
  print(l)

enpd<-enricher(gene = dl,
               TERM2GENE = gtc,  # 2-column df: cluster_id, gene_id
               universe = as.vector(gtc$gene),
               maxGSSize=5000,
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2
)

mfzs=read.table("../Gene_families/mfuzz_strat.tsv",h=T)

mfzs %>% select(strat,clust,pval) %>% pivot_wider(names_from=clust,values_from=pval) %>% column_to_rownames(var="strat") -> mfsz.m
mfsz.mr<-mfsz.m[c(2,5,7,6,4,3,1),c('9','13','8','6','10','7','11','3','1','4','12','14','2')]
pheatmap(-log10(mfsz.mr),breaks=seq(0,15,length.out=100),col=cols,cluster_rows = F, cluster_cols = F)

mfzs$padj<-p.adjust(mfzs$pval, "BH")

## Check subfonctionalisation 

tpms<-read_tsv("Pcau_summed_averaged_tpms.csv")

