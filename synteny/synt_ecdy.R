require(RColorBrewer)
require(pheatmap)
require(tidyverse)


source("~/Dropbox/Genomes/Synteny/synt-forge/synteny-functions.R")
setwd("~/Dropbox/Genomes/Synteny/synt-forge/Priapulid/")


bflo_col=c("#1F78B4","#1B9E77","#33A02C","#B2DF8A","#FB9A99","#E31A1C","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00","#FDBF6F","#BEBADA","#ffd92f")
palette(bflo_col)
par(mar=c(4,3,3,3))


# goal is to plot a matrix of comparisons against amphioxus

## priapulus [1]
pcau_bflo_synt=read.table("Pcau_Bflo_synt.txt",h=T,stringsAsFactors=T)
pcau_bflo_synt.exp<-testEnrich(pcau_bflo_synt)
plotSyntCLGrs(ordSynt(pcau_bflo_synt.exp),'Amphioxus','Priapulus')
plotSyntCLGrs((pcau_bflo_synt.exp),'Amphioxus','Priapulus')

legend('bottomright',levels(pcau_bflo_synt.exp$clg),fill=bflo_col)

## Halicryptus [2]
hspi_bflo_synt=read.table("Hspi_Bflo_syntc.txt",h=T,stringsAsFactors=T)
hspi_bflo_synt.exp<-testEnrich(hspi_bflo_synt)
plotSyntCLGrs(ordSynt(hspi_bflo_synt.exp),'Amphioxus','Halicryptus')
legend('bottomright',levels(hspi_bflo_synt.exp$clg),fill=bflo_col)
plotSyntCLGrs((hspi_bflo_synt.exp),'Amphioxus','Halicryptus')

## Gordius [3]
bflo_gspe_synt=read.table("Braflo-Gorspe_synt.txt",h=T,stringsAsFactors=T)
bflo_gspe_synt.exp<-testEnrich(bflo_gspe_synt)
plotSyntCLGrsi(ordSynt(bflo_gspe_synt.exp),'Amphioxus','Gordius')
legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)
plotSyntCLGrsi((bflo_gspe_synt.exp),'Amphioxus','Gordius')

## Caenorhabditis [4]
bflo_cele_synt=read.table("Braflo-Caeele_synt.txt",h=T,stringsAsFactors=T)
bflo_cele_synt.exp<-testEnrich(bflo_cele_synt)
plotSyntCLGrsi((bflo_cele_synt.exp),'Amphioxus','Caenorhabditis')
legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Meloidogyne [5]
bflo_mhap_synt=read.table("Braflo-Melhap_synt.txt",h=T,stringsAsFactors=T)
bflo_mhap_synt.exp<-testEnrich(bflo_mhap_synt)
plotSyntCLGrsi((bflo_mhap_synt.exp),'Amphioxus','Meloidogyne')
legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Ramazzottius [6]
bflo_rvar_synt=read.table("Braflo-Ramvar_synt.txt",h=T,stringsAsFactors=T)
bflo_rvar_synt.exp<-testEnrich(bflo_rvar_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrsi((bflo_rvar_synt.exp),'Amphioxus','Ramazzottius')
legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Pycnogonum [7]
bflo_plit_synt=read.table("Braflo-Pyclit_synt.txt",h=T,stringsAsFactors=T)
bflo_plit_synt.exp<-testEnrich(bflo_plit_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrsi((bflo_plit_synt.exp),'Amphioxus','Pycnogonum')
legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Trichonephila Tant-Bflo_synt.txt [8]
tant_bflo_synt=read.table("Tant-Bflo_synt.txt",h=T,stringsAsFactors=T)
tant_bflo_synt.exp<-testEnrich(tant_bflo_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrs((tant_bflo_synt.exp),'Amphioxus','Trichonephila')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Carcinoscorpius  horseshoe crab [9]
bflo_crot_synt=read.table("Braflo-Carrot_synt.txt",h=T,stringsAsFactors=T)
bflo_crot_synt.exp<-testEnrich(bflo_crot_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrsi((bflo_crot_synt.exp),'Amphioxus','Carcinoscorpius')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Myriapods Strigamia [10]
sacu_bflo_synt=read.table("Sacu_Bflo_synt.txt",h=T,stringsAsFactors=T)
sacu_bflo_synt.exp<-testEnrich(sacu_bflo_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrs((sacu_bflo_synt.exp),'Amphioxus','Strigamia')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Myriapods Trigoniulus [11]
tcor_bflo_synt=read.table("Tcor_Bflo_synt.txt",h=T,stringsAsFactors=T)
tcor_bflo_synt.exp<-testEnrich(tcor_bflo_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrs((tcor_bflo_synt.exp),'Amphioxus','Trigoniulus')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)


## Crab Scylla crab [12]
spar_bflo_synt=read.table("Spar-Bflo_synt.txt",h=T,stringsAsFactors=T)
spar_bflo_synt.exp<-testEnrich(spar_bflo_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrs((spar_bflo_synt.exp),'Amphioxus','Scylla')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Tribolium [13]
bflo_tcas_synt=read.table("Braflo-Tricas_synt.txt",h=T,stringsAsFactors=T)
bflo_tcas_synt.exp<-testEnrich(bflo_tcas_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrsi((bflo_tcas_synt.exp),'Amphioxus','Tribolium')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Drosophila [14]
bflo_dmel_synt=read.table("Dmel-Bflo_synt.txt",h=T,stringsAsFactors=T)
bflo_dmel_synt.exp<-testEnrich(bflo_dmel_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrsi((bflo_dmel_synt.exp),'Amphioxus','Drosophila')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)

## Idodex [15]
bflo_ipac_synt=read.table("Braflo-Ixopac_synt.txt",h=T,stringsAsFactors=T)
bflo_ipac_synt.exp<-testEnrich(bflo_ipac_synt)
#bflo_rvar_synt %>% group_by(s2chr) %>% summarise(n=n()) %>% arrange(-n) %>% print(n=25)
plotSyntCLGrsi((bflo_ipac_synt.exp),'Amphioxus','Ixodes')
#legend('bottomright',levels(bflo_gspe_synt.exp$clg),fill=bflo_col)
par(mfrow=c(5,3))
par(mar=c(3,3,3,2))
par(mgp=c(1.75, 0.75, 0))

plotSyntCLGrs((pcau_bflo_synt.exp),'Amphioxus','Priapulus')
plotSyntCLGrs((hspi_bflo_synt.exp),'Amphioxus','Halicryptus')
plotSyntCLGrsi((bflo_gspe_synt.exp),'Amphioxus','Gordius')‡
plotSyntCLGrsi((bflo_cele_synt.exp),'Amphioxus','Caenorhabditis')
plotSyntCLGrsi((bflo_mhap_synt.exp),'Amphioxus','Meloidogyne')
plotSyntCLGrsi((bflo_rvar_synt.exp),'Amphioxus','Ramazzottius')
plotSyntCLGrsi((bflo_ipac_synt.exp),'Amphioxus','Ixodes')
plotSyntCLGrs((tant_bflo_synt.exp),'Amphioxus','Trichonephila')
plotSyntCLGrsi((bflo_crot_synt.exp),'Amphioxus','Carcinoscorpius')
plotSyntCLGrs((sacu_bflo_synt.exp),'Amphioxus','Strigamia')
plotSyntCLGrs((tcor_bflo_synt.exp),'Amphioxus','Trigoniulus')
plotSyntCLGrs((spar_bflo_synt.exp),'Amphioxus','Scylla')
plotSyntCLGrsi((bflo_tcas_synt.exp),'Amphioxus','Tribolium')
plotSyntCLGrs((bflo_dmel_synt.exp),'Amphioxus','Drosophila')

plotSyntCLGrs((pcau_bflo_synt.exp),'Priapulus','Priapulus')
plotSyntCLGrs((hspi_bflo_synt.exp),'Halicryptus','Halicryptus')
plotSyntCLGrsi((bflo_gspe_synt.exp),'Gordius','Gordius')
plotSyntCLGrsi((bflo_cele_synt.exp),'Caenorhabditis','Caenorhabditis')
plotSyntCLGrsi((bflo_mhap_synt.exp),'Meloidogyne','Meloidogyne')
plotSyntCLGrsi((bflo_rvar_synt.exp),'Ramazzottius','Ramazzottius')
plotSyntCLGrsi((bflo_plit_synt.exp),'Pycnogonum','Pycnogonum')
plotSyntCLGrs((tant_bflo_synt.exp),'Trichonephila','Trichonephila')
plotSyntCLGrsi((bflo_crot_synt.exp),'Carcinoscorpius','Carcinoscorpius')
plotSyntCLGrs((sacu_bflo_synt.exp),'Strigamia','Strigamia')
plotSyntCLGrs((tcor_bflo_synt.exp),'Trigoniulus','Trigoniulus')
plotSyntCLGrs((spar_bflo_synt.exp),'Scylla','Scylla')
plotSyntCLGrsi((bflo_tcas_synt.exp),'Tribolium','Tribolium')
plotSyntCLGrsi((bflo_dmel_synt.exp),'Amphioxus','Drosophila')



ordSyntCl(pcau_bflo_synt.exp) -> pcau_bflo_synt.ord
ordSyntCl(hspi_bflo_synt.exp) -> hspi_bflo_synt.ord
ordSyntRw(bflo_gspe_synt.exp) -> bflo_gspe_synt.ord
ordSyntRw(bflo_cele_synt.exp) -> bflo_cele_synt.ord
ordSyntRw(bflo_mhap_synt.exp) -> bflo_mhap_synt.ord
ordSyntRw(bflo_rvar_synt.exp) -> bflo_rvar_synt.ord
ordSyntRw(bflo_plit_synt.exp) -> bflo_plit_synt.ord
ordSyntCl(tant_bflo_synt.exp) -> tant_bflo_synt.ord
ordSyntRw(bflo_crot_synt.exp) -> bflo_crot_synt.ord
ordSyntRw(bflo_ipac_synt.exp) -> bflo_ipac_synt.ord
ordSyntCl(sacu_bflo_synt.exp) -> sacu_bflo_synt.ord
ordSyntCl(tcor_bflo_synt.exp) -> tcor_bflo_synt.ord
ordSyntCl(spar_bflo_synt.exp) -> spar_bflo_synt.ord
ordSyntRw(bflo_tcas_synt.exp) -> bflo_tcas_synt.ord
ordSyntRw(bflo_dmel_synt.exp) -> bflo_mel_synt.ord


pcau_bflo_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Priapulus') -> pcau_bflo_synt.mpa
hspi_bflo_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Halicryptus') -> hspi_bflo_synt.mpa
bflo_gspe_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Gordius') -> bflo_gspe_synt.mpa
bflo_cele_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Caenorhabditis') -> bflo_cele_synt.mpa
bflo_mhap_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Meloidogyne') -> bflo_mhap_synt.mpa
bflo_rvar_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Ramazzottius') -> bflo_rvar_synt.mpa
bflo_plit_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Pycnogonum') -> bflo_plit_synt.mpa
tant_bflo_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Trichonephila') -> tant_bflo_synt.mpa
bflo_crot_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Carcinoscorpius') -> bflo_crot_synt.mpa
bflo_ipac_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Ixodes') -> bflo_ipac_synt.mpa
sacu_bflo_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Strigamia') -> sacu_bflo_synt.mpa
tcor_bflo_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Trigoniulus') -> tcor_bflo_synt.mpa
spar_bflo_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Scylla') -> spar_bflo_synt.mpa
bflo_tcas_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Tribolium') -> bflo_tcas_synt.mpa
bflo_dmel_synt.exp %>% group_by(clg) %>% summarise(mpa=-log10(min(padj))) %>% mutate(species='Drosophila') -> bflo_dmel_synt.mpa

rbind(pcau_bflo_synt.mpa,hspi_bflo_synt.mpa,bflo_gspe_synt.mpa, bflo_cele_synt.mpa,bflo_mhap_synt.mpa,bflo_rvar_synt.mpa,bflo_plit_synt.mpa,
      tant_bflo_synt.mpa,bflo_crot_synt.mpa,bflo_ipac_synt.mpa,sacu_bflo_synt.mpa,tcor_bflo_synt.mpa,spar_bflo_synt.mpa,bflo_tcas_synt.mpa,bflo_dmel_synt.mpa) -> all.mpa

rbind(pcau_bflo_synt.mpa,hspi_bflo_synt.mpa,bflo_gspe_synt.mpa, bflo_cele_synt.mpa,bflo_mhap_synt.mpa,bflo_rvar_synt.mpa,tant_bflo_synt.mpa,bflo_crot_synt.mpa,bflo_ipac_synt.mpa,sacu_bflo_synt.mpa,tcor_bflo_synt.mpa,spar_bflo_synt.mpa,bflo_tcas_synt.mpa,bflo_dmel_synt.mpa) -> all.mpa


all.mpa %>% pivot_wider(names_from='clg',values_from='mpa',values_fill=1) %>% column_to_rownames(var='species') %>% as.matrix() -> all.mpa.mat
#all.map.mat=read.table("all_min_pvals.txt",h=T,row.names=1)
write.table(all.mpa.mat,file="all_min_pvals.txt",quote=F,sep='\t')
all.mpa.mat=all.mpa.mat[rev(1:nrow(all.mpa.mat)), ]
bpal=colorRampPalette(brewer.pal(name='PuBu',n=9))(100)

pheatmap(all.mpa.mat,cluster_rows = F,cluster_cols=F,col=bpal,breaks = seq(0,150, length.out=100))




    