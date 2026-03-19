require(RColorBrewer)
require(pheatmap)
require(tidyverse)
require(wesanderson)
require(stringr)
#require(ggalluvial)

source("~/Dropbox/Genomes/Synteny/synt-forge/synteny-functions.R")


setwd("~/Dropbox/Genomes/Synteny/synt-forge/Priapulid/")

pmax_col=c("#A6CEE3","#1F78B4","#1B9E77","#33A02C","#B2DF8A","#FB9A99","#E31A1C","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#6A3D9A","#E6AB02","#A6761D","#666666","#FF7F00","#BEBADA","#e9521e")
arub_col=c("#A6CEE3","#1F78B4","#1B9E77","#E31A1C","#FB9A99","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00","#FDBF6F","#BEBADA","#ffd92f","#e9521e")


univ_col=c("#A6CEE3","#1F78B4","#1B9E77","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FFFF99","#B15928","#D95F02",
           "#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00",
           "#FDBF6F","#BEBADA","#ffd92f","#e9521e","lightgrey")

bflo_col=c("#1F78B4","#1B9E77","#33A02C","#B2DF8A","#FB9A99","#E31A1C","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00","#FDBF6F","#BEBADA","#ffd92f")



par(mar=c(5,5,10,5))

pcau_pmax_synt=read.table("Pcau_Pmax_synt.txt",h=T,stringsAsFactors=T)
pcau_bflo_synt=read.table("Pcau_Bflo_synt.txt",h=T,stringsAsFactors=T)
pcau_arub_synt=read.table("Pcau_Arub_synt.txt",h=T,stringsAsFactors=T)

###### CRAB EXCERPT
atan_spar_synt=read.table("Atan-Spar_synt.txt",h=T,stringsAsFactors=T)
plotSynt(ordSynt(atan_spar_synt),'Afruca','Scylla')

###


#plotSyntCLGr(pgot_pmax_synt,'Paraspadella','Pecten')
pcau_pmax_synt.exp<-testEnrich(pcau_pmax_synt)
pcau_bflo_synt.exp<-testEnrich(pcau_bflo_synt)
pcau_arub_synt.exp<-testEnrich(pcau_arub_synt)

pcau.unified<-unify(pcau_bflo_synt.exp,pcau_arub_synt.exp,pcau_pmax_synt.exp)
write.table(pcau.unified,file="Pcau_unified_ALG.tsv",sep="\t",quote=F,row.names=F)

palette(bflo_col)

plotSyntCLGrs(pcau_bflo_synt.exp,'Priapulus','Amphioxus')
legend('bottomright',levels(pcau_bflo_synt.exp$clg),fill=bflo_col)

palette(pmax_col)
plotSyntCLGrs(pcau_pmax_synt.exp,'Priapulus','Pecten')
legend('bottomright',levels(pcau_pmax_synt.exp$clg),fill=pmax_col)

palette(arub_col)
plotSyntCLGrs(pcau_arub_synt.exp,'Priapulus','Asterias')
legend('bottomright',levels(pcau_arub_synt.exp$clg),fill=arub_col)



unify <- function(bflo_synt,arub_synt,pmax_synt){
  bflo_synt$sp<-'Bflo'
  arub_synt$sp<-'Arub'
  pmax_synt$sp<-'Pmax'

  arub_clg=c('A1','A2','D','E','F','G','H','I','J1','J2','K','L','M','N','P','Q','R')
  pmax_clg=c('C1','C2','O1')
  bflo_clg=c('B1','B2','B3','O2')

  pmax_synt %>% filter(sp=='Pmax' & clg %in% pmax_clg) -> pmax_synt.r
  arub_synt %>% filter(sp=='Arub' & clg %in% arub_clg) -> arub_synt.r
  bflo_synt %>% filter(sp=='Bflo' & clg %in% bflo_clg) -> bflo_synt.r

  unified.exp<-rbind(pmax_synt.r,arub_synt.r,bflo_synt.r)
  unified.exp$col<-ifelse(unified.exp$scol=='S',as.character(unified.exp$clg),'Z')
  algz=c('A1','A2','B1','B2','B3','C1','C2','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P','Q','R','Z')
  unified.exp$col<-factor(unified.exp$col, levels = algz)
  
  return(unified.exp)
}
#unified.all.exp %>% filter((sp='Pmax' & clg %in% pmax_clg) | (sp=='Arub' & clg %in% arub_clg) | (sp=='Bflo' & clg %in% bflo_clg)) -> unified.exp
unified<-unify(pcau_bflo_synt.exp,pcau_arub_synt.exp,pcau_pmax_synt.exp)
unified$clg<-factor(unified$clg, levels = algs)
unified$s1chr<-factor(unified$s1chr, levels = c('chr1','chr2','chr3','chr6','chr9','chr11','chr7','chr5','chr8','chr12','chr13','chr10','chr4','chr14'))
unified.exp$s1chr<-factor(unified.exp$s1chr, levels = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','scaf1'))

ggplot(unified, aes(clg,s1chr))+ geom_point(aes(size=n,color=col))+
  theme_light()+ theme(axis.text.x =  element_text(angle = 90))+
  xlab('Priapulus')+ylab('CLG')+scale_color_manual(values = univ_col,breaks = algz)+theme_minimal()
unified.exp$col
names(univ_col)=c('A1','A2','B1','B2','B3','C1','C2','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P','Q','R')
algz=c('A1','A2','B1','B2','B3','C1','C2','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P','Q','R','Z')
algs=c('A1','A2','B1','B2','B3','C1','C2','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P','Q','R')
names(univ_col)<-algz

hspi_pmax_synt=read.table("Hspi_Pmax_synt.txt",h=T,stringsAsFactors=T)
hspi_bflo_synt=read.table("Hspi_Bflo_synt.txt",h=T,stringsAsFactors=T)
hspi_arub_synt=read.table("Hspi_Arub_synt.txt",h=T,stringsAsFactors=T)

hspi_pmax_synt.exp<-testEnrich(hspi_pmax_synt)
hspi_bflo_synt.exp<-testEnrich(hspi_bflo_synt)
hspi_arub_synt.exp<-testEnrich(hspi_arub_synt)

palette(bflo_col)
plotSyntCLGrs(hspi_bflo_synt.exp,'Priapulus','Amphioxus')
legend('bottomright',levels(hspi_bflo_synt.exp$clg),fill=bflo_col)
)
palette(pmax_col)
plotSyntCLGrs(hspi_pmax_synt.exp,'Priapulus','Pecten')
legend('bottomright',levels(hspi_pmax_synt.exp$clg),fill=pmax_col)

palette(arub_col)
plotSyntCLGrs(hspi_arub_synt.exp,'Priapulus','Asterias')
legend('bottomright',levels(hspi_arub_synt.exp$clg),fill=arub_col)


hspi.unified<-unify(hspi_bflo_synt.exp,hspi_arub_synt.exp,hspi_pmax_synt.exp)
hspi.unified$clg<-factor(hspi.unified$clg, levels = algs)

write.table(hspi.unified,file="Hspi_unified_ALG.tsv",sep="\t",quote=F,row.names=F)
#unified$s1chr<-factor(unified$s1chr, levels = c('chr1','chr2','chr3','chr6','chr9','chr11','chr7','chr5','chr8','chr12','chr13','chr10','chr4','chr14'))
unified.exp$s1chr<-factor(unified.exp$s1chr, levels = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','scaf1'))

ggplot(hspi.unified, aes(clg,s1chr))+ geom_point(aes(size=n,color=col))+
  theme_light()+ theme(axis.text.x =  element_text(angle = 90))+
  xlab('CLG')+ylab('Halicryptus')+scale_color_manual(values = univ_col,breaks = algz)+theme_minimal()

hspi.unified %>% mutate(s1chr = gsub("\\|arrow", "", s1chr)) %>% filter(s1chr != 'ScKMIYO_7') -> hspi.unified
head(DT)
sacu_pmax_synt=read.table("Sacu_Pmax_synt.txt",h=T,stringsAsFactors=T)
sacu_bflo_synt=read.table("Sacu_Bflo_synt.txt",h=T,stringsAsFactors=T)
sacu_arub_synt=read.table("Sacu_Arub_synt.txt",h=T,stringsAsFactors=T)

sacu_pmax_synt.exp<-testEnrich(sacu_pmax_synt)
sacu_bflo_synt.exp<-testEnrich(sacu_bflo_synt)
sacu_arub_synt.exp<-testEnrich(sacu_arub_synt)

palette(bflo_col)
plotSyntCLGrs(ordSynt(sacu_bflo_synt.exp),'Amphioxus','Strigamia')
legend('bottomright',levels(sacu_bflo_synt.exp$clg),fill=bflo_col)

palette(pmax_col)
plotSyntCLGrs(ordSynt(sacu_pmax_synt.exp),'Pecten','Strigamia')
legend('bottomright',levels(sacu_pmax_synt.exp$clg),fill=pmax_col)

palette(arub_col)
plotSyntCLGrs(sacu_arub_synt.exp,'Asterias','Strigamia')
legend('bottomright',levels(sacu_arub_synt.exp$clg),fill=arub_col)


sacu.unified<-unify(sacu_bflo_synt.exp,sacu_arub_synt.exp,sacu_pmax_synt.exp)
write.table(sacu.unified,file="Sacu_unified_ALG.tsv",sep="\t",quote=F,row.names=F)
sacu.unified$clg<-factor(sacu.unified$clg, levels = algs)
#unified$s1chr<-factor(unified$s1chr, levels = c('chr1','chr2','chr3','chr6','chr9','chr11','chr7','chr5','chr8','chr12','chr13','chr10','chr4','chr14'))
unified.exp$s1chr<-factor(unified.exp$s1chr, levels = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','scaf1'))

ggplot(sacu.unified, aes(clg,s1chr))+ geom_point(aes(size=n,color=col))+
  theme_light()+ theme(axis.text.x =  element_text(angle = 90))+
  xlab('CLG')+ylab('Strigamia')+scale_color_manual(values = univ_col,breaks = algz)+theme_minimal()

tcor_bflo_synt=read.table("Tcor_Bflo_synt.txt",h=T,stringsAsFactors=T)
tcor_bflo_synt.exp<-testEnrich(tcor_bflo_synt)

sacu.unified<-unify(sacu_bflo_synt.exp,sacu_arub_synt.exp,sacu_pmax_synt.exp)
write.table(sacu.unified,file="Sacu_unified_ALG.tsv",sep="\t",quote=F,row.names=F)
sacu.unified$clg<-factor(sacu.unified$clg, levels = algs)

palette(bflo_col)
plotSyntCLGrs(ordSynt(tcor_bflo_synt.exp),'Amphioxus','Trigoniulus')
legend('bottomright',levels(sacu_bflo_synt.exp$clg),fill=bflo_col)




plotSynt(spar_pmax_synt.exp,'Scylla','Pecten')

tant_pmax_synt=read.table("Tant-Pmax_synt.txt",h=T,stringsAsFactors=T)
tant_bflo_synt=read.table("Tant-Bflo_synt.txt",h=T,stringsAsFactors=T)
tant_arub_synt=read.table("Tant-Arub_synt.txt",h=T,stringsAsFactors=T)

tant_pmax_synt.exp<-testEnrich(tant_pmax_synt)
tant_bflo_synt.exp<-testEnrich(tant_bflo_synt)
tant_arub_synt.exp<-testEnrich(tant_arub_synt)

tant.unified<-unify(tant_bflo_synt.exp,tant_arub_synt.exp,tant_pmax_synt.exp)
write.table(tant.unified,file="Tant_unified_ALG.tsv",sep="\t",quote=F,row.names=F)
tant.unified$clg<-factor(tant.unified$clg, levels = algs)

ggplot(tant.unified, aes(clg,s1chr))+ geom_point(aes(size=n,color=col))+
  theme_light()+ theme(axis.text.x =  element_text(angle = 90))+
  xlab('CLG')+ylab('Strigamia')+scale_color_manual(values = univ_col,breaks = algz)+theme_minimal()

plotSyntCLGrs(ordSynt(tant_pmax_synt.exp),'Pecten','Trichonephila')


###
#### RETENTION RATES #####
#skChrAsgn=read.table("skate_clg_assign.txt",stringsAsFactors = T)
#names(skChrAsgn)<-c('Chr','Cat','Lab')

pcau_bflo_synt %>% #mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) %>% 
  group_by(clg) %>% tally %>%  column_to_rownames(var='clg')-> clg_ct


pcau_bflo_synt.exp %>% #mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) %>% 
  filter(scol=='S') %>% group_by(clg,s1chr) %>% tally -> sig_clg_ct
sig_clg_ct %>% mutate(ret=n/clg_ct[as.vector(clg),])-> sig_clg_ret

ggplot(sig_clg_ret,aes(x=ret,fill=clg))+geom_histogram(alpha = 0.8)+
  geom_density(alpha = 0.1)+scale_fill_manual(values=bflo_col)+
  theme_minimal()+facet_wrap(~clg,scales='free_y')

ggplot(sig_clg_ret,aes(x=ret))+geom_histogram(alpha = 0.8)+
  geom_density(alpha = 0.1)+scale_fill_manual(values=bflo_col)+
  theme_minimal()#+facet_wrap(~clg,scales='free_y')


#skChrAsgn %>% column_to_rownames(var='Chr') -> skChrAsgn
sk_clg_ct %>% mutate(cat=skChrAsgn[as.vector(s2chr),]$Cat) %>%
  separate(cat, c("l","clgf","Evt"), sep = "_") %>% select(-l) %>%
  mutate(prz=grepl(clg,clgf)) %>% filter(prz==T) -> sk_clg_ct


sk_clg_ct %>% mutate(dup=str_sub(Evt,-1,-1)) -> sk_clg_ct

sk_clg_ct %>% group_by(s2chr,clgf,Evt) %>% filter(ret==max(ret)) -> sk_clg_flt
ggplot(sk_clg_ret,aes(x=ret,fill=dup,color=dup))+geom_density(alpha = 0.5)+theme_light()
ggplot(sk_clg_ret,aes(x=ret,fill=dup,color=dup))+geom_histogram(alpha = 0.5)+facet_grid(rows=vars(dup)
                                                                                        
                                                                                        #### 



# Crab
spar_pmax_synt=read.table("Spar-Pmax_synt.txt",h=T,stringsAsFactors=T)
spar_bflo_synt=read.table("Spar-Bflo_synt.txt",h=T,stringsAsFactors=T)
spar_arub_synt=read.table("Spar-Arub_synt.txt",h=T,stringsAsFactors=T)

spar_pmax_synt.exp<-testEnrich(spar_pmax_synt)
spar_bflo_synt.exp<-testEnrich(spar_bflo_synt)
spar_arub_synt.exp<-testEnrich(spar_arub_synt)
spar.unified<-unify(spar_bflo_synt.exp,spar_arub_synt.exp,spar_pmax_synt.exp)
write.table(spar.unified,file="Spar_unified_ALG.tsv",sep="\t",quote=F,row.names=F)


plotSyntCLGrs(spar_arub_synt.exp,'Scylla','Asterias')


######## PREPARATION RIBBON plot

#### Pcau vs Arub 

arub_bflo_synt=read.table("Arub_Bflo_synt.txt",h=T,stringsAsFactors=T)
arub_bflo_synt.exp <- testEnrich(arub_bflo_synt)
arub_bflo_synt.exp %>% filter(scol=='S') %>% select(s1g,s1chr,clg) %>% mutate(s2g=s1g,sxchr=s1chr,clgb=clg) %>% select(-s1g,-clg,-s1chr) -> arub_bflo_synt.expf
inner_join(pcau_arub_synt.exp,arub_bflo_synt.expf) %>% relocate(clgb,.after=sxchr) -> pcau_arub_synt.expr
levels(pcau_arub_synt.expr$clg)
levels(pcau_arub_synt.expr$clgb)

pcau_arub_synt.expr %>% mutate(clgc=if_else(clg %in% c("B3","B2+C2","B1","C1","O1","O2+B3"),clgb,clg)) -> pcau_arub_synt.expr
write.table(pcau_arub_synt.expr,file='Pcau_Arub_syntr.txt',sep='\t',quote=F,row.names=F)

#bflo_pmax_synt.exp %>% filter(scol=='S') %>% select(s2g,s2chr,clg) %>% mutate(sxchr=s2chr,clgb=clg) %>% select(-clg,-s2chr) -> bflo_pmax_synt.expf
inner_join(pgot_pmax_synt.exp,bflo_pmax_synt.expf) %>% relocate(clgb,.after=sxchr) %>% select(-sp) -> pgot_pmax_synt.expr
pgot_pmax_synt.expr %>% mutate(clgc=if_else(clg %in% c("B2+M","H+Q","L+J","C","J","O"),clgb,clg)) -> pgot_pmax_synt.expr
write.table(pgot_pmax_synt.expr,file='Pgot_Pmax_syntr.txt',sep='\t',quote=F,row.names=F)
pcau_arub_synt.exp

#### Pcau vs Pmax
pcau_pmax_synt=read.table("Pcau_Pmax_synt.txt",h=T,stringsAsFactors=T)
pcau_pmax_synt.exp <- testEnrich(pcau_pmax_synt)
inner_join(pcau_pmax_synt.exp,bflo_pmax_synt.expf) %>% relocate(clgb,.after=sxchr) -> pcau_pmax_synt.expr
pcau_pmax_synt.expr %>% mutate(clgc=if_else(clg %in% c("B2+M","H+Q","L+J","C","J","O"),clgb,clg)) -> pcau_pmax_synt.expr
write.table(pcau_pmax_synt.expr,file='Pcau_Pmax_syntr.txt',sep='\t',quote=F,row.names=F)
pcau_pmax_synt.expr %>% select(s1g,s1chr,clgc) %>% mutate(sxchr=s1chr) %>% select(-s1chr) -> pcau_pmax_synt.expf

#### Pcau vs Hspi
pcau_hspi_synt=read.table("Pcau_Hspi_synt.txt",h=T,stringsAsFactors=T)
pcau_hspi_synt.exp <- testEnrich(pcau_hspi_synt)
hspi_uni=read.table("Hspi_unified_ALG.tsv",h=T)
hspi_uni %>% select(s1g,s1chr,col) %>%  mutate(sxchr=s1chr,s2g=s1g,clgx=col) %>% select(-s1chr,-s1g,-col) -> hspi_unir
inner_join(pcau_hspi_synt.exp,hspi_unir) %>% relocate(clgx,.after=sxchr) -> pcau_hspi_synt.expr
pcau_uni=read.table("Pcau_unified_ALG.tsv",h=T)
pcau_uni %>% select(s1g,s1chr,col) %>%  mutate(sychr=s1chr,clgy=col) %>% select(-s1chr,-col) -> pcau_unir
inner_join(pcau_hspi_synt.expr,pcau_unir) %>% relocate(clgy,.after=sychr) -> pcau_hspi_synt.exprr


#inner_join(pcau_hspi_synt.exp,pcau_pmax_synt.expf) %>% relocate(clgc,.after=sxchr) -> pcau_hspi_synt.expr
write.table(pcau_hspi_synt.exprr,file='Pcau_Hspi_syntr2.txt',sep='\t',quote=F,row.names=F)

#### Tcor vs Sacu
sacu_tcor_synt=read.table("Stracu-Tricor_synt.txt",h=T,stringsAsFactors=T)
sacu_tcor_synt.exp <- testEnrich(sacu_tcor_synt)
plotSynts(ordSynt(sacu_tcor_synt.exp),'Strigamia','Triconiulus')
tcor_bflo_synt.exp %>% select(s1g,s1chr,clg) %>% mutate(sxchr=s1chr,s2g=s1g,clgc=clg) %>% select(-s1chr,-s1g,-clg) -> tcor_bflo_synt.expf

inner_join(sacu_tcor_synt.exp,tcor_bflo_synt.expf) %>% relocate(clgc,.after=sxchr) %>% select(-clg) -> sacu_tcor_synt.expr
write.table(sacu_tcor_synt.expr,file='Sacu_Tcor_syntr.txt',sep='\t',quote=F,row.names=F)



