require(RIdeogram)
require(tidyverse)
library("wesanderson")


`%notin%` <- Negate(`%in%`)

setwd("~/Dropbox/Genomes/Synteny/synt-forge/Priapulid/")

colors=data.frame(row.names = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q'),
                  col=c("1F78B4","33A02C","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A","66A61E","6A3D9A","E6AB02",
                        "A6761D","8DD3C7","666666","FF7F00","BEBADA","ffd92f"))

colors=data.frame(row.names = c('A','A1','A2','B1','B2+M','B3','C1','C2','D','E','F','G','H+Q','I','J','K+O2','L+J2','N','O1','P'),
                  col=c("A6CEE3","A6CEE3","1F78B4","1B9E77","33A02C","B2DF8A","FB9A99","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A","66A61E","6A3D9A","E6AB02","A6761D","666666","FF7F00","BEBADA"))

colors=data.frame(row.names = c('A','A1','A2','B1','B2','B3','C1','C2','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P','Q'),
                  col=c("A6CEE3","A6CEE3","1F78B4","1B9E77","33A02C","B2DF8A","FB9A99","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A","66A61E","6A3D9A","CAB2D6","E6AB02","A6761D","8DD3C7","666666","FF7F00","FDBF6F","BEBADA","ffd92f"))

univ_col=c("#A6CEE3","#1F78B4","#1B9E77","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FFFF99","#B15928","#D95F02",
           "#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00",
           "#FDBF6F","#BEBADA","#ffd92f","#e9521e","lightgrey")



# EXAMPLE

data(karyotype_dual_comparison, package="RIdeogram")
head(karyotype_dual_comparison)
data(synteny_dual_comparison, package="RIdeogram")
head(synteny_dual_comparison)

##### Pcau vs Hspi #####
 
lg=read.table("Pcau_Hspi_syntrr.txt",h=T)

chrf=read.table("~/Dropbox/Genomes/Synteny/synt-forge/geneInfo/Hali_chr.txt")
chrn=as.vector(chrf$V3)
names(chrn)=as.vector(chrf$V1)
lg$s2chrn=chrn[lg$s2chr]
lg$s1chrn<-NULL
na.omit

lg %>% filter(s2chrn=='chr8' & scol=='S')

write.table(lg,file="Pcau_Hspi_syntrr.txt",sep='\t',quote=F)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Priapulus',fill="Priapulus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pcau

lg  %>% group_by(s2chrn) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chrn) %>% 
  mutate(Start=1,species='Halicryptus',fill="Halicryptus",size=12,color=252525) %>% relocate(Start,.after=Chr) %>% drop_na()-> k_hspi

k_pcau_ord=c("chr1","chr2","chr3","chr11","chr4","chr6","chr5","chr7","chr8","chr9","chr10","chr13","chr12","chr14")
k_hspi_ord=c("chr8","chr5","chr3","chr2","chr6","chr4","chr1","chr7","chr9","chr11","chr10","chr13","chr12")

k_pcau=k_pcau[order(match(k_pcau$Chr,k_pcau_ord)),]
k_hspi=k_hspi[order(match(k_hspi$Chr,k_hspi_ord)),]

rbind(k_pcau,k_hspi) -> chroms

lg$Species_1=match(lg$s1chr,k_pcau$Chr)
lg$Species_2=match(lg$s2chrn,k_hspi$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clgc) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgc,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgc) -> lgt

ideogram(karyotype = chroms, synteny = lgt)

##### Sacu vs Tcor
lg=read.table("Sacu_Tcor_syntr.txt",h=T)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Strigamia',fill="Strigamia",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_sacu

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Trigoniulus',fill="Trigoniulus",size=12,color=252525) %>% relocate(Start,.after=Chr) %>% drop_na()-> k_tcor

k_sacu_ord=c("OX442315.1","OX442314.1","OX442320.1","OX442319.1","OX442322.1","OX442318.1","OX442321.1","OX442317.1","OX442316.1","OX442323.1","OX442324.1")
k_tcor_ord=rev(c("Scaffold_1043","Scaffold_6978","Scaffold_1977","Scaffold_6181","Scaffold_9089","Scaffold_9407","Scaffold_4307","Scaffold_6695",
             "Scaffold_5748","Scaffold_4228","Scaffold_9389","Scaffold_7780","Scaffold_9406"))

k_sacu=k_sacu[order(match(k_sacu$Chr,k_sacu_ord)),]
k_tcor=k_tcor[order(match(k_tcor$Chr,k_tcor_ord)),]

rbind(k_sacu,k_tcor) -> chroms

lg$Species_1=match(lg$s1chr,k_sacu$Chr)
lg$Species_2=match(lg$s2chr,k_tcor$Chr)

lg %>% filter(scol=='S') %>% filter(clgc!='Z') %>% select(Species_1,s1gp,Species_2,s2gp,clgc) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgc,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgc) -> lgt
lg %>% filter(s1chr=='OX442321.1')
ideogram(karyotype = chroms, synteny = lgt)


##### Pcau vs Pmax

lg=read.table("Pcau_Pmax_syntr.txt",h=T)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Priapulus',fill="Priapulus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pcau

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Pecten',fill="Pecten",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pmax

#k_pcau_ord=c("chr1","chr2","chr9","chr6","chr3","chr4","chr5","chr7","chr8","chr10","chr11","chr12","chr14")
k_pcau_ord=c("chr1","chr2","chr3","chr11","chr9","chr6","chr4","chr5","chr7","chr8","chr10","chr12","chr14")

k_pmax_ord=c("scaffold_4","scaffold_12","scaffold_8","scaffold_6","scaffold_19","scaffold_5","scaffold_1","scaffold_14","scaffold_17","scaffold_3",
             "scaffold_13","scaffold_2","scaffold_15","scaffold_11",
             "scaffold_10","scaffold_16","scaffold_7",
             "scaffold_9","scaffold_18")

k_pcau=k_pcau[order(match(k_pcau$Chr,k_pcau_ord)),]
k_pmax=k_pmax[order(match(k_pmax$Chr,k_pmax_ord)),]

rbind(k_pcau,k_pmax) -> chroms

lg$Species_1=match(lg$s1chr,k_pcau$Chr)
lg$Species_2=match(lg$s2chr,k_pmax$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clgc) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgc,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgc) -> lgt

ideogram(karyotype = chroms, synteny = lgt)

#### Pcau vs Hali
lg=read.table("Pcau_Hspi_syntr.txt",h=T)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Priapulus',fill="Priapulus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pcau

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Halicryptus',fill="Halicryptus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_hspi

k_pcau_ord=c("chr1","chr2","chr3","chr11","chr9","chr6","chr4","chr5","chr7","chr8","chr10","chr12","chr14")
k_hspi_ord=c("ScKMIYO_426","ScKMIYO_86","ScKMIYO_1","ScKMIYO_108","ScKMIYO_414","ScKMIYO_436","ScKMIYO_412","ScKMIYO_104",
             "ScKMIYO_45","ScKMIYO_40","ScKMIYO_96","ScKMIYO_130","ScKMIYO_140")

k_pcau=k_pcau[order(match(k_pcau$Chr,k_pcau_ord)),]
k_hspi=k_hspi[order(match(k_hspi$Chr,k_hspi_ord)),]

rbind(k_pcau,k_hspi) -> chroms

lg$Species_1=match(lg$s1chr,k_pcau$Chr)
lg$Species_2=match(lg$s2chr,k_hspi$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clgc) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgc,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgc) -> lgt

ideogram(karyotype = chroms, synteny = lgt)



####### Pgot vs Pmax
lg=read.table("Pgot_Pmax_syntr.txt",h=T)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Paraspadella',fill="Paraspadella",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pgot

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Pecten',fill="Pecten",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pmax

k_pgot_ord=c("chr1","chr7","chr2","chr3","chr4","chr5","chr6","chr8","chr9","scaf1")
k_pmax_ord=c("scaffold_19","scaffold_1","scaffold_14","scaffold_17","scaffold_12","scaffold_8","scaffold_3",
             "scaffold_6","scaffold_13","scaffold_2","scaffold_15","scaffold_11",
             "scaffold_10","scaffold_5","scaffold_4","scaffold_16","scaffold_7",
             "scaffold_9","scaffold_18")

k_pgot=k_pgot[order(match(k_pgot$Chr,k_pgot_ord)),]
k_pmax=k_pmax[order(match(k_pmax$Chr,k_pmax_ord)),]

rbind(k_pgot,k_pmax) -> chroms

lg$Species_1=match(lg$s1chr,k_pgot$Chr)
lg$Species_2=match(lg$s2chr,k_pmax$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clgc) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgc,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgc) -> lgt

ideogram(karyotype = chroms, synteny = lgt)
##### Pgot vs Avag

lg=read.table("Pgot_Avag_syntr.txt",h=T)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Paraspadella',fill="Paraspadella",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pgot

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Adinetta',fill="Adinetta",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_avag

k_pgot_ord=c("chr1","chr7","chr2","chr3","chr4","chr5","chr6","chr8","chr9","scaf1")
k_avag_ord=c('Chrom_2A','Chrom_5A','Chrom_6A','Chrom_1A','Chrom_4A','Chrom_3A')

k_pgot=k_pgot[order(match(k_pgot$Chr,k_pgot_ord)),]
k_avag=k_avag[order(match(k_avag$Chr,k_avag_ord)),]

rbind(k_pgot,k_avag) -> chroms

lg$Species_1=match(lg$s1chr,k_pgot$Chr)
lg$Species_2=match(lg$s2chr,k_avag$Chr)
#filter(scol=='S') %>% 
lg %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

ideogram(karyotype = chroms, synteny = lgt)


##### Pmax vs Ofus #####

lg=read.table("Pmax_Ofus_syntr.txt",h=T)

lg %>% filter(!(s1chr!='PMA9' & clgb=='C1')) %>%
  filter(!(s1chr=='PMA10' & clgb=='I')) %>%
  filter(!(s1chr!='PMA17' & clgb=='C2')) %>%
  filter(!(s1chr!='PMA19' & clgb=='B3')) %>%
  filter(!(s1chr!='PMA15' & clgb=='B1')) %>%
  filter(!(s1chr!='PMA1' & clgb=='Q')) -> lg

lg %>% filter(s1chr=='PMA4') %>% group_by(clgb) %>% tally
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Pecten',fill="Pecten",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pmax

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Owenia',fill="owenia",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_ofus


k_ofus_ord=c('OFU1','OFU3','OFU2','OFU5','OFU7','OFU4','OFU6','OFU8','OFU11','OFU9','OFU10','OFU12')

k_pmax_ord=c('PMA4','PMA16','PMA18','PMA2','PMA10','PMA12','PMA7','PMA3','PMA13','PMA9','PMA15','PMA19','PMA17','PMA1','PMA5','PMA11','PMA8','PMA6','PMA14')

k_ofus=k_ofus[order(match(k_ofus$Chr,k_ofus_ord)),]
k_pmax=k_pmax[order(match(k_pmax$Chr,k_pmax_ord)),]

rbind(k_pmax,k_ofus) -> chroms


lg$Species_1=match(lg$s1chr,k_pmax$Chr)
lg$Species_2=match(lg$s2chr,k_ofus$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clgb) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgb,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgb) -> lgt

ideogram(karyotype = chroms, synteny = lgt)



##### Bflo vs Pmax #####

lg=read.table("Bflo_Pmax_syntr.txt",h=T)
head(lg)

lg %>% filter(!(s2chr=='PMA4' & clg=='C1')) %>%
       filter(!(s2chr=='PMA10' & clg=='B2')) %>%
       filter(!(s2chr=='PMA1' & clg=='C2')) -> lg
  
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Branchiostoma',fill='Branchiostoma',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_bflo

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Pecten',fill='Pecten',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pmax

k_bflo_ord=c('BFL15','BFL2','BFL16','BFL14','BFL17','BFL8','BFL4','BFL5','BFL1','BFL10','BFL18','BFL3','BFL13','BFL19','BFL9','BFL6','BFL11','BFL7','BFL12')
k_pmax_ord=c('PMA4','PMA16','PMA18','PMA2','PMA10','PMA12','PMA7','PMA3','PMA13','PMA9','PMA15','PMA19','PMA17','PMA1','PMA5','PMA11','PMA8','PMA6','PMA14')

k_bflo=k_bflo[order(match(k_bflo$Chr,k_bflo_ord)),]
k_pmax=k_pmax[order(match(k_pmax$Chr,k_pmax_ord)),]

rbind(k_bflo,k_pmax) -> chroms

lg$Species_1=match(lg$s1chr,k_bflo$Chr)
lg$Species_2=match(lg$s2chr,k_pmax$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

ideogram(karyotype = chroms, synteny = lgt)

#### Ofus vs Sben ####

lg=read.table("Ofus_Sben_syntr.txt",h=T)

head(lg)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Owenia',fill="Owenia",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_ofus

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Streblospio',fill="Streblospio",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_sben


k_ofus_ord=c('OFU1','OFU3','OFU2','OFU5','OFU7','OFU4','OFU6','OFU8','OFU11','OFU9','OFU10','OFU12')

k_sben_ord=c("SBE6","SBE8","SBE7","SBE9","SBE2","SBE5","SBE10","SBE4","SBE1","SBE11","SBE3")

k_ofus=k_ofus[order(match(k_ofus$Chr,k_ofus_ord)),]
k_sben=k_sben[order(match(k_sben$Chr,k_sben_ord)),]

rbind(k_ofus,k_sben) -> chroms


lg$Species_1=match(lg$s1chr,k_ofus$Chr)
lg$Species_2=match(lg$s2chr,k_sben$Chr)
#head(lg)
#levels(as.factor(lg$clgb))
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clgb) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgb,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgb) -> lgt

head(lg)

ideogram(karyotype = chroms, synteny = lgt)

##### Pmax vs Llon #####

lg=read.table("Pmax_Llon_syntr.txt",h=T)
head(lg)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Pecten',fill="Pecten",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pmax

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Lineus',fill="Lineus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_llon

k_pmax_ord=c('PMA4','PMA16','PMA18','PMA2','PMA10','PMA12','PMA7','PMA3','PMA13','PMA9','PMA15','PMA19','PMA17','PMA1','PMA5','PMA11','PMA8','PMA6','PMA14')
k_llon_ord=c('LLO4','LLO12','LLO19','LLO13','LLO17','LLO8','LLO6','LLO7','LLO3','LLO18','LLO2','LLO16','LLO15','LLO14','LLO1','LLO5','LLO11','LLO9','LLO10')

k_pmax=k_pmax[order(match(k_pmax$Chr,k_pmax_ord)),]
k_llon=k_llon[order(match(k_llon$Chr,k_llon_ord)),]

rbind(k_pmax,k_llon) -> chroms

lg$Species_1=match(lg$s1chr,k_pmax$Chr)
lg$Species_2=match(lg$s2chr,k_llon$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clgb) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clgb,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clgb) -> lgt

ideogram(karyotype = chroms, synteny = lgt)


##############
##### Leri vs Ggal #####

lg=read.table("Leri-Ggal_syntr.txt",h=T)
head(lg)
lg %>% group_by(s2chr) %>% summarise(size=max(s2gp)) -> gal.size

lg %>% mutate(s2chr=paste('GGA',s2chr,sep="")) %>%
  mutate(s1chr=gsub('Leri_','LER',gsub('C','',s1chr))) %>%
  mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) -> lg

lg %>% filter(s1chr %notin% c('LER50','LER45','LER44')) %>%
  filter(s2chr %notin% c('GGA30','GGA16','GGAW','GGA31')) -> lg 



lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Leucoraja',fill='Leucoraja',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_ler

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Gallus',fill='Gallus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_gal

k_gal_ord=c("GGA4","GGA5","GGA33","GGAZ","GGA32","GGA2","GGA3","GGA22","GGA9","GGA1","GGA26","GGA10","GGA6","GGA12","GGA7","GGA1","GGA8","GGA18","GGA23","GGA13","GGA14","GGA11","GGA15","GGA19","GGA20","GGA17","GGA21","GGA28","GGA27","GGA24","GGA25")
#k_gal_ord=c("GGA4","GGAZ","GGA2","GGA3","GGA5","GGA9","GGA31","GGA1","GGA10","GGA6","GGA12","GGA8","GGA30","GGA18","GGA23","GGA13","GGA7","GGA14","GGA11","GGA15","GGA19","GGA20","GGA17","GGA21","GGA28","GGAW","GGA26","GGA25","GGA22","GGA32","GGA27","GGA24","GGA33")
#k_ler_ord=c("LER12","LER1","LER3","LER2","LER4","LER8","LER5","LER18","LER9","LER14","LER37","LER19","LER22","LER13","LER6","LER33","LER36","LER34","LER15","LER16","LER10","LER39","LER23","LER26","LER11","LER7","LER20","LER17","LER25","LER28","LER21","LER31","LER30","LER29","LER24","LER35","LER27","LER38","LER32","LER40")
k_ler_ord=c("LER12","LER9","LER18","LER40","LER1","LER3","LER38","LER2","LER4","LER8","LER5","LER35","LER14","LER13","LER6","LER24","LER33","LER36","LER34","LER15","LER16","LER7","LER37","LER19","LER22","LER10","LER23","LER26","LER39","LER11","LER20","LER17","LER25","LER28","LER21","LER31","LER30","LER29","LER27","LER32")

k_ler=k_ler[order(match(k_ler$Chr,k_ler_ord)),]
k_gal=k_gal[order(match(k_gal$Chr,k_gal_ord)),]

rbind(k_ler,k_gal) -> chroms

lg$Species_1=match(lg$s1chr,k_ler$Chr)
lg$Species_2=match(lg$s2chr,k_gal$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

ideogram(karyotype = chroms, synteny = as.data.frame(lgt))

to_invert=c('LER1','LER10','LER2','LER9','LER22')
lg %>% mutate(inv=ifelse(s1chr %in% to_invert,'Y','N')) -> lg
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg,s1chr,inv) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  group_by(s1chr) %>% mutate(Start_1=ifelse(inv=='Y',max(Start_1)-Start_1,Start_1)) %>% ungroup() %>% select(-s1chr,-inv) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) %>%
  arrange(Species_1,Start_1,Species_2,Start_2) -> lgt

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

ideogram(karyotype = chroms, synteny = as.data.frame(lgt))




##### Locu vs Ggal #####

lg3=read.table("Locu-Ggal_syntra2.txt",h=T)
head(lg3)


#lg %>% filter(s1chr %notin% c('LER50','LER45','LER44')) %>%
#  filter(s2chr %notin% c('LOC24','LOC28')) -> lg 

lg3 %>% mutate(s2chr=paste('GGA',s2chr,sep="")) %>%
       mutate(s1chr=gsub('LG','LOC',s1chr)) %>%
       mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) %>%
       filter(s1chr %notin% c('LOC28','LOC24')) %>%
        filter(s2chr %notin% c('GGA30','GGA16','GGAW','GGA31','GGA25')) -> lg3

lg3 %>% filter(s1chr=='LOC6') %>% group_by(s2chr) %>% tally
               

lg3 %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Lepisosteus',fill='Lepisosteus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_loc

lg3  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Gallus',fill='Gallus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_gal

#k_gal_ord=c("GGA5","GGA33","GGA4","GGAZ","GGA32","GGA2","GGA22","GGA3","GGA9","GGA1","GGA26","GGA10","GGA6","GGA12","GGA8","GGA18","GGA23","GGA13","GGA7","GGA14","GGA11","GGA15","GGA19","GGA20","GGA17","GGA21","GGA28","GGA27","GGA24","GGA25","GGA16","GGAW")
#k_gal_ord=c("GGA5","GGA4","GGA33","GGAZ","GGA32","GGA2","GGA3","GGA22","GGA9","GGA1","GGA26","GGA10","GGA6","GGA12","GGA7","GGA1","GGA8","GGA18","GGA23","GGA13","GGA14","GGA11","GGA15","GGA19","GGA20","GGA17","GGA21","GGA28","GGA27","GGA24","GGA25")
k_gal_ord=c("GGA5","GGA4","GGA33","GGAZ","GGA32","GGA2","GGA3","GGA22","GGA9","GGA1","GGA26","GGA10","GGA6","GGA12","GGA7","GGA1","GGA8","GGA18","GGA23","GGA13","GGA14","GGA11","GGA15","GGA19","GGA20","GGA17","GGA21","GGA28","GGA27","GGA24")

#k_loc_ord=c('LOC27','LOC7','LOC4','LOC2','LOC11','LOC9','LOC16','LOC1','LOC14','LOC8','LOC17','LOC3','LOC5','LOC10','LOC6','LOC12','LOC13','LOC23','LOC20','LOC22','LOC18','LOC21','LOC25','LOC19','LOC15','LOC26','LOC24','LOC28')
#k_loc_ord=c("LOC27","LOC7","LOC4","LOC2","LOC9","LOC11","LOC16","LOC1","LOC14","LOC8","LOC17","LOC3","LOC5","LOC12","LOC10","LOC6","LOC13","LOC23","LOC20","LOC22","LOC18","LOC21","LOC25","LOC19","LOC15","LOC26","LOC24")
k_loc_ord=c("LOC27","LOC7","LOC4","LOC2","LOC9","LOC11","LOC16","LOC1","LOC14","LOC8","LOC17","LOC3","LOC5","LOC12","LOC10","LOC6","LOC13","LOC23","LOC20","LOC22","LOC18","LOC21","LOC25","LOC19","LOC15","LOC26")



k_gal=k_gal[order(match(k_gal$Chr,k_gal_ord)),]
k_loc=k_loc[order(match(k_loc$Chr,k_loc_ord)),]

rbind(k_loc,k_gal) -> chroms

lg3$Species_1=match(lg3$s1chr,k_loc$Chr)
lg3$Species_2=match(lg3$s2chr,k_gal$Chr)
head(lg3)

lg3 %>% mutate(scol=replace(scol,s1chr=='LOC14' & s2chr=='GGA1','S')) -> lg3

lg3 %>% filter(s1chr=='LOC6') %>% group_by(s2chr) %>% tally

lg3 %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt3

ideogram(karyotype = chroms, synteny = lgt3)

######
to_invert=c('GGA1')
lg3 %>% mutate(inv=ifelse(s2chr %in% to_invert,'Y','N')) -> lg3
lg3 %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg,s2chr,inv) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  group_by(s2chr) %>% mutate(Start_2=ifelse(inv=='Y',max(Start_2)-Start_2,Start_2)) %>% ungroup() %>% select(-s2chr,-inv) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) %>%
  arrange(Species_1,Start_1,Species_2,Start_2) -> lgt3

ideogram(karyotype = chroms, synteny = as.data.frame(lgt3))

#lg %>% select(clg)


##### sizes ####

chroms=read.table("/Users/fmarletaz/Dropbox (UCL)/Genomes/Skate/Skate_genome/data/micro_size.txt")
names(chroms)=c('chrom','size','sp')

ggplot(data=chroms, aes(x=chrom, y=size, group=sp,color=sp)) +
  geom_line()+
  geom_point()+theme_bw()





