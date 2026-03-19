require(tidyverse)
library(ggridges)

setwd("~/Dropbox/Genomes/Priapulid/Analyses/dating_graham/")


age=read_tsv("dat_res_clg_3.txt")

scld=c('Ecdysozoa','Priapulida','Priapulida/Priapulida',
       'Panarthropoda','Euarthopoda',
       'Myriapoda','Chelicerata','Crustacea')

age %>% filter(clade %in% scld) -> rage

rage$clade=factor(rage$clade,levels=scld)

ggplot(rage , aes(x = -age, y = clade)) + xlim(-700,-300)+geom_vline(xintercept = -520)+
  geom_density_ridges()+theme_minimal()


##

scld=c('Priapulida','Priapulidomorpha','Ecdysozoa','Nematoda','Tardigrada','Panarthropoda','Euarthopoda',
       'Mandibulata','Pancrustacea','Insecta','Myriapoda','Chelicerata',
       'Priapulida/Priapulida','Xiphosura/Xiphosura')


