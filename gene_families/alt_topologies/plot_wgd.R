setwd('/Users/fmarletaz/Dropbox (UCL)/Genomes/Priapulid/Analyses/Gene_families/Restricted')
tcp<-read.delim("topology_content.tsv")
tcp %>% 

library(tidyverse)

sel=c('Ecdysozoa','Priapulida','Priapulimorpha','Cryptovermes','Euarthopoda','Mandibulata','Chelicerata','Xiphosura','Myriapoda','Millipede','Pancrustacea','Altocrustacea','Vericrustacea','Insecta','Holometabola','Cru2','Cru1','Triradiata','Tardigrada','Nematoda','Nem_1','Chromadorea','Rhbdidita','Triradiata','Panathropoda')
sel=c(c('Ecdysozoa','Priapulida','Cryptovermes','Euarthopoda','Chelicerata','Myriapoda','Pancrustacea','Insecta','Tardigrada','Nematoda','Panarthropoda','Triradiata'))
m <- read.delim("topology_content_loss.tsv") %>%
  filter(node_a %in% sel,node_b %in% sel) %>%
  select(node_a, node_b, jaccard) %>%
  pivot_wider(names_from = node_b, values_from = jaccard) %>%
  column_to_rownames("node_a") %>%
  as.matrix()


## put nodes shared by both topologies on the diagonal
ord <- intersect(rownames(m), colnames(m))
m <- m[c(ord, setdiff(rownames(m), ord)), c(ord, setdiff(colnames(m), ord))]

pheatmap::pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f",
                   breaks = seq(0, 1, length.out = 101),
                   color = colorRampPalette(
                     RColorBrewer::brewer.pal(7, "GnBu"))(100),
                   fontsize = 14, fontsize_number = 12,
                   main = "Gene family content, Panthropoda vs Triradiata")

l <- read.delim("loss_convergence.tsv") %>%
  filter(relationship %in% c('self','disjoint')) %>%
  filter(node_a!='Cryptovermes',node_b!='Cryptovermes') %>% 
  select(node_a, node_b, jaccard) %>%
  pivot_wider(names_from = node_b, values_from = jaccard) %>%
  column_to_rownames("node_a") %>%
  as.matrix()

pheatmap::pheatmap(l, 
                   breaks = seq(0, 1, length.out = 101),
                   color = colorRampPalette(
                     RColorBrewer::brewer.pal(7, "GnBu"))(100),
                   fontsize = 14, fontsize_number = 12,
                   main = "Joint losses")

ord <- c("Metazoa","Planulozoa","Bilateria","Protostomia","Ecdysozoa",
         "Scalidophora","Priapulida","Priapulidae","Pricau")   # root to tip

d <- read.delim("../Revisions/wgd_placement.tsv") %>%
  filter(node %in% ord) %>%
  mutate(node = factor(node, levels = rev(ord)))

ggplot(d, aes(pct, node)) +
  geom_col(aes(x = pct_expected), fill = "grey90", width = .7) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = .18,
                 colour = "grey25", linewidth = .4) +
  geom_point(aes(size = n_tested), shape = 21, fill = "#2C6FBB",
             colour = "white", stroke = .6) +
  scale_size_area(max_size = 5, name = "pairs tested") +
  labs(x = "WGD-consistent paralogue pairs (%)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

