

go<-read.delim("OrthoPanarthLoss_GO.txt") #%>% filter(fold>1,p_enrichment<.05)
go %>% filter(Enrich=='e')


jl<-read.delim("go_joint.tsv")
jl %>% filter(p_enrichment<0.05,fold>1,GO %in% as.vector(sel)) -> jls
sel=c("GO:1990716","GO:0016101","GO:0001523","GO:0034080","GO:0034340","GO:0005179","GO:0036156","GO:0048496","GO:2001057","GO:0008063","GO:1990879","GO:1990879","GO:0046415","GO:0016101","GO:0001523")

jls
d <- jls %>%
  mutate(
         logp = -log10(pmax(p_enrichment, 1e-300)),
         name = fct_reorder(name, logp, .fun = function(x) max(x, na.rm = TRUE)))

ggplot(d, aes(set, name)) +
  geom_point(aes(size = n_joint_with, fill = logp),shape = 21, color = "black", stroke = 0.3) +
  theme_minimal(base_size = 9) +theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_viridis_c(option = "magma") +   scale_size_continuous(range = c(1, 8))

ggplot(d, aes(set, name)) +
  geom_point(aes(size = n_joint_with, fill = fold),
             shape = 21, colour = "grey35", stroke = .25, na.rm = TRUE) +
  geom_point(data = filter(d, sig), aes(size = n_joint_with),
             shape = 21, fill = NA, colour = "black", stroke = .8, na.rm = TRUE) +
  scale_fill_gradient2(low = "#3D6DA8", mid = "grey94", high = "#C1444A",
                       midpoint = 1, name = "fold\nenrichment",
                       na.value = "grey96") +
  scale_size_area(max_size = 6, name = "families lost\njointly") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       caption = sprintf("outlined, FDR < %.2f; unoutlined, tested but not significant", ALPHA)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major = element_line(colour = "grey93", linewidth = .3),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 0),
        axis.text.y = element_text(size = 7.5),
        plot.caption = element_text(size = 6, colour = "grey45"),
        legend.key.size = unit(.4, "cm"),
        legend.position = "right")

wes_palette("Zissou1")
d <- g %>%
  filter(name %in% sel) %>%
  complete(name = sel, set = unique(g$set)) %>%
  mutate(logq = -log10(pmax(q_enrichment, 1e-300)),
         sig  = !is.na(q_enrichment) & q_enrichment < ALPHA,
         name = fct_reorder(name, logq, .fun = function(x) max(x, na.rm = TRUE)))

p <- ggplot(d, aes(set, name)) +
  geom_point(aes(size = n_joint_with, fill = fold),
             shape = 21, colour = "grey35", stroke = .25, na.rm = TRUE) +
  geom_point(data = filter(d, sig), aes(size = n_joint_with),
             shape = 21, fill = NA, colour = "black", stroke = .8, na.rm = TRUE) +
  scale_fill_gradient2(low = "#3D6DA8", mid = "grey94", high = "#C1444A",
                       midpoint = 1, name = "fold\nenrichment",
                       na.value = "grey96") +
  scale_size_area(max_size = 6, name = "families lost\njointly") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       caption = sprintf("outlined, FDR < %.2f; unoutlined, tested but not significant", ALPHA)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major = element_line(colour = "grey93", linewidth = .3),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 0),
        axis.text.y = element_text(size = 7.5),
        plot.caption = element_text(size = 6, colour = "grey45"),
        legend.key.size = unit(.4, "cm"),
        legend.position = "right")
