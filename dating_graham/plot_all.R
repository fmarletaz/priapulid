library(ggtree); library(ape); library(ggplot2)

files <- list.files(pattern = "_sample\\.chronogram$")

pdf("family_chronograms.pdf", width = 8, height = 6)   # one page per family
for (f in files) {
  tr <- read.tree(f)
  
  ci  <- do.call(rbind, lapply(strsplit(tr$node.label, "_"), function(z)
    if (length(z) == 2) rev(as.numeric(z)) else c(NA, NA)))   # -> lo, hi
  dat <- data.frame(node = Ntip(tr) + seq_along(tr$node.label))
  dat$range <- Map(c, -ci[, 2], -ci[, 1])                # revts: negative ages
  
  print(
    revts(ggtree(tr, linewidth = .4)) %<+% dat +
      geom_range("range", colour = "#6f8fb5", size = 2, alpha = .6) +
      geom_tiplab(size = 2, fontface = "italic") +
      scale_x_continuous("Time (Ma)", labels = \(b) -b) +
      coord_cartesian(xlim = c(-700, 60)) +              # same scale on every page
      theme_tree2(base_size = 8) +
      ggtitle(sub("_sample\\.chronogram$", "", f))
  )
}
dev.off()
