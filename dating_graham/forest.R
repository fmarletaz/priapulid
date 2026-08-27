## Landscape A4 forest plots: one panel per node, one row per gene family.
## Ports forest_panels.py. Reads either table:
##   og  age  lo  hi  label  type  members1  members2   (readdiv post-processing)
##   og  clade  age  lo  hi  n_species                  (chronograms.py, has header)

library(dplyr); library(tidyr); library(ggplot2)

CHENGJIANG <- 518
CAMBRIAN   <- c(538.8, 485.4)
NODES      <- c("dupl:Priapulida/Priapulida", "Priapulida", "Ecdysozoa")
LABELS     <- c("dupl:Priapulida/Priapulida" = "Priapulid WGD",
                "Priapulida"                 = "Priapulida",
                "Ecdysozoa"                  = "Ecdysozoa")

## ---- read ------------------------------------------------------------------
raw <- read.delim("dat_res_clg_3.txt", header = FALSE, stringsAsFactors = FALSE)

if (ncol(raw) >= 6 && suppressWarnings(!is.na(as.numeric(raw[1, 2])))) {
  d <- raw %>%
    setNames(c("og","age","lo","hi","label","type","m1","m2")[seq_len(ncol(raw))]) %>%
    mutate(node = ifelse(type == "dupl", paste0("dupl:", label), label))
} else {                                   # chronograms.py output, with header
  d <- read.delim("chrono_nodes.tsv") %>%
    rename(node = clade)
}

## one value per family per node (a few families contribute >1 duplication)
d <- d %>%
  filter(node %in% NODES) %>%
  group_by(node, og) %>%
  summarise(age = median(age), lo = median(lo), hi = median(hi), .groups = "drop")

## ---- pooled random-effects estimate (DerSimonian-Laird) --------------------
pool_re <- function(y, lo, hi) {
  v  <- ((hi - lo) / 3.92)^2
  w  <- 1 / v
  fx <- sum(w * y) / sum(w)
  Q  <- sum(w * (y - fx)^2); df <- length(y) - 1
  C  <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)
  w2 <- 1 / (v + tau2)
  mu <- sum(w2 * y) / sum(w2); sem <- sqrt(1 / sum(w2))
  tibble(mu = mu, ci_lo = mu - 1.96 * sem, ci_hi = mu + 1.96 * sem,
         tau = sqrt(tau2), I2 = if (Q > 0) max(0, (Q - df) / Q) * 100 else 0)
}

pooled <- d %>% group_by(node) %>% group_modify(~ pool_re(.x$age, .x$lo, .x$hi)) %>%
  ungroup()

n_fam <- d %>% count(node, name = "n_fam")

## ---- rank families within each panel, and build facet titles ---------------
d <- d %>% group_by(node) %>% arrange(age, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>% ungroup()

titles <- pooled %>% left_join(n_fam, by = "node") %>%
  mutate(facet = sprintf("%s\nn = %d families   I?? = %.0f%%   ?? = %.0f Myr",
                         LABELS[node], n_fam, I2, tau))
key <- setNames(titles$facet, titles$node)

d      <- d      %>% mutate(facet = factor(key[node], levels = key[NODES]))
pooled <- pooled %>% left_join(n_fam, by = "node") %>%
  mutate(facet = factor(key[node], levels = key[NODES]),
         y = -0.06 * n_fam)

## ---- plot ------------------------------------------------------------------
p <- ggplot(d, aes(y = rank)) +
  annotate("rect", xmin = CAMBRIAN[2], xmax = CAMBRIAN[1],
           ymin = -Inf, ymax = Inf, fill = "#f2e6c9") +
  geom_vline(xintercept = CHENGJIANG, colour = "#c1272d",
             linetype = "dashed", linewidth = .35) +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = facet),
                 linewidth = .22, alpha = .75) +
  geom_point(aes(x = age), size = .25, colour = "#14283f") +
  ## pooled estimate below the families
  geom_linerange(data = pooled, aes(y = y, xmin = ci_lo, xmax = ci_hi),
                 linewidth = .9, colour = "black", inherit.aes = FALSE) +
  geom_point(data = pooled, aes(x = mu, y = y), shape = 23, size = 2,
             fill = "black", inherit.aes = FALSE) +
  geom_text(data = pooled,
            aes(x = mu, y = y * 2.1,
                label = sprintf("%.0f Ma (%.0f???%.0f)", mu, ci_lo, ci_hi)),
            size = 2.6, vjust = 1, inherit.aes = FALSE) +
  facet_wrap(~ facet, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = setNames(
    c("#C1444A", "#5E82A8", "#5E82A8")[seq_along(NODES)], key[NODES]),
    guide = "none") +
  scale_x_reverse("Age (Ma)", breaks = seq(300, 800, 50)) +
  labs(y = "Gene families, ranked by median age") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_blank(),
        strip.text         = element_text(hjust = 0, size = 8.5, lineheight = 1.2))

ggsave("forest_panels.pdf", p, width = 11.69, height = 8.27, useDingbats = FALSE)

print(left_join(pooled, n_fam, by = "node") %>%
        transmute(node = LABELS[node], n_fam, pooled = round(mu, 1),
                  ci = sprintf("%.1f-%.1f", ci_lo, ci_hi),
                  tau = round(tau, 1), I2 = round(I2)))
