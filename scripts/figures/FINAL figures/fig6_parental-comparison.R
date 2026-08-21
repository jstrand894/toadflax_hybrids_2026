# =============================================================================
# Figure 6 -- Posterior pairwise Aitchison centroid distances (all 15 pairs)
# Forest plot, ordered smallest to largest, with the D6-Y6 parental
# separation drawn as a reference line.
#
# Consumes: epred [draws, 6, 20] and `plants` from bayesian_vocs.R
# =============================================================================

source("scripts/bayes/bayesian_vocs.R")
stopifnot(exists("epred"), exists("plants"))

library(tidyverse)
library(ggtext)

# ── 1. Full posterior draws for every pair (not just summaries) ───────────────
pair_list <- combn(plants, 2, simplify = FALSE)

pair_draws <- map(pair_list, function(p) {
  i <- match(p[1], plants)
  j <- match(p[2], plants)
  sqrt(rowSums((epred[, i, ] - epred[, j, ])^2))
})
names(pair_draws) <- map_chr(pair_list, ~ paste(.x, collapse = " vs "))

# Parental reference distance, kept as full draws so the ratio carries
# uncertainty from both numerator and denominator
i0 <- match("D6", plants)
j0 <- match("Y6", plants)
d_parents <- sqrt(rowSums((epred[, i0, ] - epred[, j0, ])^2))

# ── 2. Summarise, including ratio to parental separation ─────────────────────
hybrids   <- c("D6Y6", "Y6D6")
wild      <- c("RAD", "HBR")
parentals <- c("D6", "Y6")

classify <- function(a, b) {
  s <- c(a, b)
  if (setequal(s, parentals))            return("Parental separation")
  if (setequal(s, wild))                 return("Convergence pair")
  if (setequal(s, hybrids))              return("Convergence pair")
  if (any(s %in% wild) & any(s %in% parentals)) return("Wild hybrid vs parental")
  "Other pair"
}

dist_summary <- imap_dfr(pair_draws, function(d, nm) {
  p     <- str_split(nm, " vs ", simplify = TRUE)
  ratio <- d / d_parents
  tibble(
    pair      = nm,
    group1    = p[1],
    group2    = p[2],
    dist_med  = median(d),
    dist_lo   = unname(quantile(d, 0.025)),
    dist_hi   = unname(quantile(d, 0.975)),
    ratio_med = median(ratio),
    ratio_lo  = unname(quantile(ratio, 0.025)),
    ratio_hi  = unname(quantile(ratio, 0.975)),
    p_less    = mean(d < d_parents)
  )
}) %>%
  mutate(
    category = map2_chr(group1, group2, classify),
    category = factor(category, levels = c(
      "Convergence pair", "Wild hybrid vs parental",
      "Parental separation", "Other pair")),
    pair = fct_reorder(pair, dist_med)
  )

print(dist_summary, n = Inf, width = Inf)

# ── 3. Reference line values from the parental posterior ─────────────────────
ref_med <- median(d_parents)
ref_lo  <- unname(quantile(d_parents, 0.025))
ref_hi  <- unname(quantile(d_parents, 0.975))

# ── 4. Ratio annotation for the two convergence pairs ────────────────────────
annot <- dist_summary %>%
  filter(category == "Convergence pair") %>%
  mutate(lab = sprintf("%.0f%% of parental separation", 100 * ratio_med))

# ── 5. Plot ──────────────────────────────────────────────────────────────────
pal <- c(
  "Convergence pair"        = "#0F6E56",
  "Wild hybrid vs parental" = "#D85A30",
  "Parental separation"     = "#3B3B3B",
  "Other pair"              = "#9E9E9E"
)

x_max <- max(dist_summary$dist_hi) * 1.32

fig6 <- ggplot(dist_summary, aes(x = dist_med, y = pair, color = category)) +
  annotate("rect", xmin = ref_lo, xmax = ref_hi, ymin = -Inf, ymax = Inf,
           fill = "#3B3B3B", alpha = 0.07) +
  geom_vline(xintercept = ref_med, linetype = "dashed",
             color = "#3B3B3B", linewidth = 0.5) +
  geom_linerange(aes(xmin = dist_lo, xmax = dist_hi), linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(data = annot, aes(x = dist_hi, label = lab),
            hjust = -0.08, size = 2.8, show.legend = FALSE) +
  annotate("text", x = ref_med, y = Inf, vjust = -0.6, hjust = 0.5,
           label = "D6 vs Y6 parental separation",
           size = 2.9, color = "#3B3B3B") +
  scale_color_manual(values = pal, name = NULL) +
  scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Aitchison centroid distance\n[posterior median, 95% CrI]",
    y = NULL
  ) +
  theme(
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    panel.border       = element_rect(color = "#d9d9d9", fill = NA),
    panel.grid.major.y = element_line(color = "#ececec", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks         = element_blank(),
    axis.text.y        = element_text(size = 9),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    legend.key         = element_blank(),
    plot.margin        = margin(14, 10, 6, 6)
  )

fig6

# ggsave("figures/fig6_centroid_distances.pdf", fig6,
#        width = 7, height = 5.2, device = cairo_pdf)