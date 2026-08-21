# =============================================================================
# Figure 7 -- Chemistry and behavior as one relationship, not two lists
# One point per germplasm x species. X = compositional position (PC1).
# Y = posterior choice proportion. A segment connects the two species at
# each germplasm, so the species gap is a visible distance rather than
# something the reader has to compute by eye across two panels.
#
# Consumes: pca_out, epred, plants (from 01_pca.R)
#           bayes_df (from ytubes_bayes.R)
# =============================================================================

source("scripts/bayes/01_pca.R")
source("scripts/bayes/ytubes_bayes.R")
stopifnot(exists("pca_out"), exists("epred"), exists("plants"), exists("bayes_df"))

library(tidyverse)
library(ggrepel)

# ── 1. Shared germplasm axis ─────────────────────────────────────────────────
germplasm_order <- c("D6", "D6Y6", "Y6", "Y6D6", "RAD", "HBR")

voc_to_ytube <- c(D6 = "DT6", Y6 = "YT6", D6Y6 = "D6Y6",
                  Y6D6 = "Y6D6", RAD = "RAD", HBR = "HBR")

plant_short <- c(
  "D6" = "D6", "Y6" = "Y6", "D6Y6" = "D6Y6",
  "Y6D6" = "Y6D6", "RAD" = "RAD", "HBR" = "HBR"
)

# ── 2. Compositional position -- posterior PC1 per germplasm ─────────────────
pc1_draws <- sapply(seq_along(plants), function(i) {
  centered <- sweep(epred[, i, ], 2, pca_out$center)
  as.vector(centered %*% pca_out$rotation[, 1])
})
colnames(pc1_draws) <- plants

pc1_summary <- as_tibble(pc1_draws) %>%
  pivot_longer(everything(), names_to = "Plant", values_to = "pc1") %>%
  group_by(Plant) %>%
  summarise(pc1_med = median(pc1),
            pc1_lo  = quantile(pc1, .025),
            pc1_hi  = quantile(pc1, .975),
            .groups = "drop")

# ── 3. Behavioral choice -- female data, all six germplasm ───────────────────
choice_summary <- bayes_df %>%
  filter(sex == "female") %>%
  mutate(Plant = names(voc_to_ytube)[match(as.character(pos_trt), voc_to_ytube)]) %>%
  select(species, Plant, estimate, lower.CL, upper.CL)

# ── 4. Join chemistry and behavior on germplasm ───────────────────────────────
plot_df <- choice_summary %>%
  left_join(pc1_summary, by = "Plant") %>%
  mutate(Plant = factor(Plant, levels = germplasm_order))

# One chemistry interval per germplasm; retain true PC1 position
chem_err <- pc1_summary %>%
  mutate(Plant = factor(Plant, levels = germplasm_order))

# Separate the two species slightly along x, reproducibly
pc1_span <- diff(range(chem_err$pc1_med, na.rm = TRUE))
x_dodge <- max(0.012 * pc1_span, 0.01)

species_dodge <- c(
  "janthinus" = -x_dodge,
  "janthiniformis" = x_dodge
)

plot_df <- plot_df %>%
  mutate(
    x_order = pc1_med,
    x_plot = pc1_med + unname(species_dodge[as.character(species)])
  )

segment_df <- plot_df %>%
  select(Plant, species, x_plot, estimate) %>%
  pivot_wider(
    names_from = species,
    values_from = c(x_plot, estimate),
    names_sep = "_"
  )

chem_err <- chem_err %>%
  mutate(x_order = pc1_med)


# ── 6. Plot ──────────────────────────────────────────────────────────────────
fig7 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_segment(
    data = segment_df,
    aes(
      x = x_plot_janthinus, xend = x_plot_janthiniformis,
      y = estimate_janthinus, yend = estimate_janthiniformis
    ),
    color = "grey65", linewidth = 0.5, linetype = "dotted"
  ) +
  geom_errorbar(
    data = plot_df,
    aes(x = x_plot, ymin = lower.CL, ymax = upper.CL, color = species),
    width = 0.08, linewidth = 0.6, alpha = 0.7
  ) +
  geom_point(
    data = plot_df,
    aes(x = x_plot, y = estimate, color = species),
    size = 3.4
  ) +
  geom_text(
    data = chem_err,
    aes(x = pc1_med, y = 0.965, label = plant_short[as.character(Plant)]),
    size = 3.2,
    color = "grey20",
    vjust = 0
  ) +
  scale_color_manual(values = species_colors,
                     labels = c("janthinus" = "M. janthinus",
                                "janthiniformis" = "M. janthiniformis"),
                     name = NULL) +
  scale_y_continuous(limits = c(0.25, 1.0), breaks = seq(0.3, 1.0, 0.1)) +
  annotate("text", x = min(chem_err$pc1_lo), y = 0.225,
           label = "\u2190 wild hybrid composition", hjust = 0, size = 3, color = "grey40") +
  annotate("text", x = max(chem_err$pc1_hi), y = 0.225,
           label = "parental composition \u2192", hjust = 1, size = 3, color = "grey40") +
  labs(
    x = "Compositional position (PC1, posterior median with 95% CrI)",
    y = "Proportion choosing plant arm\n(posterior median, 95% CrI)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(face = "italic", size = 9),
    axis.title       = element_text(size = 10),
    plot.margin      = margin(18, 10, 6, 6)
  )

fig7

# ggsave("figures/fig7_chemistry_behavior.pdf", fig7,
       # width = 7.5, height = 5, device = cairo_pdf)
