# =============================================================
# Figure 4 Option 1 -- Dot + CrI plot, species x compound
# Consumes: depol_long, means_marginal (from EAD_bayes.R)
# =============================================================

library(dplyr)
library(ggplot2)
library(ggtext)

compound_order <- depol_long %>%
  filter(!is.na(depol)) %>%
  group_by(compound) %>%
  summarise(grand_mean = mean(depol, na.rm = TRUE), .groups = "drop") %>%
  arrange(grand_mean) %>%
  pull(compound)

compound_labels <- c(
  "caryophyllene" = "\u03b2-caryophyllene",
  "cocimene"      = "(<i>E</i>)-\u03b2-ocimene",
  "limonene"      = "\u03b4-limonene",
  "linalool"      = "linalool",
  "sixmethyl"     = "6-methyl-5-hepten-2-one",
  "tocimene"      = "(<i>Z</i>)-\u03b2-ocimene",
  "z3hex"         = "(<i>Z</i>)-3-hexenyl acetate"
)

species_colors <- c(
  "janthinus"      = "#D85A30",
  "janthiniformis" = "#0F6E56"
)
species_labels <- c(
  "janthinus"      = "italic('M. janthinus')",
  "janthiniformis" = "italic('M. janthiniformis')"
)

fig4_dot <- means_marginal %>%
  mutate(
    compound = factor(compound, levels = compound_order),
    species  = factor(species,  levels = c("janthinus", "janthiniformis"))
  ) %>%
  ggplot(aes(x = post_med, y = compound, color = species)) +
  geom_linerange(aes(xmin = cri_lo, xmax = cri_hi),
                 position = position_dodge(width = 0.5),
                 linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3) +
  scale_y_discrete(labels = compound_labels) +
  scale_color_manual(values = species_colors,
                     labels = lapply(species_labels, function(x) parse(text = x))) +
  labs(x = "Posterior median depolarization (mV)",
       y = NULL,
       color = NULL) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#e8e8e8", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    panel.border       = element_rect(color = "#d9d9d9", fill = NA),
    axis.text.y        = element_markdown(),
    legend.position    = "bottom",
    legend.text        = element_text(size = 10)
  )

fig4_dot

ggsave("figures/fig4_option1_dotplot.pdf", fig4_dot,
       width = 6, height = 4, device = cairo_pdf)

# =============================================================
# Figure 4 Option 2 -- Heatmap, species x compound
# Consumes: depol_long, means_marginal (from EAD_bayes.R)
# =============================================================

library(dplyr)
library(ggplot2)
library(ggtext)

compound_order <- depol_long %>%
  filter(!is.na(depol)) %>%
  group_by(compound) %>%
  summarise(grand_mean = mean(depol, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(grand_mean)) %>%  # descending so highest is at top
  pull(compound)

compound_labels <- c(
  "caryophyllene" = "\u03b2-caryophyllene",
  "cocimene"      = "(<i>E</i>)-\u03b2-ocimene",
  "limonene"      = "\u03b4-limonene",
  "linalool"      = "linalool",
  "sixmethyl"     = "6-methyl-5-hepten-2-one",
  "tocimene"      = "(<i>Z</i>)-\u03b2-ocimene",
  "z3hex"         = "(<i>Z</i>)-3-hexenyl acetate"
)

species_labels <- c(
  "janthinus"      = "italic('M. janthinus')",
  "janthiniformis" = "italic('M. janthiniformis')"
)

fig4_heat <- means_marginal %>%
  mutate(
    compound = factor(compound, levels = compound_order),
    species  = factor(species,  levels = c("janthinus", "janthiniformis"))
  ) %>%
  ggplot(aes(x = species, y = compound, fill = post_med)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", post_med, cri_lo, cri_hi)),
            size = 2.8, color = "white", lineheight = 1.1) +
  scale_fill_gradient(low = "#e8f4f0", high = "#0F6E56",
                      name = "Posterior median\ndepolarization (mV)") +
  scale_x_discrete(labels = lapply(species_labels, function(x) parse(text = x))) +
  scale_y_discrete(labels = compound_labels) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    axis.text.x      = element_text(size = 10, face = "italic"),
    axis.text.y      = element_markdown(size = 10),
    axis.ticks       = element_blank(),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 9)
  )

fig4_heat

ggsave("figures/fig4_option2_heatmap.pdf", fig4_heat,
       width = 6, height = 4, device = cairo_pdf)
