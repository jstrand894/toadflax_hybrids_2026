source("scripts/01_volatile_proc.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)    # column_to_rownames

# ── Shared color palette ───────────────────────────────────────────────────────
plant_colors <- c(
  "YT6"  = "#D85A30",
  "Y6D6" = "#F0997B",
  "DT6"  = "#0F6E56",
  "D6Y6" = "#5DCAA5",
  "RAD"  = "#BA7517",
  "HBR"  = "#EF9F27"
)

plant_levels <- c("DT6", "D6Y6", "YT6", "Y6D6", "RAD", "HBR")

# ── Top compounds to display ───────────────────────────────────────────────────
# Selected by EAD activity and/or high abundance across germplasm
top_compounds <- c(
  "cis-b-ocimene",
  "z-3-hex acetate - 10.46",
  "caryophyllene - 22.42",
  "decanal - 16.55",
  "nonanal - 13.51",
  "trans-b-ociemene - 11.42",
  "linalool - 13.41",
  "beta-pinene - 9.91",
  "dodecanal - 22.07"
)

# Display labels (with proper unicode characters for axis text)
compound_labels <- c(
  "cis-b-ocimene"            = "(<i>Z</i>)-\u03b2-ocimene",
  "z-3-hex acetate - 10.46"  = "(<i>Z</i>)-3-hexenyl acetate",
  "caryophyllene - 22.42"    = "\u03b2-caryophyllene",
  "decanal - 16.55"          = "decanal",
  "nonanal - 13.51"          = "nonanal",
  "trans-b-ociemene - 11.42" = "(<i>E</i>)-\u03b2-ocimene",
  "linalool - 13.41"         = "linalool",
  "beta-pinene - 9.91"       = "\u03b2-pinene",
  "dodecanal - 22.07"        = "dodecanal"
)

# ── Figure 2A: Mean emission rates (bar chart) ─────────────────────────────────
# Summarise raw concentrations to per-germplasm means ± SE
bar_df <- compounds_mat %>%
  as_tibble() %>%
  bind_cols(meta %>% select(Plant), .) %>%
  pivot_longer(-Plant, names_to = "compound", values_to = "value") %>%
  filter(compound %in% top_compounds) %>%
  group_by(Plant, compound) %>%
  summarise(
    mean_conc = mean(value, na.rm = TRUE),
    se        = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  ) %>%
  mutate(
    Plant    = recode(Plant, "D6" = "DT6", "Y6" = "YT6"),
    compound = factor(compound, levels = top_compounds),
    Plant    = factor(Plant, levels = plant_levels)
  )

mean_volatiles <- 
  ggplot(bar_df, aes(x = compound, y = mean_conc, fill = Plant)) +
  geom_bar(position = position_dodge(width = 0.75), 
           width = 0.7,
           color = "black",
           linewidth = 0.2,
           stat = "identity",
           alpha = 0.75) +
  geom_errorbar(
    aes(ymin = mean_conc - se, ymax = mean_conc + se),
    position = position_dodge(width = 0.75),
    width    = 0.2,
    linewidth = 0.4
  ) +
  scale_fill_manual(values = plant_colors, name = "Toadflax\ngermplasm") +
  scale_x_discrete(labels = compound_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = NULL,
    y = expression("Volatile concentration (ng g"^{-1} ~ "hr"^{-1} * ")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x        = element_markdown(angle = 35, hjust = 1),
    axis.line          = element_line(linewidth = 0.4),
    axis.ticks         = element_line(linewidth = 0.4),
    legend.position    = "right",
    legend.key.size    = unit(0.45, "cm"),
    legend.text        = element_text(size = 10),
    legend.title       = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

mean_volatiles

ggsave("figures/fig1_mean-volatiles.pdf", mean_volatiles,
       width = 7, height = 4, device = cairo_pdf)


