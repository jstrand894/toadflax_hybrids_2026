library(dplyr); library(ggplot2); library(forcats); library(patchwork); library(tidyr)

blend <- raw_df %>%
  select(Plant, Day, all_of(ead_cols)) %>%
  rowwise() %>%
  mutate(tot = sum(c_across(all_of(ead_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(tot > 0) %>%
  mutate(across(all_of(ead_cols), ~ .x / tot * 100)) %>%
  select(-tot) %>%
  pivot_longer(-c(Plant, Day), names_to = "compound", values_to = "pct")
lab <- c("cis-b-ocimene"               = "(Z)-\u03b2-ocimene",
         "trans-b-ociemene - 11.42"    = "(E)-\u03b2-ocimene",
         "z-3-hex acetate - 10.46"     = "(Z)-3-hexenyl acetate",
         "caryophyllene - 22.42"       = "\u03b2-caryophyllene",
         "linalool - 13.41"            = "linalool",
         "6-methyl-3-heptanone - 9.81" = "6-methyl-5-hepten-2-one",
         "limonene - 11.07"            = "limonene")
# ^ isomer labels here follow Table 1. Fix once you've settled the flip.

ord <- c("(Z)-\u03b2-ocimene", "(E)-\u03b2-ocimene", "(Z)-3-hexenyl acetate",
         "\u03b2-caryophyllene", "linalool", "6-methyl-5-hepten-2-one",
         "limonene")

plant_ord <- c("D6", "D6Y6", "Y6D6", "Y6", "HBR", "RAD")

pal <- c("(Z)-\u03b2-ocimene"        = "#2E5E4E",
         "(E)-\u03b2-ocimene"        = "#5C8D77",
         "(Z)-3-hexenyl acetate"     = "#8FBBA3",
         "\u03b2-caryophyllene"      = "#C6D9C8",
         "linalool"                  = "#E8D9A0",
         "6-methyl-5-hepten-2-one"   = "#D8A657",
         "limonene"                  = "#B5432F")

bl <- blend %>%
  mutate(cmpd = factor(unname(lab[compound]), levels = ord),
         Plant = factor(Plant, levels = plant_ord),
         grp = if_else(Plant %in% c("HBR", "RAD"), "Wild hybrid", "Synthetic"))

pA_dat <- bl %>%
  group_by(Plant, grp, cmpd) %>%
  summarise(pct = mean(pct, na.rm = TRUE), .groups = "drop")

pA <- ggplot(pA_dat, aes(Plant, pct, fill = cmpd)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
  facet_grid(~ grp, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "% of EAD-active blend", tag = "A") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 9),
        panel.spacing = unit(1.1, "lines"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

lim_dat <- bl %>% filter(cmpd == "limonene")
lim_mean <- lim_dat %>%
  group_by(Plant, grp) %>%
  summarise(m = mean(pct, na.rm = TRUE),
            se = sd(pct, na.rm = TRUE) / sqrt(sum(!is.na(pct))),
            .groups = "drop")

pB <- ggplot(lim_dat, aes(Plant, pct)) +
  geom_jitter(width = 0.13, height = 0, size = 1.5, alpha = 0.35,
              colour = "#B5432F") +
  geom_errorbar(data = lim_mean, aes(Plant, y = m, ymin = m - se, ymax = m + se),
                width = 0.14, linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(data = lim_mean, aes(Plant, m), size = 2.6,
             colour = "#B5432F", inherit.aes = FALSE) +
  facet_grid(~ grp, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = NULL, y = "Limonene, % of blend", tag = "B") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(1.1, "lines"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

fig <- pA / pB + plot_layout(heights = c(1.6, 1))
print(fig)

# ggsave("fig3_blend.png", fig, width = 7.2, height = 6.4, dpi = 400)
# ggsave("fig3_blend.pdf", fig, width = 7.2, height = 6.4)