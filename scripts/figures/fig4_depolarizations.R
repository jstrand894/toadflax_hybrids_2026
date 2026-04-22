library(dplyr)
library(readxl)
library(ggtext)


depolarizations <-
read_excel("data/EAG/depolarizations.xlsx", sheet = "compound response")

compound_order <- depolarizations %>%
  dplyr::select(species, plant, compound, run1:run5) %>%
  pivot_longer(run1:run5, names_to = "run", values_to = "value") %>%
  group_by(compound) %>%
  summarise(grand_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(grand_mean) %>%
  pull(compound)

facet_labels <- c(
  "janthinus"      = "italic('M. janthinus')",
  "janthiniformis" = "italic('M. janthiniformis')"
)

plant_labels <- c(
  "dt"   = "DT6",
  "d6y6" = "D6Y6",
  "yt"   = "YT6",
  "y6d6" = "Y6D6",
  "hbr"  = "HBR",
  "rad"  = "RAD"
)

compound_labels <- c(
  "caryophyllene" = "\u03b2-caryophyllene",
  "cocimene"      = "(<i>E</i>)-\u03b2-ocimene",
  "limonene"      = "\u03b4-limonene",
  "linalool"      = "linalool",
  "sixmethyl"     = "6-methyl-5-hepten-2-one",
  "tocimene"      = "(<i>Z</i>)-\u03b2-ocimene",
  "z3hex"         = "(<i>Z</i>)-3-hexenyl acetate"
)

plant_colors <- c(
  "yt"   = "#D85A30",
  "y6d6" = "#F0997B",
  "dt"   = "#0F6E56",
  "d6y6" = "#5DCAA5",
  "rad"  = "#BA7517",
  "hbr"  = "#EF9F27"
)

depolarizations_plot <-
depolarizations %>% 
  dplyr::select(species, plant, compound, run1:run5) %>% 
  pivot_longer(run1:run5, names_to = "run", values_to = "value") %>%
  group_by(species, plant, compound) %>%
  summarise(
    mean = mean(value),
    se   = sd(value) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  mutate(
    plant    = factor(plant, levels = c("dt", "d6y6", "yt", "y6d6", "hbr", "rad")),
    species  = factor(species, levels = c("janthinus", "janthiniformis")),
    compound = factor(compound, levels = compound_order)
  ) %>%
  ggplot() +
  geom_bar(aes(x = compound, y = mean, fill = plant),
           stat = "identity",
           color = "black", linewidth = 0.2,
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = compound, ymin = mean - se, ymax = mean + se, group = plant),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  facet_wrap(~species, labeller = labeller(species = as_labeller(facet_labels, label_parsed)), 
             ncol = 1)  +
  labs(y = "Mean depolarization (mV)",
       x = NULL,
       fill = "Toadflax\ngermplasm") +
  theme(
    panel.grid = element_blank(),       
    panel.background = element_blank(),  
    plot.background = element_blank(),
    legend.title       = element_text(size = 10, face = "bold"),
    panel.border = element_rect(color = "#d9d9d9", fill = NA),
    axis.text.x = element_markdown(angle = 25, hjust = 1)
  ) +
  scale_x_discrete(labels = compound_labels) +
  scale_fill_manual(values = plant_colors, labels = plant_labels)

depolarizations_plot

ggsave("figures/fig4_depolarizations.pdf", depolarizations_plot,
       width = 7, height = 4, device = cairo_pdf)


