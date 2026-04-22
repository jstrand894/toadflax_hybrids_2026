source("scripts/ytubes.R")


# ── Plot ──────────────────────────────────────────────────────────────────────
# Faceted by sex; x-axis = plant genotype; color = weevil species.
# Males only have dt6 and yt6 data.

p1 <- ggplot(binom_df %>% filter(sex == "female"),
             aes(x = pos_trt, y = estimate, color = species)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.15) +
  geom_point(size = 3) +
  geom_text(aes(label = sig, y = upper.CL + 0.02),
            size = 5,
            show.legend = FALSE) +
  scale_x_discrete(limits = c("DT6", "D6Y6", "YT6", "Y6D6", "RAD", "HBR")) +
  scale_y_continuous(limits = c(0.25, 1.00),
                     labels = scales::number_format(accuracy = 0.01)) +
  facet_wrap(~ species, ncol = 2,
             labeller = labeller(species = c(
               "janthiniformis" = "M. janthiniformis",
               "janthinus"      = "M. janthinus"))) +
  labs(
    x       = "\nToadflax germplasm",
    y       = "Proportion choosing plant arm",
  ) +
  theme_classic() +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "italic"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position   = "none",
    axis.line         = element_blank()
  ) +
  scale_color_manual(values = c(
    "janthinus"      = "#D85A30",
    "janthiniformis" = "#0F6E56"
  ))

p1


ggsave("figures/fig5_ytubes.pdf", p1,
       width = 7, height = 4, device = cairo_pdf)
