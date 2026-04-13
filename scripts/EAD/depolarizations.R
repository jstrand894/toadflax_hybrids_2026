# =============================================================
# EAG Depolarization Analysis
# Mecinus janthinus vs. M. janthiniformis compound response
# Linear mixed model with pairwise emmeans contrasts
# =============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggsignif)

# -------------------------------------------------------------
# 1. Load data
# -------------------------------------------------------------

depol_data <-
  read_excel("data/EAG/depolarizations.xlsx",
             sheet = "compound response") %>%
  dplyr::select(species, plant, compound, run1:run5)

# -------------------------------------------------------------
# 2. Clean and reshape
# Each row is one species x plant x compound combination.
# run1:run5 are individual depolarization recordings (mV above
# baseline); pivot to long format for modeling.
# -------------------------------------------------------------

depol_long <- depol_data %>%
  filter(!is.na(species), !is.na(compound),
         species != "test", compound != "test") %>%
  pivot_longer(cols      = starts_with("run"),
               names_to  = "run",
               values_to = "depol") %>%
  mutate(species   = as.factor(species),
         compound  = as.factor(compound),
         plant     = as.factor(plant),
         weevil_id = as.factor(paste(species, plant, run, sep = "_")))
depol_long %>%
  distinct(weevil_id) %>%
  nrow()
depol_long %>%
  distinct(species, plant, run) %>%
  count(species, plant)
# -------------------------------------------------------------
# 3. Fit linear mixed model
# Fixed effects: species * compound (full factorial)
# Random effect: plant (individual replicate, not a treatment)
# -------------------------------------------------------------

mod <- lmer(depol ~ species * compound + (1 | weevil_id), 
            data = depol_long)


summary(mod)
VarCorr(mod)

# Type III F-tests via Satterthwaite approximation
anova(mod)

# -------------------------------------------------------------
# 4. Pairwise emmeans contrasts
# Compare species within each compound; Holm adjustment for
# 7 simultaneous tests (one per compound)
# -------------------------------------------------------------

emm <- emmeans(mod, pairwise ~ species | compound)

summary(emm$contrasts, infer = TRUE, adjust = "holm")

# -------------------------------------------------------------
# 5. Build annotation data frame for significance brackets
# -------------------------------------------------------------

emm_means      <- as.data.frame(emm$emmeans)
compound_levels <- levels(emm_means$compound)

annot_df <- summary(emm$contrasts, infer = TRUE, adjust = "holm") %>%
  as.data.frame() %>%
  mutate(
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  filter(sig_label != "ns") %>%
  left_join(
    emm_means %>%
      group_by(compound) %>%
      summarise(y_pos = max(upper.CL) + 0.3, .groups = "drop"),
    by = "compound"
  ) %>%
  mutate(x_num = as.numeric(factor(compound, levels = compound_levels)))

# -------------------------------------------------------------
# 6. Plot estimated marginal means with significance brackets
# -------------------------------------------------------------

ggplot(emm_means,
       aes(x = compound, y = emmean, color = species, group = species)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(0.4), width = 0.2) +
  geom_signif(
    xmin        = annot_df$x_num - 0.15,
    xmax        = annot_df$x_num + 0.15,
    annotations = annot_df$sig_label,
    y_position  = annot_df$y_pos,
    tip_length  = 0.01,
    textsize    = 4,
    color       = "black"
  ) +
  labs(y = "Estimated mean depolarization (mV)",
       x = "Compound") +
  theme_classic()