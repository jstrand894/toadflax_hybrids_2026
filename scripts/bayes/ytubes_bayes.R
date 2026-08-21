library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

# ── Data ──────────────────────────────────────────────────────────────────────
# New spreadsheet structure:
#   species   = weevil species (janthinus / janthiniformis)
#   Sex       = weevil sex (female / male) -- analyzed separately
#   Pos Trt   = plant genotype in the positive arm: yt6, dt6, hbr, rad, d6y6, y6d6
#               d6y6 = hybrid (dt6 parent first), y6d6 = hybrid (yt6 parent first)
#               Neg Trt is always blank (bl)
#   # Pos     = beetles choosing plant arm
#   # Neg     = beetles choosing blank arm
#   # NR      = non-responders (excluded from denominator)
#
# Each row is a single test day. Rows are aggregated across days within each
# species x sex x plant cell before testing.
#
# Male data exists only for dt6 and yt6; all other plant genotypes are female-only.

df_raw <- read_excel("data/ytubes/ytubes_results.xlsx",
                     sheet = "all") %>%
  rename(
    season    = Season,
    date      = Date,
    sex       = Sex,
    pos_trt   = 'Pos Trt',
    neg_trt   = 'Neg Trt',
    pos_count = '# Pos',
      neg_count = '# Neg',
      nr_count  = '# NR'
  ) %>%
  dplyr::select(species, season, date, sex, pos_trt, neg_trt,
                pos_count, neg_count, nr_count) %>%
  filter(!is.na(species)) %>%
  mutate(
    species = factor(species),
    sex     = factor(tolower(sex)),
    pos_trt = factor(pos_trt, levels = c("YT6", "DT6", "Y6D6", 
                                         "D6Y6", "RAD", "HBR"))
  )

# Aggregate across days within each species x sex x plant cell.
# NR excluded from total (denominator = choosers only).
df <- df_raw %>%
  group_by(species, sex, pos_trt) %>%
  summarise(
    pos_count = sum(pos_count, na.rm = TRUE),
    neg_count = sum(neg_count, na.rm = TRUE),
    nr_count  = sum(nr_count,  na.rm = TRUE),
    n_days    = n(),
    .groups   = "drop"
  ) %>%
  mutate(total = pos_count + neg_count)

# ── Bayesian Beta-Binomial per cell ──────────────────────────────────────────
# Beta(1, 1) prior (uniform) on choice proportion p.
# Posterior is Beta(pos_count + 1, neg_count + 1) analytically.
# Report: posterior median, 95% CrI, P(p > 0.5 | data).

bayes_df <- df %>%
  mutate(
    alpha_post = pos_count + 1,   # Beta(1,1) prior + successes
    beta_post  = neg_count + 1,   # Beta(1,1) prior + failures
    estimate   = qbeta(0.500, alpha_post, beta_post),
    lower.CL   = qbeta(0.025, alpha_post, beta_post),
    upper.CL   = qbeta(0.975, alpha_post, beta_post),
    p_gt_half  = pbeta(0.5, alpha_post, beta_post, lower.tail = FALSE)
  ) %>%
  mutate(
    sig = ifelse(p_gt_half >= 0.95, "*", "")
  )

# Print results table
bayes_df %>%
  dplyr::select(species, sex, pos_trt, pos_count, neg_count, nr_count,
                total, n_days, estimate, lower.CL, upper.CL, p_gt_half, sig) %>%
  arrange(species, sex, pos_trt) %>%
  print(n = Inf)