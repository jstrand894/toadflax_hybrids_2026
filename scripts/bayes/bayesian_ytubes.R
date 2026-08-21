# =============================================================================
# Y-tube olfactometry -- Bayesian binomial analysis
# Beta-Binomial posterior for each species x sex x germplasm cell
# Prior: Beta(1, 1) -- uniform on [0,1], equivalent to one success and one
#        failure of prior information; weakly informative and standard for
#        proportion estimation with no strong prior expectation.
# Posterior: Beta(alpha + pos, beta + neg) analytically -- no MCMC needed.
# =============================================================================

library(dplyr)
library(tidyr)
library(readxl)

# ── Data ─────────────────────────────────────────────────────────────────────
df_raw <- read_excel("data/ytubes/ytubes_results.xlsx", sheet = "all") %>%
  rename(
    season    = Season,
    date      = Date,
    sex       = Sex,
    pos_trt   = `Pos Trt`,
    neg_trt   = `Neg Trt`,
    pos_count = `# Pos`,
    neg_count = `# Neg`,
    nr_count  = `# NR`
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

# Aggregate across days within each species x sex x germplasm cell.
# NR excluded from denominator.
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

# ── Bayesian Beta-Binomial posterior per cell ─────────────────────────────────
# Beta(1,1) prior + Binomial(pos, total) likelihood
# => Beta(1 + pos, 1 + neg) posterior
#
# Quantities reported:
#   post_med   -- posterior median proportion choosing plant arm
#   cri_lo     -- 2.5th percentile of posterior (lower 95% CrI)
#   cri_hi     -- 97.5th percentile of posterior (upper 95% CrI)
#   p_gt_half  -- P(proportion > 0.5 | data): posterior probability that
#                 weevils prefer the plant arm over chance

bayes_df <- df %>%
  rowwise() %>%
  mutate(
    alpha_post = 1 + pos_count,   # Beta posterior shape 1
    beta_post  = 1 + neg_count,   # Beta posterior shape 2
    post_med   = qbeta(0.50, alpha_post, beta_post),
    cri_lo     = qbeta(0.025, alpha_post, beta_post),
    cri_hi     = qbeta(0.975, alpha_post, beta_post),
    p_gt_half  = pbeta(0.5, alpha_post, beta_post, lower.tail = FALSE)
  ) %>%
  ungroup() %>%
  dplyr::select(-alpha_post, -beta_post)

# Print results
bayes_df %>%
  dplyr::select(species, sex, pos_trt, pos_count, neg_count, nr_count,
                total, n_days, post_med, cri_lo, cri_hi, p_gt_half) %>%
  arrange(species, sex, pos_trt) %>%
  print(n = Inf)