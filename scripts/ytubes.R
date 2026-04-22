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

# ── Per-cell binomial tests ───────────────────────────────────────────────────
# Two-sided exact binomial tests for CIs; one-sided (greater) for p-values.
# Tests whether proportion choosing plant arm exceeds 0.5.
# BH correction applied across all 16 simultaneous tests.

binom_df <- df %>%
  rowwise() %>%
  mutate(
    bt_two   = list(binom.test(pos_count, total, p = 0.5, alternative = "two.sided")),
    bt_one   = list(binom.test(pos_count, total, p = 0.5, alternative = "greater")),
    estimate = bt_two$estimate,
    lower.CL = bt_two$conf.int[1],
    upper.CL = bt_two$conf.int[2],
    p.value  = bt_one$p.value
  ) %>%
  ungroup() %>%
  dplyr::select(-bt_two, -bt_one) %>%
  mutate(
    sig = ifelse(p.value < 0.05, "*", "")
  )

# Print results table
binom_df %>%
  dplyr::select(species, sex, pos_trt, pos_count, neg_count, nr_count,
                total, n_days, estimate, lower.CL, upper.CL, p.value, sig) %>%
  arrange(species, sex, pos_trt) %>%
  print(n = Inf)


