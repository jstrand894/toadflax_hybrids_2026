# =============================================================
# EAD Depolarization Analysis -- Bayesian LMM
# Mecinus janthinus vs. M. janthiniformis compound response
# =============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(brms)
library(cmdstanr)

cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
options(mc.cores = parallel::detectCores())

# -------------------------------------------------------------
# 1. Load and reshape data
# -------------------------------------------------------------
depol_data <- read_excel("data/EAG/depolarizations.xlsx",
                         sheet = "compound response") %>%
  dplyr::select(species, plant, compound, run1:run5)

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

# -------------------------------------------------------------
# 2. Fit Bayesian LMM
# Fixed: species * compound + plant
# Random: weevil_id (individual antenna)
# -------------------------------------------------------------
fit_eag <- brm(
  depol ~ species * compound + plant + (1 | weevil_id),
  data    = depol_long,
  family  = gaussian(),
  prior   = c(
    prior(normal(0, 5),   class = b),
    prior(normal(0, 5),   class = Intercept),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains  = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.95),
  seed    = 42, backend = "cmdstanr"
)

print(fit_eag)  # check Rhat, ESS

# -------------------------------------------------------------
# 3. Posterior means per species x plant x compound
# -------------------------------------------------------------
nd <- depol_long %>%
  distinct(species, plant, compound) %>%
  mutate(weevil_id = NA)

epred_eag <- posterior_epred(fit_eag, newdata = nd, re_formula = NA)
# epred_eag is [draws x nrow(nd)]

means_df <- nd %>%
  mutate(
    post_med = apply(epred_eag, 2, median),
    cri_lo   = apply(epred_eag, 2, quantile, 0.025),
    cri_hi   = apply(epred_eag, 2, quantile, 0.975)
  )

print(means_df, n = Inf)

# -------------------------------------------------------------
# 3b. Posterior means per species x compound (marginalized over plant)
# -------------------------------------------------------------
nd_marginal <- depol_long %>%
  distinct(species, compound)

# Average posterior draws across all 6 plant levels for each species x compound
means_marginal <- map_dfr(levels(depol_long$compound), function(cmpd) {
  map_dfr(levels(depol_long$species), function(spp) {
    idx   <- which(nd$compound == cmpd & nd$species == spp)
    draws <- rowMeans(epred_eag[, idx, drop = FALSE])
    tibble(
      species  = spp,
      compound = cmpd,
      post_med = median(draws),
      cri_lo   = as.numeric(quantile(draws, 0.025)),
      cri_hi   = as.numeric(quantile(draws, 0.975))
    )
  })
})

# -------------------------------------------------------------
# 4. Pairwise species contrasts within each compound
#    Marginalized over plant by averaging contrasts across
#    all 6 plant levels
# -------------------------------------------------------------
compounds <- levels(depol_long$compound)
plants    <- levels(depol_long$plant)

contrasts_df <- map_dfr(compounds, function(cmpd) {
  
  # For each plant, get the draws for each species at this compound
  # then average the contrast across plants
  diff_by_plant <- map(plants, function(plt) {
    i_j  <- which(nd$compound == cmpd & nd$species == "janthinus"      & nd$plant == plt)
    i_jf <- which(nd$compound == cmpd & nd$species == "janthiniformis" & nd$plant == plt)
    epred_eag[, i_j] - epred_eag[, i_jf]
  })
  
  # Average the per-plant draw vectors (each is length n_draws)
  diff <- Reduce("+", diff_by_plant) / length(plants)
  
  tibble(
    compound = cmpd,
    contrast = "janthinus - janthiniformis",
    diff_med = median(diff),
    cri_lo   = quantile(diff, 0.025),
    cri_hi   = quantile(diff, 0.975),
    pd       = max(mean(diff > 0), mean(diff < 0))
  )
})

print(contrasts_df, n = Inf)
