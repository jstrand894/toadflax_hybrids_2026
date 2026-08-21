# =============================================================
# results_summary.R
# Pull and print key results from EAD, Y-tube, and VOC analyses
# =============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(brms)
library(cmdstanr)
library(tidyverse)

cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
options(mc.cores = parallel::detectCores())


# =============================================================
# SECTION 1: EAD Depolarization (Bayesian LMM)
# =============================================================
cat("\n", strrep("=", 60), "\n")
cat("SECTION 1: EAD Depolarization Results\n")
cat(strrep("=", 60), "\n\n")

# -- Load or re-fit model --
if (file.exists("models/fit_eag.rds")) {
  fit_eag <- readRDS("models/fit_eag.rds")
  cat("Loaded fit_eag from models/fit_eag.rds\n\n")
} else {
  cat("fit_eag.rds not found -- re-fitting model...\n")
  source("scripts/bayes/EAD_bayes.R")
  dir.create("models", showWarnings = FALSE)
  saveRDS(fit_eag, "models/fit_eag.rds")
}

# -- Model summary --
cat("--- Model Summary (check Rhat, ESS) ---\n")
print(fit_eag)

# -- Reconstruct nd and epred_eag --
depol_data <- read_excel("data/EAG/depolarizations.xlsx",
                         sheet = "compound response") %>%
  dplyr::select(species, plant, compound, run1:run5)

depol_long <- depol_data %>%
  filter(!is.na(species), !is.na(compound),
         species != "test", compound != "test") %>%
  pivot_longer(cols = starts_with("run"),
               names_to = "run", values_to = "depol") %>%
  mutate(species   = as.factor(species),
         compound  = as.factor(compound),
         plant     = as.factor(plant),
         weevil_id = as.factor(paste(species, plant, run, sep = "_")))

nd <- depol_long %>%
  distinct(species, plant, compound) %>%
  mutate(weevil_id = NA)

epred_eag <- posterior_epred(fit_eag, newdata = nd, re_formula = NA)

# -- Posterior means per species x plant x compound --
cat("\n--- Posterior Medians: Species x Plant x Compound ---\n")
means_df <- nd %>%
  mutate(
    post_med = apply(epred_eag, 2, median),
    cri_lo   = apply(epred_eag, 2, quantile, 0.025),
    cri_hi   = apply(epred_eag, 2, quantile, 0.975)
  )
print(means_df, n = Inf)

# -- Posterior means marginalized over plant --
cat("\n--- Posterior Medians: Species x Compound (marginalized over plant) ---\n")
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
print(means_marginal, n = Inf)

# -- Pairwise species contrasts within each compound --
cat("\n--- Pairwise Contrasts: janthinus - janthiniformis (per compound) ---\n")
compounds <- levels(depol_long$compound)
plants    <- levels(depol_long$plant)

contrasts_df <- map_dfr(compounds, function(cmpd) {
  diff_by_plant <- map(plants, function(plt) {
    i_j  <- which(nd$compound == cmpd & nd$species == "janthinus"      & nd$plant == plt)
    i_jf <- which(nd$compound == cmpd & nd$species == "janthiniformis" & nd$plant == plt)
    epred_eag[, i_j] - epred_eag[, i_jf]
  })
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


# =============================================================
# SECTION 2: Y-tube Olfactometry (Bayesian GLMM)
# =============================================================
cat("\n", strrep("=", 60), "\n")
cat("SECTION 2: Y-tube Olfactometry Results\n")
cat(strrep("=", 60), "\n\n")

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
    pos_trt = factor(pos_trt, levels = c("YT6", "DT6", "Y6D6", "D6Y6", "RAD", "HBR"))
  )

# -- Raw descriptive counts per cell (no model, no significance claims) --
df <- df_raw %>%
  group_by(species, sex, pos_trt) %>%
  summarise(
    pos_count = sum(pos_count, na.rm = TRUE),
    neg_count = sum(neg_count, na.rm = TRUE),
    nr_count  = sum(nr_count,  na.rm = TRUE),
    n_days    = n(),
    .groups   = "drop"
  ) %>%
  mutate(total = pos_count + neg_count,
         raw_prop = pos_count / total)

cat("--- Raw Choice Counts per Cell (descriptive only) ---\n")
df %>%
  arrange(species, sex, pos_trt) %>%
  print(n = Inf)

# --- Female model: full germplasm panel ---
df_female <- df_raw %>%
  filter(sex == "female") %>%
  mutate(total = pos_count + neg_count,
         date  = factor(as.character(date))) %>%
  filter(total > 0) %>%
  droplevels()

# -- Design confound check (date x species coverage) --
cat("\n--- Design: Trials per Species x Germplasm ---\n")
print(df_female %>% count(species, pos_trt, wt = total), n = Inf)

cat("\n--- Design: Date x Species Coverage ---\n")
print(xtabs(~ date + species, df_female))

# -- Load or re-fit model --
if (file.exists("models/fit_female.rds")) {
  fit_female <- readRDS("models/fit_female.rds")
  cat("Loaded fit_female from models/fit_female.rds\n\n")
} else {
  cat("fit_female.rds not found -- re-fitting model...\n")
  fit_female <- brm(
    pos_count | trials(total) ~ species * pos_trt + (1 | date),
    data = df_female,
    family = binomial(link = "logit"),
    prior = set_prior("normal(0, 2)", class = "b"),
    chains = 4, cores = 4, iter = 4000, seed = 1,
    backend = "cmdstanr"
  )
  dir.create("models", showWarnings = FALSE)
  saveRDS(fit_female, "models/fit_female.rds")
}

if (nrow(fit_female$data) != nrow(df_female)) {
  stop("fit_female cache is stale relative to df_female -- delete models/fit_female.rds")
}

# -- Wide-prior sensitivity fit (same data, normal(0, 5) on fixed effects) --
if (file.exists("models/fit_female_prior_n5.rds")) {
  fit_female_wide <- readRDS("models/fit_female_prior_n5.rds")
  cat("Loaded fit_female_wide from models/fit_female_prior_n5.rds\n\n")
} else {
  cat("fit_female_prior_n5.rds not found -- re-fitting model...\n")
  fit_female_wide <- update(fit_female,
                             prior = set_prior("normal(0, 5)", class = "b"),
                             seed = 1)
  dir.create("models", showWarnings = FALSE)
  saveRDS(fit_female_wide, "models/fit_female_prior_n5.rds")
}

# Males were tested on a reduced germplasm panel and are summarized
# descriptively only (see the raw counts table above) -- not modeled here.

cat("--- Model Summary (check Rhat, ESS) ---\n")
print(summary(fit_female))

# -- Shared helpers: cell estimates, contrasts, and interaction contrasts
#    all read off the same posterior draws so nothing can drift apart --
diffs <- function(fit, dat) {
  nd <- dat %>% distinct(species, pos_trt) %>%
    mutate(total = 1, date = levels(dat$date)[1])
  p <- posterior_epred(fit, newdata = nd, re_formula = NA)
  map(levels(dat$pos_trt), function(g) {
    p[, which(nd$pos_trt == g & nd$species == "janthinus")] -
      p[, which(nd$pos_trt == g & nd$species == "janthiniformis")]
  }) %>% set_names(levels(dat$pos_trt))
}

cell_table <- function(fit, dat) {
  nd <- dat %>% distinct(species, pos_trt) %>%
    mutate(total = 1, date = levels(dat$date)[1])
  p <- posterior_epred(fit, newdata = nd, re_formula = NA)
  nd %>%
    select(species, pos_trt) %>%
    mutate(p_med     = apply(p, 2, median),
           lo        = apply(p, 2, quantile, .025),
           hi        = apply(p, 2, quantile, .975),
           p_gt_half = apply(p, 2, function(x) mean(x > 0.5))) %>%
    arrange(species, pos_trt)
}

contrast_table <- function(fit, dat) {
  dg <- diffs(fit, dat)
  imap_dfr(dg, function(d, g) {
    tibble(pos_trt  = g,
           contrast = "janthinus - janthiniformis",
           diff_med = median(d),
           lo       = quantile(d, .025),
           hi       = quantile(d, .975),
           pd       = max(mean(d > 0), mean(d < 0)))
  })
}

interaction_table <- function(fit, dat) {
  dg <- diffs(fit, dat)
  combn(names(dg), 2, simplify = FALSE) %>%
    map_dfr(function(pp) {
      dd <- dg[[pp[1]]] - dg[[pp[2]]]
      tibble(comparison = paste(pp, collapse = " vs "),
             diff_med   = median(dd),
             lo         = quantile(dd, .025),
             hi         = quantile(dd, .975),
             pd         = max(mean(dd > 0), mean(dd < 0)))
    }) %>% arrange(desc(pd))
}

tab_cells       <- cell_table(fit_female, df_female)
tab_species     <- contrast_table(fit_female, df_female)
tab_interaction <- interaction_table(fit_female, df_female)

tab_prior <- left_join(
  contrast_table(fit_female,      df_female),
  contrast_table(fit_female_wide, df_female),
  by = "pos_trt", suffix = c("_n2", "_n5")
)

cat("\n--- Posterior Choice Proportions per Cell (fit_female) ---\n")
print(tab_cells, n = Inf)

cat("\n--- Species Contrasts per Germplasm (janthinus - janthiniformis) ---\n")
print(tab_species, n = Inf)

cat("\n--- Interaction Contrasts: Difference of Species-Differences Across Germplasm ---\n")
print(tab_interaction, n = Inf)

cat("\n--- Prior Sensitivity: normal(0,2) vs normal(0,5) on Species Contrasts ---\n")
print(tab_prior, n = Inf)

# -- Freeze result tables and provenance --
dir.create("outputs", showWarnings = FALSE)

readr::write_csv(tab_cells,       "outputs/ytube_cell_proportions.csv")
readr::write_csv(tab_species,     "outputs/ytube_species_contrasts.csv")
readr::write_csv(tab_interaction, "outputs/ytube_interaction_contrasts.csv")
readr::write_csv(tab_prior,       "outputs/ytube_prior_sensitivity.csv")

readr::write_csv(
  df_female %>% count(species, pos_trt, wt = total),
  "outputs/ytube_design_counts.csv"
)
write.csv(
  unclass(xtabs(~ date + species, df_female)),
  "outputs/ytube_design_date_by_species.csv"
)

git_sha <- tryCatch(
  system("git rev-parse HEAD", intern = TRUE),
  error = function(e) "unknown", warning = function(w) "unknown"
)

writeLines(
  c(paste("generated", Sys.time()),
    paste("git sha", git_sha[1]),
    paste("brms", as.character(packageVersion("brms"))),
    paste("cmdstan", cmdstanr::cmdstan_version()),
    paste("n rows female", nrow(df_female)),
    paste("seed", 1),
    "--- fit_female (prior: normal(0,2)) ---",
    capture.output(print(fit_female)),
    "--- fit_female_wide (prior: normal(0,5)) ---",
    capture.output(print(fit_female_wide))),
  "outputs/ytube_session_info.txt"
)


nr_summary <- df %>%
  filter(sex == "female") %>%
  mutate(nr_rate = nr_count / (total + nr_count)) %>%
  select(species, pos_trt, pos_count, neg_count, nr_count, total, nr_rate)
print(nr_summary, n = Inf)

nr_summary %>% group_by(species) %>%
  summarise(nr_rate = sum(nr_count) / sum(nr_count + total))

df_female %>% group_by(species) %>%
  summarise(first = min(as.Date(as.character(date))),
            last  = max(as.Date(as.character(date))),
            n_days = n_distinct(date))

tab_cells <- tab_cells %>%
  left_join(df_female %>%
              group_by(species, pos_trt) %>%
              summarise(pos = sum(pos_count), neg = sum(neg_count),
                        n = sum(total), .groups = "drop"),
            by = c("species", "pos_trt"))
tab_cells

rh <- brms::rhat(fit_female)
es <- brms::neff_ratio(fit_female)
cat("max Rhat", max(rh, na.rm = TRUE), " min ESS ratio", min(es, na.rm = TRUE), "\n")

# =============================================================
# SECTION 3: VOC Whole-Profile Analysis (Bayesian Multivariate)
# =============================================================
cat("\n", strrep("=", 60), "\n")
cat("SECTION 3: VOC Whole-Profile Results\n")
cat(strrep("=", 60), "\n\n")

# -- Load or re-fit model --
if (file.exists("models/fit_voc.rds")) {
  fit <- readRDS("models/fit_voc.rds")
  cat("Loaded fit_voc from models/fit_voc.rds\n\n")
} else {
  cat("fit_voc.rds not found -- re-fitting model (this will take a while)...\n")
  source("scripts/bayes/bayesian_vocs.R")
  dir.create("models", showWarnings = FALSE)
  saveRDS(fit, "models/fit_voc.rds")
}


# -- Reconstruct model_df and epred --
source("scripts/01_volatile_proc.r")
stopifnot(exists("clr_df"), nrow(clr_df) == 107)

drop_sparse <- TRUE
sparse_cut  <- 0.70

meta_sub <- clr_df %>%
  select(Plant, Day) %>%
  mutate(Plant = factor(Plant), Day = factor(Day))

zero_frac <- colMeans(as.matrix(compounds) == 0)
keep_cmpd <- if (drop_sparse) names(zero_frac)[zero_frac <= sparse_cut] else names(zero_frac)

cmpd_key <- tibble(
  safe     = paste0("cmpd", sprintf("%02d", seq_along(keep_cmpd))),
  original = keep_cmpd
)

clr_sub <- clr_df %>%
  select(all_of(keep_cmpd)) %>%
  set_names(cmpd_key$safe)

cmpd_names <- cmpd_key$safe  
model_df   <- bind_cols(clr_sub, meta_sub)

plants <- levels(model_df$Plant)
nd     <- tibble(Plant = plants, Day = levels(model_df$Day)[1])

epred <- posterior_epred(fit, newdata = nd, re_formula = NA)
sigma <- posterior_epred(fit, newdata = nd, re_formula = NA, dpar = "sigma")

# -- Posterior centroids on PC1/PC2 --
cat("--- Posterior Centroids (PC1/PC2) per Germplasm ---\n")
pca_out   <- prcomp(clr_sub, scale. = FALSE)
med_cents <- apply(epred, c(2, 3), median)
cents_pc <- scale(med_cents, center = pca_out$center, scale = FALSE) %*% pca_out$rotation

centroid_summary <- as_tibble(cents_pc[, 1:2]) %>%
  set_names(c("PC1", "PC2")) %>%
  mutate(Plant = plants)
print(centroid_summary)

# -- Pairwise Aitchison distances (all 15 pairs) --
cat("\n--- Pairwise Aitchison Centroid Distances (all 15 pairs) ---\n")
pairs <- combn(plants, 2, simplify = FALSE)

pairwise_dist <- map_dfr(pairs, function(p) {
  i <- match(p[1], plants); j <- match(p[2], plants)
  d <- sqrt(rowSums((epred[, i, ] - epred[, j, ])^2))
  tibble(group1   = p[1], group2 = p[2],
         dist_med = median(d),
         dist_lo  = quantile(d, .025),
         dist_hi  = quantile(d, .975))
}) %>% arrange(dist_med)
print(pairwise_dist, n = Inf)

# -- Convergence summaries (HBR~RAD, D6Y6~Y6D6) --
cat("\n--- Convergence Summaries: Hybrid and Ecotype Pairs ---\n")
i0 <- match("D6", plants); j0 <- match("Y6", plants)
d_parents <- sqrt(rowSums((epred[, i0, ] - epred[, j0, ])^2))

conv_pairs  <- list(c("HBR", "RAD"), c("D6Y6", "Y6D6"))
convergence <- map_dfr(conv_pairs, function(p) {
  i <- match(p[1], plants); j <- match(p[2], plants)
  d_pair <- sqrt(rowSums((epred[, i, ] - epred[, j, ])^2))
  ratio  <- d_pair / d_parents
  tibble(pair                = paste(p, collapse = " vs "),
         dist_med            = median(d_pair),
         dist_lo             = quantile(d_pair, .025),
         dist_hi             = quantile(d_pair, .975),
         ratio_med           = median(ratio),
         ratio_lo            = quantile(ratio, .025),
         ratio_hi            = quantile(ratio, .975),
         p_less_than_parents = mean(d_pair < d_parents))
})
print(convergence, width = Inf)

# -- Per-germplasm dispersion --
cat("\n--- Per-Germplasm Dispersion (Bayesian betadisper analogue) ---\n")
disp <- sqrt(apply(sigma^2, c(1, 2), sum))
colnames(disp) <- plants

dispersion <- as_tibble(disp) %>%
  pivot_longer(everything(), names_to = "Plant", values_to = "d") %>%
  group_by(Plant) %>%
  summarise(med = median(d),
            lo  = quantile(d, .025),
            hi  = quantile(d, .975),
            .groups = "drop") %>%
  arrange(desc(med))
print(dispersion)

cat("\n--- Dispersion Pairwise Contrasts (pd > 0.95 only) ---\n")
disp_pairwise <- map_dfr(pairs, function(p) {
  d <- disp[, p[1]] - disp[, p[2]]
  tibble(group1   = p[1],
         group2   = p[2],
         diff_med = median(d),
         lo       = quantile(d, .025),
         hi       = quantile(d, .975),
         pd       = max(mean(d > 0), mean(d < 0)))
})
print(filter(disp_pairwise, pd > 0.95), n = Inf)

cat("\n", strrep("=", 60), "\n")
cat("All results printed.\n")
cat(strrep("=", 60), "\n")
