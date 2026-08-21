# =============================================================================
# Toadflax VOC whole-profile analysis -- BAYESIAN (estimation, not testing)
# Multivariate Gaussian model on CLR-transformed compounds directly.
#
# Produces the quantities the manuscript needs:
#   - posterior germplasm centroids in Aitchison space (with uncertainty)
#   - pairwise centroid-distance posteriors w/ 95% CrI  (all 15 pairs)
#   - convergence summaries for HBR~RAD, D6Y6~Y6D6
#   - per-germplasm dispersion posteriors (the Y6-variability result)
# =============================================================================

source("scripts/01_volatile_proc.r")   # clr_df, clr_mat, compounds
stopifnot(exists("clr_df"), nrow(clr_df) == 107)

library(tidyverse)
library(brms)
library(cmdstanr)

cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
options(mc.cores = parallel::detectCores())

# -----------------------------------------------------------------------------
# 1. Data prep -- drop ultra-sparse compounds (>70% zeros) as before,
#    then CLR-transform. No PCA.
# -----------------------------------------------------------------------------
drop_sparse <- TRUE
sparse_cut  <- 0.70

meta_sub <- clr_df %>%
  select(Plant, Day) %>%
  mutate(Plant = factor(Plant), Day = factor(Day))

zero_frac <- colMeans(as.matrix(compounds) == 0)
keep_cmpd <- if (drop_sparse) names(zero_frac)[zero_frac <= sparse_cut] else names(zero_frac)
cat("Modeling", length(keep_cmpd), "of", ncol(compounds), "compounds\n")

# Safe column names for brms -- short alphanumeric only
cmpd_key <- tibble(
  safe     = paste0("cmpd", sprintf("%02d", seq_along(keep_cmpd))),
  original = keep_cmpd
)

clr_sub <- clr_df %>%
  select(all_of(keep_cmpd)) %>%
  set_names(cmpd_key$safe)

cmpd_names <- cmpd_key$safe

model_df <- bind_cols(clr_sub, meta_sub)

# -----------------------------------------------------------------------------
# 2. Multivariate distributional model on CLR compounds directly
#    location: cmpd_k ~ Plant + (1|Day)
#    scale:    sigma   ~ Plant
#    set_rescor(FALSE): CLR compounds ARE correlated, but modeling the full
#    residual correlation matrix across 20 responses would require an
#    LKJ-prior correlation model that is extremely slow at this dimension.
#    Dropping rescor is a standard pragmatic choice for high-dimensional
#    multivariate brms; the centroid distances and dispersions we care about
#    are computed from the *mean* surface (epred), which is not affected by
#    residual correlation structure.
# -----------------------------------------------------------------------------
forms <- lapply(cmpd_names, function(v)
  bf(as.formula(paste0("`", v, "` ~ Plant + (1|Day)")), sigma ~ Plant))
mv_formula <- Reduce(`+`, forms) + set_rescor(FALSE)

priors <- do.call(c, lapply(cmpd_names, function(v) c(
  prior_string("normal(0, 5)",   class = "b",         resp = v),
  prior_string("normal(0, 2)",   class = "b",         resp = v, dpar = "sigma"),
  prior_string("normal(0, 2)",   class = "Intercept", resp = v, dpar = "sigma"),
  prior_string("exponential(1)", class = "sd",        resp = v, lb = 0)
)))

fit <- brm(mv_formula, data = model_df, family = gaussian(),
           prior = priors, chains = 4, iter = 4000, warmup = 1000,
           control = list(adapt_delta = 0.95, max_treedepth = 12),
           seed = 42, backend = "cmdstanr")

print(fit)   # check Rhat ~ 1.00, ESS healthy

# -----------------------------------------------------------------------------
# 3. Posterior centroids per germplasm in CLR (Aitchison) space
#    re_formula = NA marginalizes over the day random effect
# -----------------------------------------------------------------------------
plants <- levels(model_df$Plant)
nd     <- tibble(Plant = plants, Day = levels(model_df$Day)[1])

epred <- posterior_epred(fit, newdata = nd, re_formula = NA)  # [draws, 6, 20]
sigma <- posterior_epred(fit, newdata = nd, re_formula = NA, dpar = "sigma")

# Centroid summary on PC1/PC2 for the ordination figure -- run PCA on the
# posterior median centroids so the figure axes match the sample-level PCA
pca_out     <- prcomp(clr_sub, scale. = FALSE)   # same PCA as before, for plotting
med_cents   <- apply(epred, c(2, 3), median)     # [6, 20] median centroid matrix
cents_pc    <- med_cents %*% pca_out$rotation    # project into PC space

centroid_summary <- as_tibble(cents_pc[, 1:2]) %>%
  set_names(c("PC1", "PC2")) %>%
  mutate(Plant = plants)
print(centroid_summary)

# -----------------------------------------------------------------------------
# 4. Pairwise Aitchison centroid distances -- all 15 pairs
# -----------------------------------------------------------------------------
pairs <- combn(plants, 2, simplify = FALSE)

pairwise_dist <- map_dfr(pairs, function(p) {
  i <- match(p[1], plants); j <- match(p[2], plants)
  d <- sqrt(rowSums((epred[, i, ] - epred[, j, ])^2))  # per-draw Aitchison distance
  tibble(group1 = p[1], group2 = p[2],
         dist_med = median(d),
         dist_lo  = quantile(d, .025),
         dist_hi  = quantile(d, .975))
}) %>% arrange(dist_med)

print(pairwise_dist, n = Inf)

# -----------------------------------------------------------------------------
# 5. Convergence summaries for the two key pairs
#    Reports ratio to parental separation + P(d_pair < d_parents)
# -----------------------------------------------------------------------------
i0 <- match("D6", plants); j0 <- match("Y6", plants)
d_parents <- sqrt(rowSums((epred[, i0, ] - epred[, j0, ])^2))

conv_pairs <- list(c("HBR", "RAD"), c("D6Y6", "Y6D6"))
convergence <- map_dfr(conv_pairs, function(p) {
  i <- match(p[1], plants); j <- match(p[2], plants)
  d_pair <- sqrt(rowSums((epred[, i, ] - epred[, j, ])^2))
  ratio  <- d_pair / d_parents
  tibble(pair              = paste(p, collapse = " vs "),
         dist_med          = median(d_pair),
         dist_lo           = quantile(d_pair, .025),
         dist_hi           = quantile(d_pair, .975),
         ratio_med         = median(ratio),
         ratio_lo          = quantile(ratio, .025),
         ratio_hi          = quantile(ratio, .975),
         p_less_than_parents = mean(d_pair < d_parents))
})
print(convergence)

# -----------------------------------------------------------------------------
# 6. Per-germplasm dispersion (Bayesian betadisper analogue)
#    Scalar within-group spread = sqrt(sum of squared sigma across compounds)
# -----------------------------------------------------------------------------
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


print(convergence, width = Inf)
