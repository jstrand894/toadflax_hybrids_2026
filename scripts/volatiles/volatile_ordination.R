# =============================================================================
# Toadflax VOC multivariate analysis
# Aitchison PCA, PERMANOVA, pairwise PERMANOVA, and dispersion testing
# =============================================================================

source("scripts/volatiles/volatile_proc.r")

stopifnot(
  exists("clr_df"),
  exists("clr_mat"),
  nrow(clr_df) == 107
)

library(vegan)

# -----------------------------------------------------------------------------
# Data preparation
# Plant is coerced to factor here to ensure all six levels are retained
# and contrasts are handled correctly downstream.
# Use all samples (clr_sub <- clr_df); swap in a filter() call here if you
# want to restrict to a subset of Type or Plant later
# -----------------------------------------------------------------------------

clr_sub  <- clr_df
meta_sub <- clr_sub %>%
  select(`GC File`, ID, Plant, Type, Day) %>%
  mutate(Plant = factor(Plant))

mat_sub <- clr_sub %>%
  select(all_of(colnames(compounds))) %>%
  as.matrix()

clr_sub %>% count(`Type`)
R.version
# -----------------------------------------------------------------------------
# Aitchison PCA
# CLR-transformed data are already variance-stabilized, so scale. = FALSE
# preserves the compositional geometry (no additional standardization needed)
# -----------------------------------------------------------------------------

pca_out    <- prcomp(mat_sub, scale. = FALSE)
summary(pca_out)

pca_scores <- as_tibble(pca_out$x[, 1:4]) %>%
  bind_cols(meta_sub, .) 

ggplot(pca_scores, aes(x = PC1, y = PC2, color = Plant, shape = Type)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Plant), type = "norm", level = 0.68) +
  theme_bw() +
  labs(title = "Aitchison PCA of toadflax VOC profiles")



library(patchwork)

# Variance explained for axis labels
pct_var <- summary(pca_out)$importance["Proportion of Variance", ] * 100
pc1_lab <- paste0("PC1 (", round(pct_var[1], 1), "% of variance)")
pc2_lab <- paste0("PC2 (", round(pct_var[2], 1), "% of variance)")

# Color and shape mappings (full set, subsetting handled per plot)
color_vals <- c(
  "D6"   = "#1D9E75",
  "Y6"   = "#D85A30",
  "D6Y6" = "#7F77DD",
  "Y6D6" = "#BA7517",
  "HBR"  = "#378ADD",
  "RAD"  = "#D4537E"
)

label_vals <- c(
  "D6"   = "D6 (Dalmatian)",
  "Y6"   = "Y6 (Yellow)",
  "D6Y6" = "D6Y6 (Dal. \u2640)",
  "Y6D6" = "Y6D6 (Yel. \u2640)",
  "HBR"  = "HBR (Boulder R.)",
  "RAD"  = "RAD (Radersburg)"
)

# Shared theme
theme_ord <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.major     = element_blank(),
      panel.grid.minor     = element_blank(),
      legend.position      = "bottom",
      legend.title         = element_blank(),
      legend.text          = element_text(size = 9),
      legend.key.size      = unit(0.5, "cm"),
      axis.text            = element_text(size = 9),
      axis.title           = element_text(size = 10),
      plot.title           = element_text(size = 10, face = "bold", hjust = 0.5),
      strip.background     = element_blank(),
      strip.text           = element_text(size = 10, face = "bold")
    )
}

# Shared geom layers -- call with a 2-element character vector of Plant names
ord_layers <- function(plants) {
  list(
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey70"),
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey70"),
    geom_point(alpha = 0.22, size = 2, stroke = 0),
    stat_ellipse(aes(group = Plant), type = "norm", level = 0.95,
                 geom = "polygon", color = NA, alpha = 0.10),
    stat_ellipse(aes(group = Plant), type = "norm", level = 0.95,
                 linewidth = 0.6),
    scale_color_manual(values = color_vals[plants], labels = label_vals[plants]),
    scale_fill_manual(values  = color_vals[plants], labels = label_vals[plants]),
    labs(x = pc1_lab, y = pc2_lab, color = NULL, fill = NULL),
    guides(
      color = guide_legend(override.aes = list(size = 3, alpha = 1)),
      fill  = "none"
    )
  )
}

# Centroids as explicit data (most reliable approach)
centroids <- pca_scores %>%
  group_by(Plant) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

# --- Plot 1: Parental lines ---
scores_p1 <- pca_scores %>% filter(Plant %in% c("D6", "Y6"))
cents_p1  <- centroids   %>% filter(Plant %in% c("D6", "Y6"))

p1 <- ggplot(scores_p1, aes(x = PC1, y = PC2, color = Plant, fill = Plant)) +
  ord_layers(c("D6", "Y6")) +
  geom_point(data = cents_p1, size = 5,   stroke = 1.8, color = "white") +
  geom_point(data = cents_p1, size = 3.5, stroke = 0) +
  ggtitle("Parental lines") +
  theme_ord()

# --- Plot 2: F1 hybrids ---
scores_p2 <- pca_scores %>% filter(Plant %in% c("D6Y6", "Y6D6"))
cents_p2  <- centroids   %>% filter(Plant %in% c("D6Y6", "Y6D6"))

p2 <- ggplot(scores_p2, aes(x = PC1, y = PC2, color = Plant, fill = Plant)) +
  ord_layers(c("D6Y6", "Y6D6")) +
  geom_point(data = cents_p2, size = 5,   stroke = 1.8, color = "white") +
  geom_point(data = cents_p2, size = 3.5, stroke = 0) +
  ggtitle(expression(bold(F[1]~hybrids))) +
  theme_ord()

# --- Plot 3: Wild hybrids ---
scores_p3 <- pca_scores %>% filter(Plant %in% c("HBR", "RAD"))
cents_p3  <- centroids   %>% filter(Plant %in% c("HBR", "RAD"))

p3 <- ggplot(scores_p3, aes(x = PC1, y = PC2, color = Plant, fill = Plant)) +
  ord_layers(c("HBR", "RAD")) +
  geom_point(data = cents_p3, size = 5,   stroke = 1.8, color = "white") +
  geom_point(data = cents_p3, size = 3.5, stroke = 0) +
  ggtitle("Wild hybrids") +
  theme_ord()

# Get shared axis limits with a little padding
ellipse_range <- function(data, level = 0.95) {
  get_ellipse <- function(df) {
    mu     <- colMeans(df[, c("PC1", "PC2")])
    cv     <- cov(df[, c("PC1", "PC2")])
    t_val  <- sqrt(qchisq(level, df = 2))
    theta  <- seq(0, 2 * pi, length.out = 200)
    circle <- cbind(cos(theta), sin(theta))
    ell    <- circle %*% (t_val * chol(cv)) +
      matrix(mu, nrow = 200, ncol = 2, byrow = TRUE)
    tibble(x = ell[, 1], y = ell[, 2])   # named x/y, not V1/V2
  }
  pts <- data %>%
    group_by(Plant) %>%
    group_modify(~ get_ellipse(.x)) %>%
    ungroup()
  list(
    x = range(pts$x),
    y = range(pts$y)
  )
}

ell_lim <- ellipse_range(pca_scores)
pad_frac <- 0.05

pc1_lim <- ell_lim$x + c(-1, 1) * diff(ell_lim$x) * pad_frac
pc2_lim <- ell_lim$y + c(-1, 1) * diff(ell_lim$y) * pad_frac

fixed_axes <- list(
  scale_x_continuous(limits = pc1_lim, expand = expansion(0)),
  scale_y_continuous(limits = pc2_lim, expand = expansion(0)),
  coord_cartesian(xlim = pc1_lim, ylim = pc2_lim, clip = "on")
)

p2 <- p2 + ggtitle(expression(bold(F[1]~hybrids)))

p1 <- p1 + fixed_axes
p2 <- p2 + fixed_axes
p3 <- p3 + fixed_axes

p1 + p2 + p3 + plot_layout(nrow = 1)

# -----------------------------------------------------------------------------
# PERMANOVA -- full dataset, unblocked (exploratory)
# Euclidean distance on CLR = Aitchison distance
# by = "margin" tests each term after accounting for the other (type III)
#
# NOTE: HBR and RAD were sampled on different days than D6, D6Y6, Y6, Y6D6,
# so Plant and Day are partially confounded in this model. Treat these results
# as exploratory; the blocked models below are the primary analyses.
# -----------------------------------------------------------------------------

set.seed(42)
adonis2(mat_sub ~ Plant + Day,
        data         = meta_sub,
        method       = "euclidean",
        permutations = 999,
        by           = "margin")

# -----------------------------------------------------------------------------
# PERMANOVA -- full dataset, Day-blocked (primary)
# Restricting permutations within days removes day-to-day variation and
# avoids inflating Plant effects due to the confounded sampling design
# -----------------------------------------------------------------------------

set.seed(42)
adonis2(mat_sub ~ Plant,
        data         = meta_sub,
        method       = "euclidean",
        permutations = how(blocks = meta_sub$Day, nperm = 999))

# -----------------------------------------------------------------------------
# PERMANOVA -- co-sampled plant types only (D6, D6Y6, Y6, Y6D6)
# These four types were collected on the same days, so day blocking is fully
# legitimate here -- permutations shuffle within days containing all four types
# -----------------------------------------------------------------------------

meta_cosamp <- meta_sub %>% filter(Plant %in% c("D6", "D6Y6", "Y6", "Y6D6"))
mat_cosamp  <- mat_sub[meta_sub$Plant %in% c("D6", "D6Y6", "Y6", "Y6D6"), ]

set.seed(42)
adonis2(mat_cosamp ~ Plant,
        data         = meta_cosamp,
        method       = "euclidean",
        permutations = how(blocks = meta_cosamp$Day, nperm = 999))

# -----------------------------------------------------------------------------
# PERMANOVA -- HBR vs RAD only
# HBR and RAD share the same sampling days (2, 4, 8), so this is a clean
# blocked comparison within that sampling group
# -----------------------------------------------------------------------------

meta_hbr_rad <- meta_sub %>% filter(Plant %in% c("HBR", "RAD"))
mat_hbr_rad  <- mat_sub[meta_sub$Plant %in% c("HBR", "RAD"), ]

set.seed(42)
adonis2(mat_hbr_rad ~ Plant,
        data         = meta_hbr_rad,
        method       = "euclidean",
        permutations = how(blocks = meta_hbr_rad$Day, nperm = 999))

# -----------------------------------------------------------------------------
# Pairwise PERMANOVA -- unblocked
# vegan does not provide built-in pairwise tests, so we loop over all pairs,
# subset the data, run adonis2 on each, and apply Bonferroni correction
# across all comparisons (n = 15 pairs for 6 plant types).
# Day is not included here, which is standard for pairwise follow-ups, but
# note that cross-block comparisons involving HBR or RAD should be interpreted
# cautiously given the confounded sampling design
# -----------------------------------------------------------------------------

plants <- unique(meta_sub$Plant)
pairs  <- combn(plants, 2, simplify = FALSE)

set.seed(42)
pairwise_results <- map_dfr(pairs, function(p) {
  idx       <- meta_sub$Plant %in% p
  mat_pair  <- mat_sub[idx, ]
  meta_pair <- meta_sub[idx, ]
  
  res <- adonis2(mat_pair ~ Plant,
                 data         = meta_pair,
                 method       = "euclidean",
                 permutations = 999)
  
  tibble(
    group1 = p[1],
    group2 = p[2],
    F_val  = res$F[1],
    R2     = res$R2[1],
    p_val  = res$`Pr(>F)`[1]
  )
})

pairwise_results <- pairwise_results %>%
  mutate(p_adj = p.adjust(p_val, method = "bonferroni"))

print(pairwise_results)

# =============================================================================
# Figure 2 -- Pairwise PERMANOVA R² heatmap (6x6)
# Bonferroni-corrected significance marked with asterisk and bold tile border
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Run pairwise PERMANOVA (copied from volatile_ordination.R for self-containment)
# Skip this block if pairwise_results already exists in your environment
# -----------------------------------------------------------------------------

plants <- levels(meta_sub$Plant)   # preserves factor order
pairs  <- combn(plants, 2, simplify = FALSE)

set.seed(42)
pairwise_results <- purrr::map_dfr(pairs, function(p) {
  idx       <- meta_sub$Plant %in% p
  mat_pair  <- mat_sub[idx, ]
  meta_pair <- meta_sub[idx, ]
  
  res <- vegan::adonis2(mat_pair ~ Plant,
                        data         = meta_pair,
                        method       = "euclidean",
                        permutations = 999)
  
  tibble(
    group1 = p[1],
    group2 = p[2],
    F_val  = res$F[1],
    R2     = res$R2[1],
    p_val  = res$`Pr(>F)`[1]
  )
}) %>%
  mutate(p_adj = p.adjust(p_val, method = "bonferroni"),
         sig   = p_adj < 0.05)



# -----------------------------------------------------------------------------
# Multivariate dispersion (betadisper)
# Tests whether within-group spread differs across plant types.
# A significant result means PERMANOVA differences may reflect dispersion
# rather than (or in addition to) centroid shifts -- check Tukey output
# to identify which specific groups are driving heterogeneity.
# Y6 showed greater dispersion than D6Y6 and RAD in Tukey follow-up;
# interpret pairwise comparisons involving Y6 accordingly
# -----------------------------------------------------------------------------
pca_out$rotation
dist_sub <- dist(mat_sub, method = "euclidean")
bd       <- betadisper(dist_sub, meta_sub$Plant)

anova(bd)                         # parametric F-test
permutest(bd, permutations = 999) # permutation-based test (prefer this one)
plot(bd)                          # PCoA of distances to group centroids
TukeyHSD(bd)                      # pairwise dispersion comparisons





