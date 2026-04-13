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
# Use all samples (clr_sub <- clr_df); swap in a filter() call here if you
# want to restrict to a subset of Type or Plant later
# -----------------------------------------------------------------------------

clr_sub  <- clr_df
meta_sub <- clr_sub %>% select(`GC File`, ID, Plant, Type, Day)
mat_sub  <- clr_sub %>% select(all_of(colnames(compounds))) %>% as.matrix()

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

# -----------------------------------------------------------------------------
# PERMANOVA (adonis2)
# Euclidean distance on CLR = Aitchison distance
# by = "margin" tests each term after accounting for the other (type III)
# -----------------------------------------------------------------------------

set.seed(42)
adonis2(mat_sub ~ Plant + Day,
        data         = meta_sub,
        method       = "euclidean",
        permutations = 999,
        by           = "margin")

# -----------------------------------------------------------------------------
# Pairwise PERMANOVA
# vegan does not provide built-in pairwise tests, so we loop over all pairs,
# subset the data, run adonis2 on each, and apply Bonferroni correction
# across all comparisons (n = 15 pairs for 6 plant types)
# -----------------------------------------------------------------------------

plants <- unique(meta_sub$Plant)
pairs  <- combn(plants, 2, simplify = FALSE)

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

# -----------------------------------------------------------------------------
# Multivariate dispersion (betadisper)
# Tests whether within-group spread differs across plant types
# A significant result means PERMANOVA differences may reflect dispersion
# rather than (or in addition to) centroid shifts -- check Tukey output
# to identify which specific groups are driving heterogeneity
# -----------------------------------------------------------------------------

dist_sub <- dist(mat_sub, method = "euclidean")
bd       <- betadisper(dist_sub, meta_sub$Plant)

anova(bd)                         # parametric F-test
permutest(bd, permutations = 999) # permutation-based test (prefer this one)
plot(bd)                          # PCoA of distances to group centroids
TukeyHSD(bd)                      # pairwise dispersion comparisons