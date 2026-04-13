source("volatile_proc.r")

stopifnot(
  exists("clr_df"),
  exists("clr_mat"),
  nrow(clr_df) == 107
)


library(vegan)

# Filter to whichever Types you want -- adjust as needed
clr_sub <- clr_df %>% filter(Type %in% c("DT", "YT"))
meta_sub <- clr_sub %>% select(`GC File`, ID, Plant, Type, Day)
mat_sub  <- clr_sub %>% select(all_of(colnames(compounds))) %>% as.matrix()

# PCA on CLR (Aitchison PCA)
pca_out <- prcomp(mat_sub, scale. = FALSE)  # CLR already variance-stabilized
summary(pca_out)

pca_scores <- as_tibble(pca_out$x[, 1:4]) %>%
  bind_cols(meta_sub, .)

# Quick plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Plant, shape = Type)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Plant), type = "norm", level = 0.68) +
  theme_bw() +
  labs(title = "Aitchison PCA of toadflax VOC profiles")

# PERMANOVA
set.seed(42)
adonis2(mat_sub ~ Plant + Day,
        data = meta_sub,
        method = "euclidean",   # Euclidean on CLR = Aitchison distance
        permutations = 999,
        by = "margin")

# Pairwise if overall is significant -- needs pairwiseAdonis or a loop
# library(pairwiseAdonis)
# pairwise.adonis2(mat_sub ~ Plant, data = meta_sub, method = "euclidean")