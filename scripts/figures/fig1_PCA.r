source("scripts/01_volatile_proc.R")

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Aitchison PCA ──────────────────────────────────────────────────────────────
# CLR-transformed data are already variance-stabilized, so scale. = FALSE
# preserves compositional geometry (no additional standardization needed)

mat_sub <- clr_df %>%
  select(all_of(colnames(compounds))) %>%
  as.matrix()

meta_sub <- clr_df %>%
  select(`GC File`, ID, Plant, Type, Day) %>%
  mutate(Plant = factor(Plant))

pca_out    <- prcomp(mat_sub, scale. = FALSE)
pca_scores <- as_tibble(pca_out$x[, 1:4]) %>%
  bind_cols(meta_sub, .)

# ── Axis labels with variance explained ───────────────────────────────────────
pct_var <- summary(pca_out)$importance["Proportion of Variance", ] * 100
pc1_lab <- paste0("PC1 (", round(pct_var[1], 1), "% of variance)")
pc2_lab <- paste0("PC2 (", round(pct_var[2], 1), "% of variance)")

# ── Color and legend label mappings ───────────────────────────────────────────
color_vals <- c(
  "D6"   = "#0F6E56",
  "Y6"   = "#D85A30",
  "D6Y6" = "#5DCAA5",
  "Y6D6" = "#F0997B",
  "HBR"  = "#EF9F27",
  "RAD"  = "#BA7517"
)

label_vals <- c(
  "D6"   = "DT6 (Dalmatian)",
  "Y6"   = "YT6 (Yellow)",
  "D6Y6" = "D6Y6 (Dal. \u2640)",
  "Y6D6" = "Y6D6 (Yel. \u2640)",
  "HBR"  = "HBR (Boulder R.)",
  "RAD"  = "RAD (Radersburg)"
)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_ord <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 9),
      legend.key.size  = unit(0.5, "cm"),
      axis.text        = element_text(size = 9),
      axis.title       = element_text(size = 10),
      plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text       = element_text(size = 10, face = "bold")
    )
}

# ── Shared geom layers ────────────────────────────────────────────────────────
# `plants` is a character vector of the Plant levels shown in a given panel
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

# ── Centroids (mean PC1/PC2 per germplasm) ────────────────────────────────────
centroids <- pca_scores %>%
  group_by(Plant) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

# ── Shared axis limits across all three panels ────────────────────────────────
# Computed from the full 95% ellipse range so all panels are directly comparable

ellipse_range <- function(data, level = 0.95) {
  get_ellipse <- function(df) {
    mu    <- colMeans(df[, c("PC1", "PC2")])
    cv    <- cov(df[, c("PC1", "PC2")])
    t_val <- sqrt(qchisq(level, df = 2))
    theta <- seq(0, 2 * pi, length.out = 200)
    ell   <- cbind(cos(theta), sin(theta)) %*% (t_val * chol(cv)) +
      matrix(mu, nrow = 200, ncol = 2, byrow = TRUE)
    tibble(x = ell[, 1], y = ell[, 2])
  }
  pts <- data %>%
    group_by(Plant) %>%
    group_modify(~ get_ellipse(.x)) %>%
    ungroup()
  list(x = range(pts$x), y = range(pts$y))
}

ell_lim  <- ellipse_range(pca_scores)
pad_frac <- 0.05
pc1_lim  <- ell_lim$x + c(-1, 1) * diff(ell_lim$x) * pad_frac
pc2_lim  <- ell_lim$y + c(-1, 1) * diff(ell_lim$y) * pad_frac

fixed_axes <- list(
  scale_x_continuous(limits = pc1_lim, expand = expansion(0)),
  scale_y_continuous(limits = pc2_lim, expand = expansion(0)),
  coord_cartesian(xlim = pc1_lim, ylim = pc2_lim, clip = "on")
)

# ── Panel 1: Parental lines ───────────────────────────────────────────────────
scores_p1 <- pca_scores %>% filter(Plant %in% c("D6", "Y6"))
cents_p1  <- centroids   %>% filter(Plant %in% c("D6", "Y6"))

p1 <- ggplot(scores_p1, aes(x = PC1, y = PC2, color = Plant, fill = Plant)) +
  ord_layers(c("D6", "Y6")) +
  geom_point(data = cents_p1, size = 5,   stroke = 1.8, color = "white") +
  geom_point(data = cents_p1, size = 3.5, stroke = 0) +
  ggtitle("Parental lines") +
  fixed_axes +
  theme_ord()

# ── Panel 2: F1 hybrids ───────────────────────────────────────────────────────
scores_p2 <- pca_scores %>% filter(Plant %in% c("D6Y6", "Y6D6"))
cents_p2  <- centroids   %>% filter(Plant %in% c("D6Y6", "Y6D6"))

p2 <- ggplot(scores_p2, aes(x = PC1, y = PC2, color = Plant, fill = Plant)) +
  ord_layers(c("D6Y6", "Y6D6")) +
  geom_point(data = cents_p2, size = 5,   stroke = 1.8, color = "white") +
  geom_point(data = cents_p2, size = 3.5, stroke = 0) +
  ggtitle(expression(bold(F[1] ~ hybrids))) +
  fixed_axes +
  theme_ord()

# ── Panel 3: Wild hybrids ─────────────────────────────────────────────────────
scores_p3 <- pca_scores %>% filter(Plant %in% c("HBR", "RAD"))
cents_p3  <- centroids   %>% filter(Plant %in% c("HBR", "RAD"))

p3 <- ggplot(scores_p3, aes(x = PC1, y = PC2, color = Plant, fill = Plant)) +
  ord_layers(c("HBR", "RAD")) +
  geom_point(data = cents_p3, size = 5,   stroke = 1.8, color = "white") +
  geom_point(data = cents_p3, size = 3.5, stroke = 0) +
  ggtitle("Wild hybrids") +
  fixed_axes +
  theme_ord()

# ── Combined figure ───────────────────────────────────────────────────────────
fig1 <- p1 + p2 + p3 + plot_layout(nrow = 1)
fig1

ggsave("figures/fig1_PCA.pdf", fig1,
       width = 10, height = 4, device = cairo_pdf)
