# =============================================================================
# Toadflax VOC ordination figure (Figure 2)
# Aitchison PCA on the 20-compound modeled set, with model-based centroids
# Plots the same compound set and the same posterior-median centroids
# reported in the Results, so figure and text come from one analysis.
# =============================================================================

source("scripts/01_volatile_proc.r")

stopifnot(
  exists("clr_df"),
  exists("compounds"),
  nrow(clr_df) == 107
)

library(tidyverse)
library(patchwork)
library(brms)

# -----------------------------------------------------------------------------
# 1. Reproduce the modeled compound set (>70% zero compounds dropped),
#    using the identical safe-name scheme as bayesian_vocs.R so the PCA
#    columns line up with the brms response order.
# -----------------------------------------------------------------------------
drop_sparse <- TRUE
sparse_cut  <- 0.70

zero_frac <- colMeans(as.matrix(compounds) == 0)
keep_cmpd <- if (drop_sparse) names(zero_frac)[zero_frac <= sparse_cut] else names(zero_frac)

cmpd_key <- tibble(
  safe     = paste0("cmpd", sprintf("%02d", seq_along(keep_cmpd))),
  original = keep_cmpd
)

clr_sub <- clr_df %>%
  select(all_of(keep_cmpd)) %>%
  set_names(cmpd_key$safe)

meta_sub <- clr_df %>%
  select(Plant, Type) %>%
  mutate(Plant = factor(Plant))

# -----------------------------------------------------------------------------
# 2. Aitchison PCA on the 20-compound set (same matrix the model used)
#    prcomp centers by default; that center is needed below to project the
#    model centroids into the same coordinate frame as the sample scores.
# -----------------------------------------------------------------------------
pca_out <- prcomp(clr_sub, scale. = FALSE)

# Flip PC1 sign so positive PC1 = high monoterpene content, matching the
# Results loadings text (Y6 positive, wild hybrids negative). PCA sign is
# arbitrary; negate scores and rotation together to keep them consistent.
pca_out$x[, 1]        <- -pca_out$x[, 1]
pca_out$rotation[, 1] <- -pca_out$rotation[, 1]

pct_var <- summary(pca_out)$importance["Proportion of Variance", ] * 100
pc1_lab <- paste0("PC1 (", round(pct_var[1], 1), "% of variance)")
pc2_lab <- paste0("PC2 (", round(pct_var[2], 1), "% of variance)")

pca_scores <- as_tibble(pca_out$x[, 1:2]) %>%
  bind_cols(meta_sub, .)

# -----------------------------------------------------------------------------
# 3. Model-based centroids (posterior medians), projected into PC space.
#    These are the SAME centroids reported in the Results table, so the
#    plotted dot equals the reported number. Must subtract pca_out$center
#    before projecting, since prcomp scores are on centered data.
# -----------------------------------------------------------------------------
fit <- readRDS("models/fit_voc.rds")

plants <- levels(meta_sub$Plant)
nd     <- tibble(Plant = plants, Day = levels(factor(clr_df$Day))[1])

epred     <- posterior_epred(fit, newdata = nd, re_formula = NA)  # [draws, 6, 20]
med_cents <- apply(epred, c(2, 3), median)                        # [6, 20]

cents_pc <- sweep(med_cents, 2, pca_out$center) %*% pca_out$rotation
centroids <- as_tibble(cents_pc[, 1:2]) %>%
  set_names(c("PC1", "PC2")) %>%
  mutate(Plant = plants)

# -----------------------------------------------------------------------------
# 4. Aesthetics
# -----------------------------------------------------------------------------
color_vals <- c(
  "D6"   = "#1D9E75", "Y6"   = "#D85A30",
  "D6Y6" = "#7F77DD", "Y6D6" = "#BA7517",
  "HBR"  = "#378ADD", "RAD"  = "#D4537E"
)

# NOTE: confirm maternal-parent notation matches the Methods crossing scheme.
# This asserts D6Y6 is Dalmatian-maternal and Y6D6 is yellow-maternal.
label_vals <- c(
  "D6"   = "D6 (Dalmatian)",     "Y6"   = "Y6 (Yellow)",
  "D6Y6" = "D6Y6 (Dal. \u2640)", "Y6D6" = "Y6D6 (Yel. \u2640)",
  "HBR"  = "HBR (Boulder R.)",   "RAD"  = "RAD (Radersburg)"
)

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

# Ellipse level set to 0.68 to match the figure caption.
ell_level <- 0.68

ord_layers <- function(plants) {
  list(
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey70"),
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey70"),
    geom_point(alpha = 0.22, size = 2, stroke = 0),
    stat_ellipse(aes(group = Plant), type = "norm", level = ell_level,
                 geom = "polygon", color = NA, alpha = 0.10),
    stat_ellipse(aes(group = Plant), type = "norm", level = ell_level,
                 linewidth = 0.6),
    scale_color_manual(values = color_vals[plants], labels = label_vals[plants]),
    scale_fill_manual(values  = color_vals[plants], labels = label_vals[plants]),
    labs(x = pc1_lab, y = pc2_lab, color = NULL, fill = NULL),
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
           fill  = "none")
  )
}

build_panel <- function(grp, title) {
  scores <- pca_scores %>% filter(Plant %in% grp)
  cents  <- centroids  %>% filter(Plant %in% grp)
  ggplot(scores, aes(x = PC1, y = PC2, color = Plant, fill = Plant)) +
    ord_layers(grp) +
    geom_point(data = cents, size = 5,   stroke = 1.8, color = "white") +
    geom_point(data = cents, size = 3.5, stroke = 0) +
    ggtitle(title) +
    theme_ord()
}

p1 <- build_panel(c("D6", "Y6"),     "Parental lines")
p2 <- build_panel(c("D6Y6", "Y6D6"), expression(bold(F[1]~hybrids)))
p3 <- build_panel(c("HBR", "RAD"),   "Wild hybrids")

# -----------------------------------------------------------------------------
# 5. Shared axis limits across panels, computed from the 0.68 ellipses of all
#    groups so the wild-hybrid displacement is readable across panels.
#    Limits applied via coord_cartesian only (zoom, not exclusion), so no
#    extreme point is silently dropped.
# -----------------------------------------------------------------------------
ellipse_range <- function(data, level = ell_level) {
  get_ellipse <- function(df) {
    mu     <- colMeans(df[, c("PC1", "PC2")])
    cv     <- cov(df[, c("PC1", "PC2")])
    t_val  <- sqrt(qchisq(level, df = 2))
    theta  <- seq(0, 2 * pi, length.out = 200)
    circle <- cbind(cos(theta), sin(theta))
    ell    <- circle %*% (t_val * chol(cv)) +
      matrix(mu, nrow = 200, ncol = 2, byrow = TRUE)
    tibble(x = ell[, 1], y = ell[, 2])
  }
  pts <- data %>% group_by(Plant) %>% group_modify(~ get_ellipse(.x)) %>% ungroup()
  list(x = range(pts$x), y = range(pts$y))
}

ell_lim  <- ellipse_range(pca_scores)
pad_frac <- 0.05
pc1_lim  <- ell_lim$x + c(-1, 1) * diff(ell_lim$x) * pad_frac
pc2_lim  <- ell_lim$y + c(-1, 1) * diff(ell_lim$y) * pad_frac

fixed_axes <- coord_cartesian(xlim = pc1_lim, ylim = pc2_lim, clip = "on")

p1 <- p1 + fixed_axes
p2 <- p2 + fixed_axes
p3 <- p3 + fixed_axes

fig2 <- p1 + p2 + p3 + plot_layout(nrow = 1)
fig2

ggsave("figures/fig2_ordination.pdf", 
       fig2, width = 11, height = 4,
       device = cairo_pdf)
