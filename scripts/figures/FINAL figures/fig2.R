# =============================================================
# figure_parental_axis.R
# Position of each germplasm relative to the Y6 -> D6 parental
# axis in full CLR space, showing both the on-axis coordinate
# and the orthogonal displacement off that axis.
#
# Assumes results_summary.R Section 3 has run, so these exist
# epred, clr_sub, model_df, plants, pairwise_dist
# =============================================================

library(tidyverse)
library(patchwork)

stopifnot(exists("epred"), exists("clr_sub"),
          exists("model_df"), exists("plants"), exists("pairwise_dist"))

X <- as.matrix(clr_sub)
stopifnot(ncol(X) == dim(epred)[3])   # sample columns must match response order

iD6 <- match("D6", plants)
iY6 <- match("Y6", plants)
stopifnot(!is.na(iD6), !is.na(iY6))

# scale_by_parents = TRUE puts both axes in units of the parental
# separation, so the reader compares off-axis displacement directly
# against the D6 to Y6 distance. FALSE leaves the y axis in raw
# Aitchison units.
scale_by_parents <- TRUE

pal <- c(Y6 = "#E8B84B", D6 = "#4C7FA8", Y6D6 = "#9BB47A", D6Y6 = "#7FA05C",
         HBR = "#B85C38", RAD = "#8C4A6B")

# -------------------------------------------------------------
# 1. Decompose each posterior centroid into an on-axis coordinate
#    and an orthogonal residual distance, one value per draw.
#    Axis runs Y6 -> D6, so 0 is the Y6 centroid and 1 is the D6
#    centroid. Values past 1 sit beyond the Dalmatian parent in
#    the direction away from yellow.
# -------------------------------------------------------------
anchor     <- epred[, iY6, ]                 # draws x compounds
axis_draws <- epred[, iD6, ] - anchor        # draws x compounds
axis_len   <- sqrt(rowSums(axis_draws^2))    # parental separation per draw

decompose <- function(g) {
  d      <- epred[, g, ] - anchor
  t      <- rowSums(d * axis_draws) / rowSums(axis_draws^2)
  resid  <- d - t * axis_draws
  r      <- sqrt(rowSums(resid^2))
  if (scale_by_parents) r <- r / axis_len
  tibble(t = t, r = r)
}

cent_draws <- map_dfr(seq_along(plants),
                      ~ decompose(.x) %>% mutate(Plant = plants[.x]))

cent_summary <- cent_draws %>%
  group_by(Plant) %>%
  summarise(t_med = median(t),
            t_lo  = quantile(t, .025),
            t_hi  = quantile(t, .975),
            r_med = median(r),
            r_lo  = quantile(r, .025),
            r_hi  = quantile(r, .975),
            p_beyond_D6 = mean(t > 1),
            .groups = "drop") %>%
  mutate(anchor_pt = Plant %in% c("D6", "Y6"))

cat("\n--- Posterior position relative to the Y6 -> D6 parental axis ---\n")
cat("t = on-axis coordinate, 0 = Y6 centroid, 1 = D6 centroid\n")
cat("r = orthogonal distance off the axis",
    if (scale_by_parents) "in units of the parental separation\n" else
      "in Aitchison units\n")
cat("p_beyond_D6 = Pr(t > 1). D6 and Y6 are anchors and are fixed by",
    "construction\n\n")
print(cent_summary %>% arrange(t_med), n = Inf)

# fraction of total displacement that is orthogonal, for the text
frac_off <- cent_draws %>%
  left_join(tibble(Plant = plants), by = "Plant") %>%
  group_by(Plant) %>%
  summarise(frac_med = median(r / sqrt(r^2 + (t * if (scale_by_parents) 1 else 1)^2 *
                                         1)), .groups = "drop")
cat("\n--- Note ---\n")
cat("If you want the orthogonal fraction of displacement, compute it inside",
    "decompose() where both components share units, rather than from the",
    "scaled summary above.\n")

# -------------------------------------------------------------
# 2. Sample-level coordinates, using the posterior median axis.
#    Fixing the axis keeps sample scatter interpretable as
#    among-plant variation rather than axis uncertainty.
# -------------------------------------------------------------
v_med   <- apply(axis_draws, 2, median)
a_med   <- apply(anchor,     2, median)
len_med <- sqrt(sum(v_med^2))

Xc <- sweep(X, 2, a_med, "-")
t_samp <- as.numeric(Xc %*% v_med) / sum(v_med^2)
resid_samp <- Xc - outer(t_samp, v_med)
r_samp <- sqrt(rowSums(resid_samp^2))
if (scale_by_parents) r_samp <- r_samp / len_med

sample_scores <- model_df %>%
  select(Plant, Day) %>%
  mutate(t = t_samp, r = r_samp)

lvl <- cent_summary %>% arrange(t_med) %>% pull(Plant)
sample_scores$Plant <- factor(sample_scores$Plant, levels = lvl)
cent_summary$Plant  <- factor(cent_summary$Plant,  levels = lvl)

# per-germplasm label offsets. dx and dy are in data units, so
# tune dy to the actual range of r once you see the plot.
lab_pos <- tribble(
  ~Plant,  ~dx,    ~dy,    ~hjust, ~vjust,
  "D6Y6",  -0.035, -0.012,  1,      1,      # lower left
  "HBR",   -0.035, -0.012,  1,      1,      # lower left
  "Y6D6",   0.035,  0.012,  0,      0,      # upper right
  "RAD",    0.035,  0.012,  0,      0,      # upper right
  "D6",     0.035,  0.012,  0,      0,
  "Y6",    -0.035, -0.012,  1,      1
)

cent_labs <- cent_summary %>%
  left_join(lab_pos, by = "Plant") %>%
  mutate(x_lab = t_med + dx, y_lab_pos = r_med + dy)

pA <- ggplot() +
  geom_vline(xintercept = c(0, 1), linetype = "dashed",
             colour = "grey55", linewidth = .4) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = .4) +
  geom_point(data = sample_scores,
             aes(t, r, colour = Plant),
             alpha = .28, size = 1.5, show.legend = FALSE) +
  geom_linerange(data = filter(cent_summary, !anchor_pt),
                 aes(x = t_med, ymin = r_lo, ymax = r_hi),
                 linewidth = .8, colour = "grey15") +
  geom_linerange(data = filter(cent_summary, !anchor_pt),
                 aes(y = r_med, xmin = t_lo, xmax = t_hi),
                 linewidth = .8, colour = "grey15") +
  geom_point(data = filter(cent_summary, !anchor_pt),
             aes(t_med, r_med, fill = Plant),
             size = 3.6, shape = 21, colour = "grey15", stroke = .9) +
  geom_point(data = filter(cent_summary, anchor_pt),
             aes(t_med, r_med, fill = Plant),
             size = 3.6, shape = 23, colour = "grey15", stroke = .9) +
  geom_text(data = cent_labs,
            aes(x_lab, y_lab_pos, label = Plant,
                hjust = hjust, vjust = vjust),
            size = 3.2, colour = "grey15") +
  annotate("text", x = 0, y = Inf, label = "Y6 centroid", vjust = 1.6,
           size = 3, colour = "grey35") +
  annotate("text", x = 1, y = Inf, label = "D6 centroid", vjust = 1.6,
           size = 3, colour = "grey35") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = "Position on the Y6 \u2192 D6 parental axis",
       y = "Off-axis displacement\n(parental separation)", tag = "A") +
  theme_classic(base_size = 11) +
  theme(plot.margin = margin(18, 12, 6, 6))

# -------------------------------------------------------------
# 3. Forest plot of the 15 pairwise Aitchison centroid distances,
#    with the parental separation as a reference band and the
#    wild-hybrid-to-yellow distances broken out.
# -------------------------------------------------------------
d_parents_summ <- {
  d <- sqrt(rowSums((epred[, iD6, ] - epred[, iY6, ])^2))
  tibble(med = median(d), lo = quantile(d, .025), hi = quantile(d, .975))
}

is_pair <- function(g1, g2, a, b) (g1 == a & g2 == b) | (g1 == b & g2 == a)

pw <- pairwise_dist %>%
  mutate(pair = paste(group1, group2, sep = " \u2013 "),
         type = case_when(
           is_pair(group1, group2, "HBR", "RAD")   ~ "Convergence pair",
           is_pair(group1, group2, "D6Y6", "Y6D6") ~ "Convergence pair",
           is_pair(group1, group2, "D6", "Y6")     ~ "Parental separation",
           is_pair(group1, group2, "HBR", "Y6")    ~ "Wild hybrid vs yellow",
           is_pair(group1, group2, "RAD", "Y6")    ~ "Wild hybrid vs yellow",
           TRUE                                    ~ "Other"),
         pair = fct_reorder(pair, dist_med))

pB <- ggplot(pw, aes(dist_med, pair, colour = type)) +
  annotate("rect", xmin = d_parents_summ$lo, xmax = d_parents_summ$hi,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = .5) +
  geom_vline(xintercept = d_parents_summ$med, linetype = "dashed",
             colour = "grey35", linewidth = .4) +
  geom_linerange(aes(xmin = dist_lo, xmax = dist_hi), linewidth = .7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("Convergence pair"      = "#B85C38",
                                 "Parental separation"   = "#4C7FA8",
                                 "Wild hybrid vs yellow" = "#8C4A6B",
                                 "Other"                 = "grey55"),
                      name = NULL) +
  labs(x = "Posterior Aitchison centroid distance (95% CrI)",
       y = NULL, tag = "B") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8.5))

fig <- pA / pB + plot_layout(heights = c(1.15, 1))

ggsave("figure2_parental_axis.png", fig,
       width = 7.2, height = 8.4, dpi = 400, bg = "white")
ggsave("figure2_parental_axis.pdf", fig,
       width = 7.2, height = 8.4, device = cairo_pdf)

fig