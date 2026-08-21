# =============================================================
# figure_parental_axis.R
# Projection onto the D6 -> Y6 parental axis in full CLR space
# plus a forest plot of the 15 pairwise centroid distances
#
# Assumes results_summary.R Section 3 has already run, so these
# exist: fit, epred, clr_sub, model_df, plants, pairwise_dist
# =============================================================

library(tidyverse)
library(patchwork)

stopifnot(exists("epred"), exists("clr_sub"), exists("model_df"), exists("plants"))

X <- as.matrix(clr_sub)
stopifnot(ncol(X) == dim(epred)[3])   # sample columns must match response order

iD6 <- match("D6", plants)
iY6 <- match("Y6", plants)

# -------------------------------------------------------------
# 1. Build the parental axis, one per posterior draw
#    Axis runs Y6 -> D6. Scores are scaled so that
#    0 = Y6 centroid and 1 = D6 centroid. Anything > 1 sits
#    beyond the Dalmatian parent in the direction away from Y6.
# -------------------------------------------------------------
anchor     <- epred[, iY6, ]                 # draws x compounds
axis_draws <- epred[, iD6, ] - anchor        # draws x compounds

# posterior centroid score for every germplasm, per draw
cent_draws <- sapply(seq_along(plants), function(g) {
  rowSums((epred[, g, ] - anchor) * axis_draws) / rowSums(axis_draws^2)
})
colnames(cent_draws) <- plants

cent_summary <- as_tibble(cent_draws) %>%
  pivot_longer(everything(), names_to = "Plant", values_to = "score") %>%
  group_by(Plant) %>%
  summarise(med = median(score),
            lo  = quantile(score, .025),
            hi  = quantile(score, .975),
            p_beyond_D6 = mean(score > cent_draws[, iD6]),
            .groups = "drop")

cat("\n--- Posterior centroid position on the Y6 -> D6 axis ---\n")
cat("0 = Y6 centroid, 1 = D6 centroid, p_beyond_D6 = Pr(score > D6 score)\n\n")
print(cent_summary %>% arrange(med), n = Inf)

# how much of each centroid's displacement is actually on-axis
off_axis <- sapply(seq_along(plants), function(g) {
  d      <- epred[, g, ] - anchor
  onaxis <- (rowSums(d * axis_draws) / rowSums(axis_draws^2))
  resid  <- d - onaxis * axis_draws
  sqrt(rowSums(resid^2)) / sqrt(rowSums(d^2))
})
colnames(off_axis) <- plants
cat("\n--- Median fraction of displacement orthogonal to the parental axis ---\n")
print(round(apply(off_axis, 2, median), 3))

# -------------------------------------------------------------
# 2. Sample-level scores, using the posterior median axis
#    (fixing the axis keeps sample scatter interpretable as
#    among-plant variation rather than axis uncertainty)
# -------------------------------------------------------------
v_med <- apply(axis_draws, 2, median)
a_med <- apply(anchor,     2, median)

sample_scores <- model_df %>%
  select(Plant, Day) %>%
  mutate(score = as.numeric(sweep(X, 2, a_med, "-") %*% v_med / sum(v_med^2)))

# order germplasm by posterior centroid position
lvl <- cent_summary %>% arrange(med) %>% pull(Plant)
sample_scores$Plant <- factor(sample_scores$Plant, levels = lvl)
cent_summary$Plant  <- factor(cent_summary$Plant,  levels = lvl)

pal <- c(Y6 = "#E8B84B", D6 = "#4C7FA8", Y6D6 = "#9BB47A", D6Y6 = "#7FA05C",
         HBR = "#B85C38", RAD = "#8C4A6B")

pA <- ggplot(sample_scores, aes(x = score, y = Plant, colour = Plant)) +
  geom_vline(xintercept = c(0, 1), linetype = c("dashed", "dashed"),
             colour = "grey40", linewidth = .4) +
  geom_jitter(height = .14, width = 0, alpha = .35, size = 1.6, show.legend = FALSE) +
  geom_linerange(data = cent_summary,
                 aes(xmin = lo, xmax = hi, y = Plant),
                 inherit.aes = FALSE, linewidth = .9, colour = "grey15") +
  geom_point(data = cent_summary, aes(x = med, y = Plant),
             inherit.aes = FALSE, size = 3, shape = 21,
             fill = "white", colour = "grey15", stroke = .9) +
  annotate("text", x = 0, y = length(lvl) + .55, label = "Y6 centroid",
           size = 3, colour = "grey30", hjust = .5) +
  annotate("text", x = 1, y = length(lvl) + .55, label = "D6 centroid",
           size = 3, colour = "grey30", hjust = .5) +
  scale_colour_manual(values = pal) +
  coord_cartesian(clip = "off",
                  ylim = c(.5, length(lvl) + .9)) +
  labs(x = "Position on the Y6 \u2192 D6 parental axis (full CLR space)",
       y = NULL, tag = "A") +
  theme_classic(base_size = 11) +
  theme(plot.margin = margin(14, 10, 6, 6))

# -------------------------------------------------------------
# 3. Forest plot of the 15 pairwise Aitchison centroid distances
#    with the D6 to Y6 parental distance as a reference line
# -------------------------------------------------------------
d_parents_summ <- {
  d <- sqrt(rowSums((epred[, iD6, ] - epred[, iY6, ])^2))
  tibble(med = median(d), lo = quantile(d, .025), hi = quantile(d, .975))
}

pw <- pairwise_dist %>%
  mutate(pair = paste(group1, group2, sep = " \u2013 "),
         type = case_when(
           pair %in% c("HBR \u2013 RAD", "RAD \u2013 HBR")     ~ "Convergence pair",
           pair %in% c("D6Y6 \u2013 Y6D6", "Y6D6 \u2013 D6Y6") ~ "Convergence pair",
           (group1 == "D6" & group2 == "Y6") |
             (group1 == "Y6" & group2 == "D6")                   ~ "Parental",
           TRUE                                                ~ "Other"),
         pair = fct_reorder(pair, dist_med))

pB <- ggplot(pw, aes(x = dist_med, y = pair, colour = type)) +
  annotate("rect",
           xmin = d_parents_summ$lo, xmax = d_parents_summ$hi,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = .5) +
  geom_vline(xintercept = d_parents_summ$med, linetype = "dashed",
             colour = "grey35", linewidth = .4) +
  geom_linerange(aes(xmin = dist_lo, xmax = dist_hi), linewidth = .7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("Convergence pair" = "#B85C38",
                                 "Parental"         = "#4C7FA8",
                                 "Other"            = "grey45"),
                      name = NULL) +
  labs(x = "Posterior Aitchison centroid distance (95% CrI)",
       y = NULL, tag = "B") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

fig