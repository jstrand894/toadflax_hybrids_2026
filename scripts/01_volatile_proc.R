# =============================================================================
# Toadflax Hybrid VOC Data: Import, Zero Replacement, and CLR Transformation
# =============================================================================

library(readxl)
library(tidyverse)
library(compositions) 
select    <- dplyr::select
filter    <- dplyr::filter
summarise <- dplyr::summarise


# -----------------------------------------------------------------------------
# 1. Import data
# -----------------------------------------------------------------------------
# The spreadsheet has metadata columns at the front, compound area columns in
# the middle, and a Day column appended at the far right.

raw <- read_excel(
  "data/volatiles/Toadflax_2023_Parent_and_Hybrid_Volatiles.xlsx",
  sheet = "compounds"
)

# Separate metadata, compound data, and the trailing Day column
meta <- bind_cols(
  raw %>% select(`GC File`, ID, Plant, Type),
  raw %>% select(Day = last_col())
)

compounds <- raw %>%
  select(-`GC File`, -ID, -Plant, -Type, -last_col())

# Quick structure check
cat("Samples:     ", nrow(compounds), "\n")
cat("Compounds:   ", ncol(compounds), "\n")
cat("Plant types: ", paste(unique(meta$Plant), collapse = ", "), "\n")
cat("Sample types:", paste(unique(meta$Type),  collapse = ", "), "\n")

raw_df <- bind_cols(meta, as_tibble(compounds_mat))
# -----------------------------------------------------------------------------
# 2. Zero replacement
# -----------------------------------------------------------------------------
# CLR requires strictly positive values. We replace zeros with half the
# smallest observed non-zero value for that compound (detection limit proxy).
# For a more principled approach, swap in zCompositions::cmultRepl().

compounds_mat <- as.matrix(compounds)

replace_zeros <- function(mat) {
  apply(mat, 2, function(col) {
    min_nonzero   <- min(col[col > 0], na.rm = TRUE)
    col[col == 0] <- min_nonzero / 2
    col
  })
}

compounds_nz <- replace_zeros(compounds_mat)


# -----------------------------------------------------------------------------
# 3. CLR transformation
# -----------------------------------------------------------------------------
# For each sample (row): CLR(x_i) = log(x_i / geometric_mean(x))
# This maps the composition to an unconstrained real space suitable for
# standard multivariate methods (PERMANOVA, PCA, etc.).

clr_mat <- t(apply(compounds_nz, 1, function(row) {
  log(row / exp(mean(log(row))))
}))

colnames(clr_mat) <- colnames(compounds)
rownames(clr_mat) <- meta$ID

# Bind CLR values back to metadata for downstream use
clr_df <- bind_cols(meta, as_tibble(clr_mat))


# -----------------------------------------------------------------------------
# 4. Sanity checks
# -----------------------------------------------------------------------------

# CLR rows must sum to ~0 (exact zero up to floating-point error)
row_sums <- apply(clr_mat, 1, sum)
cat("CLR row sum range:", paste(round(range(row_sums), 10), collapse = " to "), "\n")

# Flag any all-zero samples (would indicate missing/failed GC runs)
zero_rows <- apply(compounds_mat, 1, function(r) all(r == 0))
cat("All-zero samples:", sum(zero_rows), "\n")
if (any(zero_rows)) print(meta[zero_rows, ])

# Mean CLR per compound by plant type — useful for spotting gross patterns
clr_df %>%
  pivot_longer(
    cols      = all_of(colnames(compounds)),
    names_to  = "compound",
    values_to = "clr_value"
  ) %>%
  group_by(Plant, compound) %>%
  summarise(mean_clr = mean(clr_value), .groups = "drop") %>%
  pivot_wider(names_from = Plant, values_from = mean_clr) %>%
  print(n = Inf)