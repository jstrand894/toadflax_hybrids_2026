library(readxl)
library(tidyverse)
library(compositions)   # for clr()
# optional but preferred for zero imputation:
# library(zCompositions)

# ---- 1. Import the compounds tab ----
raw <- read_excel(
  "data/volatiles/Toadflax_2023_Parent_and_Hybrid_Volatiles.xlsx",
  sheet = "compounds"
)

# The Day column is appended at the end; move it into metadata
# Metadata cols: GC File, ID, Plant, Type, Day (last col)
# Compound cols: everything in between

meta_front <- raw %>% select(`GC File`, ID, Plant, Type)
day_col    <- raw %>% select(Day = last_col())
compounds  <- raw %>% select(-`GC File`, -ID, -Plant, -Type, -last_col())

meta <- bind_cols(meta_front, day_col)

# Confirm structure
cat("Samples:", nrow(compounds), "\n")
cat("Compounds:", ncol(compounds), "\n")
cat("Plant types:", paste(unique(meta$Plant), collapse = ", "), "\n")
cat("Sample types:", paste(unique(meta$Type), collapse = ", "), "\n")


# ---- 2. CLR transformation ----
# Replace zeros with half the minimum non-zero value (per compound)
# This is simple and defensible; swap for zCompositions::cmultRepl() if you
# want a more rigorous Bayesian multiplicative replacement

compounds_mat <- as.matrix(compounds)

zero_replace <- function(mat) {
  apply(mat, 2, function(col) {
    min_nz <- min(col[col > 0], na.rm = TRUE)
    col[col == 0] <- min_nz / 2
    col
  })
}

compounds_nz <- zero_replace(compounds_mat)

# CLR: log(x_i / geometric_mean(x)) for each sample row
clr_mat <- t(apply(compounds_nz, 1, function(row) {
  log(row / exp(mean(log(row))))
}))

colnames(clr_mat) <- colnames(compounds)
rownames(clr_mat) <- meta$ID

# Or equivalently using the compositions package:
# clr_mat <- clr(compounds_nz)   # returns a 'clr' class object
# clr_mat <- unclass(clr_mat)    # strip class for use in vegan etc.

clr_df <- as_tibble(clr_mat) %>%
  bind_cols(meta, .)


# ---- 3. Quick sanity checks ----
# CLR rows should sum to ~0 (floating point noise only)
row_sums <- apply(clr_mat, 1, sum)
cat("CLR row sum range:", round(range(row_sums), 10), "\n")  # should be ~0

# Check compound distributions by plant type
clr_long <- clr_df %>%
  pivot_longer(cols = all_of(colnames(compounds)),
               names_to = "compound",
               values_to = "clr_value")

clr_long %>%
  group_by(Plant, compound) %>%
  summarise(mean_clr = mean(clr_value), .groups = "drop") %>%
  pivot_wider(names_from = Plant, values_from = mean_clr) %>%
  print(n = 26)




zero_rows <- apply(compounds_mat, 1, function(r) all(r == 0))
cat("All-zero samples:", sum(zero_rows), "\n")
meta[zero_rows, ]
