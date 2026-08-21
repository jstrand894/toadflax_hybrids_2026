# =============================================================
# ytube_raw_counts_table.R
# Table X. Raw plant-arm and blank-arm counts per cell
# =============================================================

library(dplyr)
library(tidyr)
library(readxl)

df_raw <- read_excel("data/ytubes/ytubes_results.xlsx", sheet = "all") %>%
  rename(
    season    = Season,
    date      = Date,
    sex       = Sex,
    pos_trt   = `Pos Trt`,
    neg_trt   = `Neg Trt`,
    pos_count = `# Pos`,
    neg_count = `# Neg`,
    nr_count  = `# NR`
  ) %>%
  dplyr::select(species, season, date, sex, pos_trt, neg_trt,
                pos_count, neg_count, nr_count) %>%
  filter(!is.na(species)) %>%
  mutate(
    species = factor(species),
    sex     = factor(tolower(sex)),
    pos_trt = factor(pos_trt, levels = c("YT6", "DT6", "Y6D6", "D6Y6", "RAD", "HBR"))
  )

check <- df_raw %>%
  filter(sex == "female") %>%
  mutate(total = pos_count + neg_count) %>%
  filter(total > 0) %>%
  group_by(species, pos_trt) %>%
  summarise(model_responders = sum(total), .groups = "drop")



tableX <- df_raw %>%
  group_by(species, sex, pos_trt) %>%
  summarise(
    n_trials   = n(),
    n_days     = n_distinct(date),
    plant_arm  = sum(pos_count, na.rm = TRUE),
    blank_arm  = sum(neg_count, na.rm = TRUE),
    nonresp    = sum(nr_count,  na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(
    n_responders = plant_arm + blank_arm,
    n_tested     = n_responders + nonresp,
    prop_plant   = round(plant_arm / n_responders, 3),
    nr_rate      = round(nonresp / n_tested, 3)
  ) %>%
  arrange(species, sex, pos_trt)

tableX %>%
  filter(sex == "female") %>%
  select(species, pos_trt, n_responders) %>%
  left_join(check, by = c("species", "pos_trt")) %>%
  mutate(match = n_responders == model_responders) %>%
  print(n = Inf)

tableX %>%
  group_by(species, sex) %>%
  summarise(
    nonresp  = sum(nonresp),
    n_tested = sum(n_tested),
    nr_rate  = round(nonresp / n_tested, 3),
    .groups  = "drop"
  ) %>%
  print(n = Inf)

cat("\nTable X. Raw y-tube choice counts per cell.\n")
cat("Denominator for prop_plant is responders only (plant_arm + blank_arm).\n\n")
print(as.data.frame(tableX), row.names = FALSE)

nr_by_species <- tableX %>%
  group_by(species, sex) %>%
  summarise(
    nonresp  = sum(nonresp),
    n_tested = sum(n_tested),
    nr_rate  = round(nonresp / n_tested, 3),
    .groups  = "drop"
  )

model_check <- tableX %>%
  filter(sex == "female") %>%
  select(species, pos_trt, n_responders) %>%
  left_join(check, by = c("species", "pos_trt")) %>%
  mutate(match = n_responders == model_responders)


library(writexl)

dir.create("outputs", showWarnings = FALSE)

write_xlsx(
  list(
    `Table X` = tableX,
    `NR by species` = nr_by_species,
    `Model check` = model_check
  ),
  path = "figures/suptab1_ytube_raw_counts.xlsx"
)
