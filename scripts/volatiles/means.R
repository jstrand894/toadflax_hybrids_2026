source("scripts/volatiles/volatile_proc.r")

# ---- 4. Mean raw emissions by Plant type ----
mean_raw_by_plant <- compounds_mat %>%
  as_tibble() %>%
  bind_cols(meta %>% select(Plant), .) %>%
  group_by(Plant) %>%
  reframe(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = -Plant,
               names_to = "compound",
               values_to = "mean_raw") %>%
  pivot_wider(names_from = Plant,
              values_from = mean_raw)

print(mean_raw_by_plant, n = Inf)
