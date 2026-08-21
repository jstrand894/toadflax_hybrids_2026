library(dplyr); library(tidyr)

supp <- out %>%
  mutate(Compound = factor(unname(lab[compound]), levels = ord),
         cell = sprintf("%.2f ± %.2f", mean, se)) %>%
  select(Compound, Plant, cell) %>%
  pivot_wider(names_from = Plant, values_from = cell) %>%
  arrange(Compound) %>%
  select(Compound, D6, Y6, D6Y6, Y6D6, HBR, RAD)

write.csv(supp, "table_S_ead_blend.csv", row.names = FALSE)
write.table(supp, sep = "\t", row.names = FALSE, quote = FALSE)