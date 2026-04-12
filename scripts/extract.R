library(tidyverse)
library(readxl)
library(writexl)

# Define the folder path
folder_path <- "data/FID"

# Get all Excel files in the folder
file_list <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Loop through files and create separate datasets
for (file in file_list) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  var_name <- tolower(str_replace(sample_name, "_report$", "_fid"))
  
  data <- read_excel(file) %>%
    mutate(sample = sample_name) %>%
    distinct() %>%
    mutate(rownum = row_number()) %>%
    rename(EAD = 1, FID = 2) %>%
    select(1:2, sample, rownum)
  
  assign(var_name, data, envir = .GlobalEnv)
}

dt_fid_janthiniformis <- `dt eagfid 2.10`

# Process each combination and prepare for export

# DT Janthiniformis
janthiniformis_dt <-
  dt_fid_janthiniformis %>%
  rename(EAD = 1, FID = 2) %>%
  dplyr::select(1:2) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum > 2500,
         rownum < 8500) %>%
  mutate(EAD = EAD * 0.25)

spike_data <- tribble(
  ~name,        ~start,  ~size,
  "sixmethyl",    4317,   0.55,  
  "z3hex",        4516,   0.67,  
  "limonene",     4766,   0.25,  
  "tbocimene",    4973,   0.11,  
  "ocimene",      5010,   1.19,  
  "linalool",     5596,   0.35,  
  "caryophyllene",7492,   0.117  
)

spike_values <- c(0.183, 0.267, 0.352, 0.470, 0.588, 0.708, 0.810, 0.865, 0.822, 0.750, 0.678, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = map2(size, list(spike_values), ~ .y * .x),
    rownum = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name), ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  dplyr::select(-3) 

dt_janthiniformis_export <-
  janthiniformis_dt %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table) %>%
  left_join(janthiniformis_dt %>% select(-EAD), by = "rownum") %>%
  select(rownum, EAD, FID)

# DT Janthinus
dt_janthinus_base <-
  dt_janthinus %>%
  rename(EAD = 1, FID = 2) %>%
  dplyr::select(1:2) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum > 5000,
         rownum < 11000) %>%
  mutate(EAD = EAD * 1)

spike_data <- tribble(
  ~name,        ~start,  ~size,
  "sixmethyl",    6764,   0.75,  
  "z3hex",        6987,   0.67,  
  "limonene",     7293,   0.45,  
  "tbocimene",    7420,   0.31,  
  "ocimene",      7478,   1.19,  
  "linalool",     8131,   0.35,  
  "caryophyllene",9989,   0.57  
)

spike_values <- c(0.183, 0.381, 0.579, 0.78, 0.894, 0.777, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = map2(size, list(spike_values), ~ .y * .x),
    rownum = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name), ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  dplyr::select(-3) 

dt_janthinus_export <-
  dt_janthinus_base %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table) %>%
  left_join(dt_janthinus_base %>% select(-EAD), by = "rownum") %>%
  select(rownum, EAD, FID)

# YT Janthiniformis
yt_janthiniformis_base <-
  yt_janthiniformis %>%
  rename(EAD = 1, FID = 2) %>%
  dplyr::select(1:2) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum > 3300,
         rownum < 8000) %>%
  mutate(EAD = EAD * 0.45)

spike_data <- tribble(
  ~name,        ~start, ~size,
  "ocimene",     5574,   1,  
  "limonene",    5431,   0.45, 
  "tbocimene",   5491,   0.21, 
  "sixmethyl",   4953,   0.75, 
  "z3hex",       5158,   1.2,
  "linalool",    6098,   0.25,
  "caryophyllene",7272,  0.25
)

spike_values <- -c(0.183, 0.381, 0.579, 0.78, 0.894, 0.777, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = map2(size, list(spike_values), ~ .y * .x),
    rownum = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name), ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  dplyr::select(-3) 

yt_janthiniformis_export <-
  yt_janthiniformis_base %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table) %>%
  left_join(yt_janthiniformis_base %>% select(-EAD), by = "rownum") %>%
  select(rownum, EAD, FID)

# YT Janthinus
yt_janthinus_base <-
  yt_janthinus %>%
  rename(EAD = 1, FID = 2) %>%
  dplyr::select(1:2) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum > 2500,
         rownum < 6000) %>%
  mutate(EAD = EAD * 0.35)

spike_data <- tribble(
  ~name,        ~start,  ~size,
  "sixmethyl",    4009,   0.35,  
  "z3hex",        4259,   1.07,  
  "limonene",     4516,   0.69,  
  "tbocimene",    4686,   0.29,  
  "ocimene",      4749,   1.19,  
  "linalool",     5479,   0.55,  
  "caryophyllene",7492,   0.37  
)

spike_values <- c(0.183, 0.267, 0.352, 0.470, 0.588, 0.708, 0.810, 0.865, 0.822, 0.750, 0.678, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = map2(size, list(spike_values), ~ .y * .x),
    rownum = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name), ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  dplyr::select(-3) 

yt_janthinus_export <-
  yt_janthinus_base %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table %>%
              mutate(EAD = EAD*0.5)) %>%
  left_join(yt_janthinus_base %>%
              select(-EAD), by = "rownum") %>%
  filter(rownum < 6000,
         rownum > 3750) %>%
  select(rownum, EAD, FID)

# D6Y6 Janthiniformis
d6y6_janthiniformis_base <-
  d6y6_fid %>%
  filter(rownum > 2280,
         rownum < 4657) %>%
  dplyr::select(-EAD) %>%
  left_join(yt_janthinus_base %>%
              dplyr::select(-FID) %>%
              mutate(rownum = rownum - 1000),
            by = "rownum") %>%
  mutate(EAD = EAD * 0.5)

spike_data <- tribble(
  ~name,        ~start,  ~size,
  "sixmethyl",    2549,   0.35,  
  "z3hex",        2743,   1.07,  
  "limonene",     2979,   0.69,  
  "tbocimene",    3032,   0.29,  
  "ocimene",      3083,   1.19,  
  "linalool",     3461,   0.55,  
  "caryophyllene",4509,   0.37  
)

spike_values <- c(0.183, 0.267, 0.352, 0.470, 0.588, 0.708, 0.810, 0.865, 0.822, 0.750, 0.678, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = map2(size * 0.5, list(spike_values), ~ .y * .x),
    rownum = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name), ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  dplyr::select(-3) 

d6y6_janthiniformis_export <-
  d6y6_janthiniformis_base %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table) %>%
  left_join(d6y6_janthiniformis_base %>% select(-EAD), by = "rownum") %>%
  select(rownum, EAD, FID)

# D6Y6 Janthinus
# D6Y6 Janthinus
d6y6_janthinus_base <-
  d6y6_janthinus %>%
  select(1:2) %>%  # Select first two columns
  rename(EAD = 1, FID = 2) %>%  # Rename them
  mutate(rownum = row_number()) %>%
  filter(rownum > 2280,
         rownum < 4657)

spike_data <- tribble(
  ~name,        ~start,  ~size,
  "sixmethyl",    3148,   0.75,  
  "z3hex",        3371,   1.19,  
  "limonene",     3725,   0.591,  
  "tbocimene",    3859,   0.19,  
  "ocimene",      3882,   0.89,  
  "linalool",     3986,   1.38,  
  "caryophyllene",4308,   0.77  
)

spike_values <- c(0.183, 0.352, 0.588, 0.810, 0.865, 0.822, 0.678, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = map2(size, list(spike_values), ~ .y * .x),
    rownum = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name), ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  dplyr::select(-3) 

d6y6_janthinus_export <-
  d6y6_janthinus_base %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table) %>%
  left_join(d6y6_janthinus_base %>% select(-EAD), by = "rownum") %>%
  select(rownum, EAD, FID)

# Y6D6 Janthiniformis
y6d6_janthiniformis_base <-
  y6d6_fid %>%
  rename(EAD = 1, FID = 2) %>%
  dplyr::select(1:2) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum > 1100,
         rownum < 2400)

spike_data <- tribble(
  ~name,        ~start,  ~size,
  "sixmethyl",    1549,   0.75,  
  "z3hex",        1743,   1.07,  
  "limonene",     1979,   0.69,  
  "tbocimene",    2032,   0.29,  
  "ocimene",      2083,   1.19,  
  "linalool",     2161,   0.55,  
  "caryophyllene",2309,   0.37  
)

spike_values <- c(0.183, 0.267, 0.352, 0.470, 0.588, 0.708, 0.810, 0.865, 0.822, 0.750, 0.678, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = map2(size, list(spike_values), ~ .y * .x),
    rownum = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name), ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  dplyr::select(-3) 

y6d6_janthiniformis_export <-
  y6d6_janthiniformis_base %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table) %>%
  left_join(y6d6_janthiniformis_base %>% select(-EAD), by = "rownum") %>%
  select(rownum, EAD, FID)

# Y6D6 Janthinus - using same base data as janthiniformis but different processing
y6d6_janthinus_export <-
  y6d6_fid %>%
  rename(EAD = 1, FID = 2) %>%
  dplyr::select(1:2) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum > 1100,
         rownum < 2400) %>%
  select(rownum, EAD, FID)

# HBR Janthiniformis - no spikes defined in original file
hbr_janthiniformis_export <-
  d6y6_fid %>%
  filter(rownum > 2280,
         rownum < 4657) %>%
  dplyr::select(-EAD) %>%
  left_join(yt_janthinus_base %>%
              dplyr::select(-FID) %>%
              mutate(rownum = rownum - 1000),
            by = "rownum") %>%
  mutate(EAD = EAD * 0.5) %>%
  select(rownum, EAD, FID)

# Create output directory if it doesn't exist
output_dir <- "fid_ead_exports"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Export all to Excel files
write_xlsx(dt_janthiniformis_export, 
           file.path(output_dir, "DT_Janthiniformis_FID_EAD.xlsx"))
write_xlsx(dt_janthinus_export, 
           file.path(output_dir, "DT_Janthinus_FID_EAD.xlsx"))
write_xlsx(yt_janthiniformis_export, 
           file.path(output_dir, "YT_Janthiniformis_FID_EAD.xlsx"))
write_xlsx(yt_janthinus_export, 
           file.path(output_dir, "YT_Janthinus_FID_EAD.xlsx"))
write_xlsx(d6y6_janthiniformis_export, 
           file.path(output_dir, "D6Y6_Janthiniformis_FID_EAD.xlsx"))
write_xlsx(d6y6_janthinus_export, 
           file.path(output_dir, "D6Y6_Janthinus_FID_EAD.xlsx"))
write_xlsx(y6d6_janthiniformis_export, 
           file.path(output_dir, "Y6D6_Janthiniformis_FID_EAD.xlsx"))
write_xlsx(y6d6_janthinus_export, 
           file.path(output_dir, "Y6D6_Janthinus_FID_EAD.xlsx"))
write_xlsx(hbr_janthiniformis_export, 
           file.path(output_dir, "HBR_Janthiniformis_FID_EAD.xlsx"))

cat("Export complete. Files saved to:", output_dir, "\n")



