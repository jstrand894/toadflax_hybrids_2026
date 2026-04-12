# Export datasets to Excel files
library(writexl)

# Create a function to prepare data for export
prepare_for_export <- function(data, dataset_name) {
  # Get the original FID data and merge with spikes
  result <- data %>%
    select(rownum, FID) %>%
    mutate(time = rownum) %>%  # You can adjust time conversion if needed
    select(time, FID) %>%
    filter(!is.na(FID))

  # Write to Excel
  write_xlsx(result, paste0(dataset_name, ".xlsx"))
  cat("Exported:", dataset_name, ".xlsx\n")
}

# Export each dataset
prepare_for_export(janthiniformis_dt_spikes %>%
                   left_join(janthiniformis_dt %>% select(rownum, FID), by = "rownum"),
                   "janthiniformis_dt_with_spikes")

prepare_for_export(dt_janthinus_spikes %>%
                   left_join(dt_janthinus %>% select(rownum, FID), by = "rownum"),
                   "dt_janthinus_with_spikes")

prepare_for_export(yt_janthiniformis_spikes %>%
                   left_join(yt_janthiniformis %>% select(rownum, FID), by = "rownum"),
                   "yt_janthiniformis_with_spikes")

prepare_for_export(yt_janthinus_spikes, "yt_janthinus_with_spikes")

prepare_for_export(d6y6_janthiniformis_spikes, "d6y6_janthiniformis_with_spikes")

prepare_for_export(d6y6_janthinus_spikes, "d6y6_janthinus_with_spikes")

prepare_for_export(y6d6_janthiniformis_spikes, "y6d6_janthiniformis_with_spikes")

prepare_for_export(fid_line_spikes, "baseline_with_spikes")
