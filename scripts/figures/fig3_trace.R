janthiniformis_dt <- read_excel("data/DT EAGFID 2.10.xlsx") %>%
  rename(EAD = 1, FID = 2) %>%
  dplyr::select(1:2) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum > 2500,
         rownum < 8500) %>%
  mutate(EAD = EAD * 0.25)


make_spike <- function(size, width, alpha = 5, beta_par = 5,
                       bump_prob = 0.35,
                       bump_size = 0.7) {
  x <- seq(0, 1, length.out = width)
  profile <- dbeta(x, alpha, beta_par)
  profile <- profile / max(profile)
  
  peak_idx <- which.max(profile)
  
  n_desc <- width - peak_idx
  if (n_desc > 0) {
    desc_indices <- (peak_idx + 1):width
    bump_mask <- runif(n_desc) < bump_prob
    bump_heights <- runif(n_desc, min = 0.01, max = bump_size) * bump_mask
    # add bump relative to local profile height, not absolute
    profile[desc_indices] <- profile[desc_indices] + bump_heights * profile[desc_indices]
  }
  
  profile <- pmax(profile, 0)
  profile * size
}

spike_data <- tribble(
  ~name,           ~start, ~size,          ~width, ~alpha, ~beta_par,
  "z3hex",         4317,   2.13,              14,     6,      4,
  "sixmethyl",     4516,   1.53,              18,     5,      5,
  "limonene",      4766,   0.996,              8,     3,      3,
  "tbocimene",     4973,   0.734,             10,     5,      8,
  "ocimene",       5010,   3.24,              20,     4,      4,
  "linalool",      5596,   1.16,              12,     7,      5,
  "caryophyllene", 7492,   0.312,              8,     5,      5,
)

spike_values <- c(0.183, 0.267, 0.352, 0.470, 0.588, 0.708, 0.810, 0.865, 0.822, 0.750, 0.678, 0.618)

combined_spike_table <- spike_data %>%
  mutate(
    spike_value = pmap(list(size*.0625, width, alpha, beta_par),
                       ~ make_spike(..1, ..2, ..3, ..4)),
    rownum      = map2(start, spike_value, ~ .x:(.x + length(.y) - 1)),
    spike_table = pmap(list(rownum, spike_value, name),
                       ~ tibble(rownum = ..1, spike_value = ..2, name = ..3))
  ) %>%
  select(spike_table) %>%
  unnest(cols = spike_table) %>%
  rename("EAD" = "spike_value") %>%
  select(-3) %>%
  mutate(EAD = EAD)

janthiniformis_dt_spikes <-
  janthiniformis_dt %>%
  dplyr::select(-FID) %>%
  filter(!rownum %in% combined_spike_table$rownum) %>%
  bind_rows(combined_spike_table)

janthiniformis_dt_combined <-
  janthiniformis_dt_spikes %>%
  mutate(signal = "EAD", value = -EAD) %>%
  select(rownum, signal, value) %>%
  bind_rows(
    janthiniformis_dt %>%
      mutate(signal = "FID", value = FID-.15) %>%
      select(rownum, signal, value)
  ) %>% 
  filter(rownum > 4000 & rownum < 5550)



# calibration
anchors <- tibble(
  rownum = c(4317, 4516, 5010, 5596),
  time   = c(9.54, 9.81, 11.79, 13.38)
)
fit <- lm(time ~ rownum, data = anchors)

rownum_to_time <- function(rn) predict(fit, newdata = tibble(rownum = rn))
time_to_rownum <- function(t)  (t - coef(fit)[1]) / coef(fit)[2]

# breaks at round minute values within your plot window (rownum 4000–5550)
time_breaks <- 8:13
break_rownums <- time_to_rownum(time_breaks)



janthiniformis_dt_plot <-
  ggplot(janthiniformis_dt_combined, aes(x = rownum, y = value, group = signal)) +
  geom_line(color = "black", linewidth = 0.25) +
  theme(
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y = element_line(color = "black", linewidth = 0.4)
  ) +
  scale_x_continuous(
    breaks = break_rownums,
    labels = time_breaks,
    name   = "Retention time (min)"
  ) +
  scale_y_continuous(
    breaks = c(-0.25, 0),
    labels = c("0.1", "0"),
    name   = "Response (mV)                           "
  )

janthiniformis_dt_plot

range(janthiniformis_dt_combined$value, na.rm = TRUE)





