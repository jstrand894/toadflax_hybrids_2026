# DT_EADFIDs.R
# Dalmatian Toadflax — EAD-FID traces
# Samples: DT × Janthiniformis, DT × Janthinus

source("helpers_eag_fid.R")


# DT Janthiniformis --------------------------------------------------------

# EAG : data/EAG/janthiniformis_dt_with_spikes.xlsx  (5,999 rows, time 2501–8499)
# FID : data/FID/DT EAGFID 2.10.xlsx                 (13,479 rows)
# Original EAD scale: 0.25 | direction: inverted

dt_janthiniformis <- load_eag_fid(
  eag_path = "../data/EAG/janthiniformis_dt_with_spikes.xlsx",
  fid_path = "../data/FID/DT EAGFID 2.10.xlsx"
)

p_dt_janthiniformis <- plot_eag_fid(
  data      = dt_janthiniformis,
  ead_scale = 0.25,
  ead_flip  = TRUE,
  title     = "DT × Janthiniformis"
)

p_dt_janthiniformis
ggplotly(p_dt_janthiniformis)


# DT Janthinus -------------------------------------------------------------

# EAG : data/EAG/dt_janthinus_with_spikes.xlsx  (5,999 rows, time 5001–10999)
# FID : data/FID/DT_Janthinus.xlsx              (13,810 rows)
# Original EAD scale: 1.0 | direction: inverted

dt_janthinus <- load_eag_fid(
  eag_path = "../data/EAG/dt_janthinus_with_spikes.xlsx",
  fid_path = "../data/FID/DT_Janthinus.xlsx"
)

p_dt_janthinus <- plot_eag_fid(
  data      = dt_janthinus,
  ead_scale = 1.0,
  ead_flip  = TRUE,
  title     = "DT × Janthinus"
)

p_dt_janthinus
ggplotly(p_dt_janthinus)
