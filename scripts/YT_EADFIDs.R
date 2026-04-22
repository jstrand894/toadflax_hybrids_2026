# YT_EADFIDs.R
# Yellow Toadflax — EAD-FID traces
# Samples: YT × Janthiniformis, YT × Janthinus

source("scripts/helpers_eag_fid.R")

# YT Janthiniformis --------------------------------------------------------

# EAG : data/EAG/yt_janthiniformis_with_spikes.xlsx  (4,699 rows, time 3301–7999)
# FID : data/FID/YT_Janthiniformis.xlsx              (11,735 rows)
# Original EAD scale: 0.45 | direction: NOT inverted
# Note: spike_values were negated in original code, so the EAG file already
#       encodes downward spike direction — ead_flip = FALSE preserves this.

yt_janthiniformis <- load_eag_fid(
  eag_path = "data/EAG/yt_janthiniformis_with_spikes.xlsx",
  fid_path = "data/FID/YT_Janthiniformis.xlsx"
)

p_yt_janthiniformis <- plot_eag_fid(
  data      = yt_janthiniformis,
  ead_scale = 0.45,
  ead_flip  = FALSE,
  title     = "YT × Janthiniformis"
)

p_yt_janthiniformis
ggplotly(p_yt_janthiniformis)


# YT Janthinus -------------------------------------------------------------

# EAG : data/EAG/yt_janthinus_with_spikes.xlsx  (2,249 rows, time 3751–5999)
# FID : data/FID/YT_Janthinus.xlsx              (11,393 rows)
# Original EAD scale: 0.35 | direction: inverted
#
# BUG FIX: original plot used yt_janthinus_spikes for BOTH geom_lines.
#           FID line now correctly uses the joined data.

yt_janthinus <- load_eag_fid(
  eag_path = "../data/EAG/yt_janthinus_with_spikes.xlsx",
  fid_path = "../data/FID/YT_Janthinus.xlsx"
)

p_yt_janthinus <- plot_eag_fid(
  data      = yt_janthinus,
  ead_scale = 0.35,
  ead_flip  = TRUE,
  title     = "YT × Janthinus"
)

p_yt_janthinus
ggplotly(p_yt_janthinus)
