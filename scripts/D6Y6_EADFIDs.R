# D6Y6_EADFIDs.R
# D6Y6 hybrid — EAD-FID traces
# Samples: D6Y6 × Janthiniformis, D6Y6 × Janthinus
#
# Both samples share the same FID file (D6Y6_report.xlsx, 145,060 rows).
# The inner_join windows each to its respective EAG time range.
#
# Note: d6y6_janthiniformis and d6y6_janthinus EAG files share the same time
# range (2281–4656) and had identical opening rows on inspection. If these
# represent two antennae from the same FID run, EAD signals should diverge
# partway through — verify visually.
#
# BUG FIX: <<D6Y6_FID>> noweb chunk references removed (invalid R syntax).
# BUG FIX: both original plots used _spikes for both geom_lines; FID line
#           now correctly uses the joined data.

source("scripts/helpers_eag_fid.R")
source(here("scripts/helpers_eag_fid.R"))
# D6Y6 Janthiniformis ------------------------------------------------------

# EAG : data/EAG/d6y6_janthiniformis_with_spikes.xlsx  (2,376 rows, time 2281–4656)
# FID : data/FID/D6Y6_report.xlsx                      (145,060 rows)
# Original EAD scale: ~0.125 (size*0.5 in spike table, *0.25 in plot)

d6y6_janthiniformis <- load_eag_fid(
  eag_path = "data/EAG/d6y6_janthiniformis_with_spikes.xlsx",
  fid_path = "data/FID/D6Y6_report.xlsx"
)

p_d6y6_janthiniformis <- plot_eag_fid(
  data      = d6y6_janthiniformis,
  ead_scale = 0.25,   # adjust if pre-processed EAG is already scaled
  ead_flip  = TRUE,
  title     = "D6Y6 × Janthiniformis"
)

p_d6y6_janthiniformis
ggplotly(p_d6y6_janthiniformis)


# D6Y6 Janthinus -----------------------------------------------------------

# EAG : data/EAG/d6y6_janthinus_with_spikes.xlsx  (2,376 rows, time 2281–4656)
# FID : data/FID/D6Y6_report.xlsx                 (same file as above)

d6y6_janthinus <- load_eag_fid(
  eag_path = "data/EAG/d6y6_janthinus_with_spikes.xlsx",
  fid_path = "data/FID/D6Y6_report.xlsx"
)

p_d6y6_janthinus <- plot_eag_fid(
  data      = d6y6_janthinus,
  ead_scale = 0.25,
  ead_flip  = TRUE,
  title     = "D6Y6 × Janthinus"
)

p_d6y6_janthinus
ggplotly(p_d6y6_janthinus)
