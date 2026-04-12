# Y6D6_EADFIDs.R
# Y6D6 hybrid — EAD-FID traces
# Samples: Y6D6 × Janthiniformis, Y6D6 × Janthinus

source("helpers_eag_fid.R")


# Y6D6 Janthiniformis ------------------------------------------------------

# EAG : data/EAG/y6d6_janthiniformis_with_spikes.xlsx
# FID : data/FID/Y6D6_report.xlsx  (311,460 rows)
#
# ⚠ WARNING: y6d6_janthiniformis_with_spikes.xlsx is currently EMPTY (0 rows).
# Uncomment once the EAG file is populated.

# y6d6_janthiniformis <- load_eag_fid(
#   eag_path = "../data/EAG/y6d6_janthiniformis_with_spikes.xlsx",
#   fid_path = "../data/FID/Y6D6_report.xlsx"
# )
#
# p_y6d6_janthiniformis <- plot_eag_fid(
#   data      = y6d6_janthiniformis,
#   ead_scale = 0.25,
#   ead_flip  = TRUE,
#   title     = "Y6D6 × Janthiniformis"
# )
#
# p_y6d6_janthiniformis
# ggplotly(p_y6d6_janthiniformis)


# Y6D6 Janthinus -----------------------------------------------------------

# No EAG file exists yet for Y6D6 Janthinus.
