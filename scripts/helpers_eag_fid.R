# helpers_eag_fid.R
# Shared helper functions for EAD-FID loading and plotting.
# Source this file at the top of each plant-type script.
#
# CHANGE LOG (vs. original All EADFIDs.Rmd):
#   - EAG data now read from data/EAG/ (pre-processed, spikes already embedded)
#   - FID data now read from data/FID/
#   - The `FID` column in EAG files is actually the EAD signal (instrument
#     mislabeling) — renamed to EAD on load
#   - rownum added to FID files via row_number() (new files lack this column)
#   - Trimming is implicit: inner_join on time/rownum replaces manual
#     filter(rownum > X, rownum < Y) blocks
#   - All spike_data tribbles and combined_spike_table construction removed
#   - mutate(EAD*0.5) no-op bugs removed
#   - <<noweb>> chunk references removed (invalid in standard Rmd)

library(tidyverse)
library(readxl)
library(plotly)
library(cowplot)


# ── load_eag_fid() ────────────────────────────────────────────────────────────
# Reads one pre-processed EAG file and one FID file, inner-joins on the shared
# time/rownum index. The join naturally trims the FID to the EAG window.
#
# Args:
#   eag_path : path to EAG xlsx  — columns: time (int), FID* (float)
#   fid_path : path to FID xlsx  — columns: EAD AVE, FID AVE, EAD 1, FID 1
#
# * EAG "FID" column is the EAD/depolarization signal (instrument mislabeling).
#
# Returns: tibble with columns: time, EAD, FID
# ─────────────────────────────────────────────────────────────────────────────

load_eag_fid <- function(eag_path, fid_path) {

  eag <- read_excel(eag_path) %>%
    rename(time = 1, EAD = 2)          # col 2 labeled "FID" but is EAD signal

  fid <- read_excel(fid_path) %>%
    mutate(rownum = row_number()) %>%   # new FID files lack a rownum column
    select(rownum, FID = `FID AVE`)

  eag %>%
    inner_join(fid, by = c("time" = "rownum"))
}


# ── plot_eag_fid() ────────────────────────────────────────────────────────────
# Dual-trace plot: FID (black) + scaled EAD (blue), clean chromatogram style.
#
# Args:
#   data      : tibble from load_eag_fid()
#   ead_scale : visual scale multiplier for EAD trace
#   ead_flip  : if TRUE, EAD plotted as negative (trace drops below baseline)
#   title     : optional plot title string
# ─────────────────────────────────────────────────────────────────────────────

plot_eag_fid <- function(data, ead_scale = 1, ead_flip = TRUE, title = NULL) {

  ead_sign <- if (ead_flip) -1 else 1

  p <- ggplot(data) +
    geom_line(aes(x = time, y = FID),
              linewidth = 0.4) +
    geom_line(aes(x = time, y = ead_sign * EAD * ead_scale),
              color = "steelblue", linewidth = 0.4) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      panel.grid       = element_blank()
    )

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}
