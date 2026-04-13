# Toadflax Hybrids 2026

**Volatile organic compound profiles of hybrid toadflax (*Linaria* spp.) and olfactory host preference in two specialist biocontrol weevils (*Mecinus janthinus* and *M. janthiniformis*)**

J. R. Strand, S. E. Sing, E. J. Friesenhahn, M. J. Hofland, D. K. Weaver

---

## Overview

Yellow toadflax (*Linaria vulgaris*) and Dalmatian toadflax (*L. dalmatica*) are highly invasive plants across western North America. Two cryptic weevil species — *Mecinus janthinus* (specialist on Dalmatian toadflax) and *M. janthiniformis* (specialist on yellow toadflax) — have been introduced as biological control agents. Natural hybridization between the two toadflax species produces plants with volatile profiles that differ from either parent, potentially placing hybrids outside the host-recognition window of both weevils and creating a refuge from biocontrol pressure.

This project integrates three experimental approaches across six toadflax germplasm types (two parental species, two synthetic F1 hybrids, two wild hybrid populations):

- **GC-MS** — volatile organic compound (VOC) profiling of plant headspace
- **EAD-FID** — electroantennography coupled with flame ionization detection to identify electrophysiologically active compounds
- **Y-tube olfactometry** — behavioral choice assays measuring weevil host preference

---

## Germplasm Types

| Code | Description |
|------|-------------|
| YT | Yellow toadflax (*L. vulgaris*) — natal host of *M. janthiniformis* |
| DT | Dalmatian toadflax (*L. dalmatica*) — natal host of *M. janthinus* |
| Y6D6 | Synthetic F1 hybrid (YT × DT) |
| D6Y6 | Synthetic F1 hybrid (DT × YT) |
| HBR | Wild hybrid population (Boulder River) |
| RAD | Wild hybrid population (Radersburg) |

---

## Repository Structure

```
toadflax_hybrids_2026/
│
├── data/
│   ├── GCMS/                        # Raw and processed GC-MS volatile data
│   │   ├── all_runs.xlsx            # Combined data across all GC-MS runs
│   │   └── GCMS Runs/               # Individual run files by germplasm type
│   │
│   ├── EAG/                         # EAD-FID electrophysiology data
│   │   ├── depolarizations.xlsx     # Compiled antennal depolarization events
│   │   └── *_with_spikes.xlsx       # Per-run files with spike annotations
│   │
│   ├── FID/                         # Flame ionization detector outputs
│   │   └── *_report.xlsx            # Per-germplasm FID reports
│   │
│   ├── volatiles/                   # Processed volatile compound tables
│   │   └── Toadflax_2023_Parent_and_Hybrid_Volatiles.xlsx
│   │
│   ├── ytubes/                      # Y-tube behavioral assay results
│   │   ├── Janthiniformis Y-tubes 2023.xlsx
│   │   ├── Janthinus Y-tubes 2023.xlsx
│   │   └── ytubes_results.xlsx
│   │
│   └── Janthiniformis EAG/          # Legacy/working EAD-FID files and .Rmd notebooks
│
├── scripts/
│   ├── helpers_eag_fid.R            # Shared helper functions for EAD-FID processing
│   ├── extract.R                    # Data extraction utilities
│   ├── *_EADFIDs.R                  # Per-germplasm EAD-FID analysis scripts
│   ├── EAD/
│   │   └── depolarizations.R        # Depolarization event processing
│   └── volatiles/
│       ├── volatile_proc.R          # GC-MS data processing pipeline
│       ├── volatile_ordination.R    # Multivariate ordination (PCA/NMDS)
│       └── means.R                  # Summary statistics for volatile compounds
│
├── fid_ead_exports/                 # Cleaned, export-ready EAD-FID datasets
│   └── *_FID_EAD.xlsx               # One file per species × germplasm combination
│
├── Hybrids_mecinus_behavior_MS.pdf  # Manuscript draft
├── JRS_NFBC2025_4-9-25.pdf          # Conference presentation (NFBC 2025)
└── toadflax_hybrids_2026.Rproj      # RStudio project file
```

---

## Analysis Pipeline

1. **GC-MS processing** (`scripts/volatiles/volatile_proc.R`) — combine raw runs, align retention times, and identify compounds
2. **Volatile ordination** (`scripts/volatiles/volatile_ordination.R`) — PCA/NMDS to compare volatile profiles across germplasm types
3. **EAD-FID processing** (`scripts/*_EADFIDs.R`) — extract depolarization events and match to FID peaks
4. **Depolarization summary** (`scripts/EAD/depolarizations.R`) — compile antennal responses across compounds and species
5. **Y-tube analysis** (`scripts/ytubes.R`) — test for significant host preference using binomial/proportion tests

---

## Dependencies

All analyses are run in **R**. Key packages used include `tidyverse`, `vegan` (ordination), `readxl`/`writexl` (data I/O), and `ggplot2` (visualization).
