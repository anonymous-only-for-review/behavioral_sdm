# R Code and data for the MS: "Behavioural thermoregulation and microclimate reshape climate-driven range forecasts"

Code and data supporting analyses of behavioural buffering and hybrid species distribution models under climate warming.
### Abstract
Microclimate heterogeneity and behavioural thermoregulation can buffer ectotherms against environmental warming, yet these processes are rarely captured in models forecasting species’ responses to climate change. Here we use microclimate and biophysical models validated against field data separated by 25 years of warming (1997-2022) to quantify how behavioural thermoregulation modifies body temperature in a widespread lizard species, and how this influences its geographic distribution. We show that body temperatures increased by only 18-29% of the rise in environmental temperatures, showing the role of behavioural buffering. When incorporated into species distribution models, this buffering mechanism emerged as a key predictor of current distributions and led to markedly different projections under future warming compared to conventional climate-based models. These differences revealed cryptic refugia in cooler regions, where behavioural buffering reduces climate impacts, and false refugia in warmer regions, where buffering capacity is exceeded and persistence is overestimated by correlative approaches. Our results show that accounting for how organisms interact with microclimates can fundamentally alter forecasts of species’ responses to climate change, highlighting the importance of integrating behavioural and mechanistic processes into predictive models.

## Repository content

### R Codes
**Microclimate and biophysical model**
- `RCode_microclimate_model.R`: generates microclimatic conditions at the study area using NicheMapR (Kearney et al. 2020) and microclima (Maclean et al. 2019)
- `RCode_biophysical_model.R`: validates the biophysical model and generates mechanistic layers (thermoregulatory inaccuracy and thermoregulatory window)

**Species distribution models (correlative vs hybrid)**
- `RCode_SDM_helpers.R`: shared library — package management, project constants, plotting palettes, and statistical / raster utilities used by all SDM scripts (`select07`, `evalSDM`, `crossvalSDM`, MESS with driver fallback, niche-overlap metrics, future-mechanistic-layer loader)
- `RCode_SDM_correlative.R`: Steps 1–10 of the SDM pipeline — climate data acquisition (WorldClim 2.1 + CMIP6), GBIF spatial thinning, `M_core`/`M_bg` accessible-area buffers, bias-corrected background sampling, `select07`+VIF variable selection, BRT model fitting under repeated 5-fold spatial cross-validation (10 repeats × 5 folds = 50 fits per model), threshold selection (P10 and maxSSS), and Wilcoxon signed-rank test on paired AUC and TSS
- `RCode_SDM_integrated.R`: robustness analyses for the integrated (correlative + behavioural mechanistic) SDM — Section A: prediction uncertainty (CV maps, MESS analysis); Section B: null model test with three null-hybrid configurations (topographic, quadratic-climate, Gaussian random fields); Section C: controlled comparison via a matched-correlative model that uses only the bioclimatic variables present in the hybrid
- `RCode_SDM_projections_future.R`: project the trained ensembles onto current and future (CMIP6 SSP2-4.5, MPI-ESM1-2-HR, 2041–2060) climate, and compute MESS extrapolation, niche overlap (Schoener's D, Pearson r) and range-change metrics
- `RCode_SDM_make_figures.R`: publication-ready figures (Paul Tol CVD-safe palette, MESS overlay, soft agreement maps, latitudinal profiles, niche-divergence panels), supplementary figures (S2.2 PDP response curves; S2.3 model-performance boxplots) and raster export bundle (GeoTIFF + ASCII with manifest, README and per-class area statistics)

### Microclimate data
- **Microclimate model output**
    - `ElPardo_X_metout.csv`: microclimate model output, year X (1997-2022) — sun-exposed conditions
    - `ElPardo_X_shadmet.csv`: microclimate model output, year X (1997-2022) — shaded conditions (75% shade)
    - `ElPardo_X_soil.csv`: microclimate model output, year X (1997-2022) — soil temperature (sun)
    - `ElPardo_X_shadsoil.csv`: microclimate model output, year X (1997-2022) — soil temperature (shade)
    - `field_Te.csv`: empirical operative temperature data (Diaz et al. 2022 Funct Ecol)
    - `field_Tb.csv`: empirical field body temperature data (Diaz et al. 2022 Funct Ecol)

### SDM data
- `Psammodromus_GBIF_cleaned.RData`: cleaned GBIF occurrences for *Psammodromus algirus* (38,279 raw → cleaned subset after geographic clipping, deduplication and quality filters)
- `Psammodromus_GBIF_DOI.txt`: DOI and filters of the GBIF download (DOI: 10.15468/dl.b23yf5)

### Sources
**Simulation results**
  - `dataMAY.R`, `dataJUNE.R`: body temperatures, operative temperatures, inaccuracy and thermoregulation window (May / June, current)
  - `dataMAY_warm.R`, `dataJUNE_warm.R`: same for the warm scenario

**Raster layers**
  - `xy.data.RData`: coordinates for raster layers
  - `cells.RData`: cell IDs
  - `map.grd/.gri`: base raster
  - `buffer_map.grd/.gri`: predicted buffer (delta Tb / delta Te)
  - `meandb_map.grd/.gri`: thermoregulatory inaccuracy (current)
  - `meanActivity_map.grd/.gri`: thermoregulatory window (current)
  - `meandb_warm_map.grd/.gri`: thermoregulatory inaccuracy (future warm scenario)
  - `meanActivity_warm_map.grd/.gri`: thermoregulatory window (future warm scenario)

### Supplementary R Codes
GCM/SSP sensitivity, decision and verification scripts that support the SDM pipeline:
- `13_sensitivity_GCM.R`: sensitivity of future projections across five alternative GCMs
- `15_range_change_uncertainty.R`: multi-GCM range-change uncertainty analysis
- `17_GCM_SSP_decision.R`: diagnostic analyses leading to the choice of MPI-ESM1-2-HR / SSP2-4.5 for the main paper
- `17b_SI_multi_GCM_table.R`: SI table summarising multi-GCM performance under SSP585
- `18_verify_MPI_extension.R`: verification of the MPI-ESM1-2-HR extension at the study area
- `18b_verify_MPI_report.R`: companion report for the MPI verification
- `19_GBIF_download_submit.R`, `20_GBIF_download_fetch.R`: reproducible GBIF download (submit + fetch)

### ODMAP
`ODMAP_Rubalcaba_etal.csv`: ODMAP protocol for species distribution modelling (Zurell et al. 2020).

## How to run the SDM analyses
1. Install the R packages declared in `RCode_SDM_helpers.R` (or restore the renv lockfile when available: `renv::restore()`).
2. The first run downloads ~6 GB of climate data via `geodata::worldclim_global()` and `geodata::cmip6_world()` into `data/`. This directory is gitignored.
3. Run the scripts in this order from the repo root (`setwd()` to the cloned repository):
   ```r
   source("RCode_SDM_correlative.R")
   source("RCode_SDM_projections_future.R")
   source("RCode_SDM_integrated.R")
   source("RCode_SDM_make_figures.R")
   ```
4. The default output directory is `output/` (override via `OUTDIR <- "..."` before sourcing the helpers).

## Code & Data Availability
**Code availability.** All custom code used in this study, including the species distribution models, biophysical model, and microclimate simulations, is available at https://github.com/JRubalcaba/behavioral-buffering-sdm. A permanent snapshot will be archived on Zenodo upon acceptance, and the DOI added to the final repository release.

**Data availability.** Species occurrence records were downloaded from GBIF.org (DOI: 10.15468/dl.b23yf5). Climate data are from WorldClim v2.1 (Fick & Hijmans 2017) at 2.5 arc-min resolution, and future projections from the CMIP6 MPI-ESM1-2-HR model (SSP2-4.5, 2041–2060) accessed through the `geodata` R package. Microclimate simulations were generated with NicheMapR (Kearney et al. 2020) and microclima (Maclean et al. 2019). Mechanistic layers (thermoregulatory inaccuracy, thermoregulation window, behavioural buffer) generated for this study are provided under `Sources/`. GCM/SSP-specific future mechanistic layers and intermediate model outputs (fitted BRT models, prediction rasters, performance metrics) will be deposited on Zenodo upon acceptance with a citable DOI.

# Session Info
**Microclimate and biophysical model**
R version 4.4.1 (2024-06-14 ucrt). Platform: x86_64-w64-mingw32/x64. Running under: Windows 11 x64 (build 26200).

Attached packages:
 [1] ggplot2_3.5.2      RColorBrewer_1.1-3 RNetCDF_2.9-2      RNCEP_1.0.10       maps_3.4.2.1       terra_1.8-60
 [7] lubridate_1.9.3    stringr_1.5.1      dplyr_1.2.0        tidyr_1.3.1        raster_3.6-32      sp_2.2-0
[13] microclima_0.1.0   NicheMapR_3.3.2

**SDM pipeline**
Exact package versions are recorded in `renv.lock` (added with the SDM scripts). Per-script `sessionInfo()` is also written to `output/session_info_<script>.txt` after each run.
