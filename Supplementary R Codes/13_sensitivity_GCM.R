################################################################################
## 13_sensitivity_gcm_v2.R — Multi-GCM sensitivity (SSP5-8.5, 2041-2060)
## Psammodromus algirus SDM: Correlative vs Hybrid
##
## V2 NOTE: Unlike v1, each GCM is projected with its OWN scenario-specific
## mechanistic layers (meandb_<GCM>_ssp585, meanActivity_<GCM>_ssp585) instead
## of the constant "_warm_map" layers. Outputs go to sensitivity_gcm_v2/ to
## allow side-by-side comparison with the v1 analysis.
##
## Projects trained BRT ensembles onto 5 GCMs under SSP5-8.5 to assess
## inter-GCM variability. All 5 GCMs are projected using the same ensemble
## models to ensure consistent baseline areas.
##
## GCMs selected following Cos et al. (2022, Earth System Dynamics) and
## Brands et al. (2013) for Iberian Peninsula performance:
##   - MIROC6 (reference GCM in main pipeline)
##   - MPI-ESM1-2-HR
##   - CNRM-CM6-1
##   - EC-Earth3-Veg
##   - UKESM1-0-LL
##
## Outputs:
##   - Per-GCM: continuous/binary rasters, range change CSVs
##   - Summary: multi-GCM comparison table and figures
##
## Requires: 00_config.R, trained ensembles, current climate rasters
################################################################################

# ============================================================================ #
# 0. CONFIGURATION ----
# ============================================================================ #

source("RCode_SDM_helpers.R")

# Optional smoke-test mode: set V2_SMOKE=1 to run a reduced pipeline
SMOKE <- isTRUE(Sys.getenv("V2_SMOKE") == "1")

OUTDIR_SENS <- file.path(OUTDIR, if (SMOKE) "sensitivity_gcm_v2_smoke" else "sensitivity_gcm_v2")
dir.create(file.path(OUTDIR_SENS, "figures"), showWarnings = FALSE, recursive = TRUE)

SSP_CODE  <- "585"
SSP_LABEL <- "ssp585"
TIME_PERIOD <- "2041-2060"

# All 5 GCMs to iterate (MIROC6 included for consistent baseline)
GCMS <- c(
  "MIROC6",
  "MPI-ESM1-2-HR",
  "CNRM-CM6-1",
  "EC-Earth3-Veg",
  "UKESM1-0-LL"
)
if (SMOKE) {
  GCMS <- GCMS[1]
  cat("[V2_SMOKE=1] Reduced to single GCM:", GCMS, "\n")
}

cat("\n========================================================================\n")
cat(" MULTI-GCM SENSITIVITY ANALYSIS (SSP5-8.5, 2041-2060)\n")
cat("========================================================================\n\n")

# ============================================================================ #
# 1. LOAD TRAINED MODELS AND EXISTING DATA ----
# ============================================================================ #

cat("=== STEP 1: Loading trained models & reference data ===\n")

# Ensemble models
mod_corr <- readRDS(file.path(OUTDIR, "models/ensemble_correlative_50.rds"))
mod_hyb  <- readRDS(file.path(OUTDIR, "models/ensemble_hybrid_50.rds"))
cat(sprintf("  Correlative ensemble: %d models, vars = %s\n",
            length(mod_corr$models), paste(mod_corr$predictors, collapse = ", ")))
cat(sprintf("  Hybrid ensemble:      %d models, vars = %s\n",
            length(mod_hyb$models),  paste(mod_hyb$predictors, collapse = ", ")))

# Selected variables — use ensemble predictors (may differ from variable selection RDS
# if ensemble was retrained with a different subset)
vars_correlativo <- mod_corr$predictors
vars_hibrido     <- mod_hyb$predictors

# Land mask
land_mask <- rast(file.path(OUTDIR, "rasters/land_mask.tif"))

# Thresholds
# Thresholds (from 04_evaluation.RData — single source of truth)
load(file.path(OUTDIR, "models/04_evaluation.RData"))
thr_corr_p10 <- thr_results_corr_p10$mean_threshold
thr_hyb_p10  <- thr_results_hyb_p10$mean_threshold

cat(sprintf("  Thresholds P10:    corr=%.4f | hyb=%.4f\n", thr_corr_p10, thr_hyb_p10))
cat(sprintf("  Thresholds maxSSS: corr=%.4f | hyb=%.4f\n", thr_corr_maxSSS, thr_hyb_maxSSS))

# Reference raster for alignment
meandb2_ref <- rast("data/mechanistic layers/mechanistic_layers_v4/meandb2_map.grd")

# Mechanistic layers — current (baseline, unchanged between v1 and v2)
meandb2_curr       <- rast("data/mechanistic layers/mechanistic_layers_v4/meandb2_map.grd")
names(meandb2_curr) <- "deviation_mean"
meanActivity2_curr <- rast("data/mechanistic layers/mechanistic_layers_v4/meanActivity2_map.grd")
names(meanActivity2_curr) <- "Activity"

# V2: future mechanistic layers are GCM/SSP-specific and loaded inside the loop
# (see load_mech_future() in helpers_mech_v2.R)

# --- Build current predictions using ensemble models (consistent baseline) ---
cat("  Building current predictions with ensemble models...\n")

bio_curr_path <- list.files("data/climate/wc2.1_2.5m", pattern = "wc2.1_2.5m_bio_.*\\.tif$",
                             full.names = TRUE)
bio_curr_stack <- rast(bio_curr_path)
layer_nums <- as.integer(gsub(".*bio_(\\d+).*", "\\1", basename(bio_curr_path)))
names(bio_curr_stack) <- ifelse(layer_nums < 10,
                                 paste0("bio0", layer_nums),
                                 paste0("bio", layer_nums))
bio_curr_stack <- resample(bio_curr_stack, meandb2_ref) %>% mask(meandb2_ref)

# Correlative current stack
pred_corr_curr_stack <- bio_curr_stack[[vars_correlativo]]
pred_corr_curr_stack <- mask(pred_corr_curr_stack, land_mask, maskvalues = 0)

# Hybrid current stack
vars_hybrido_climate <- setdiff(vars_hibrido, c("deviation_mean", "Activity"))
pred_hyb_curr_stack <- c(bio_curr_stack[[vars_hybrido_climate]], meandb2_curr, meanActivity2_curr)
pred_hyb_curr_stack <- mask(pred_hyb_curr_stack, land_mask, maskvalues = 0)

cat("  All reference data loaded\n")

# ============================================================================ #
# 2. HELPER FUNCTIONS ----
# ============================================================================ #

predict_ensemble <- function(models, pred_stack, return_sd = TRUE) {
  model_vars <- models[[1]]$var.names
  stack_vars <- names(pred_stack)
  missing_vars <- setdiff(model_vars, stack_vars)
  if (length(missing_vars) > 0)
    stop(sprintf("Missing vars in stack: %s", paste(missing_vars, collapse = ", ")))
  pred_stack <- pred_stack[[model_vars]]
  preds <- lapply(seq_along(models), function(i) {
    m <- models[[i]]
    terra::predict(pred_stack, m, n.trees = m$best_ntree, type = "response")
  })
  pred_stack_all <- rast(preds)
  r_mean <- app(pred_stack_all, fun = mean, na.rm = TRUE)
  names(r_mean) <- "suitability_mean"
  if (!return_sd) return(r_mean)
  r_sd <- app(pred_stack_all, fun = sd, na.rm = TRUE)
  names(r_sd) <- "suitability_sd"
  list(r_mean = r_mean, r_sd = r_sd, r_stack = pred_stack_all)
}

calculate_range_metrics <- function(current, future, model_name) {
  cell_area_km2 <- prod(res(current)) * 111^2
  stable <- current & future
  gained <- (!current) & future
  lost   <- current & (!future)
  area_current <- global(current, "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_future  <- global(future,  "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_stable  <- global(stable,  "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_gained  <- global(gained,  "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_lost    <- global(lost,    "sum", na.rm = TRUE)[[1]] * cell_area_km2
  jaccard <- area_stable / (area_current + area_future - area_stable)
  data.frame(
    Model          = model_name,
    Current_km2    = area_current,
    Future_km2     = area_future,
    Stable_km2     = area_stable,
    Gained_km2     = area_gained,
    Lost_km2       = area_lost,
    Net_change_km2 = area_future - area_current,
    Percent_change = 100 * (area_future - area_current) / area_current,
    Jaccard        = jaccard
  )
}

# ============================================================================ #
# 3. CURRENT PREDICTIONS & BINARY MAPS (consistent baseline) ----
# ============================================================================ #

cat("\n=== STEP 3: Current predictions with ensemble models ===\n")

pred_corr_curr <- predict_ensemble(mod_corr$models, pred_corr_curr_stack)
pred_corr_curr$r_mean <- mask(pred_corr_curr$r_mean, land_mask)

pred_hyb_curr <- predict_ensemble(mod_hyb$models, pred_hyb_curr_stack)
pred_hyb_curr$r_mean <- mask(pred_hyb_curr$r_mean, land_mask)

# Binarize current (same thresholds, consistent models)
bin_corr_curr_p10    <- ifel(is.na(pred_corr_curr$r_mean), NA,
                              pred_corr_curr$r_mean >= thr_corr_p10)
bin_hyb_curr_p10     <- ifel(is.na(pred_hyb_curr$r_mean), NA,
                              pred_hyb_curr$r_mean  >= thr_hyb_p10)
bin_corr_curr_maxSSS <- ifel(is.na(pred_corr_curr$r_mean), NA,
                              pred_corr_curr$r_mean >= thr_corr_maxSSS)
bin_hyb_curr_maxSSS  <- ifel(is.na(pred_hyb_curr$r_mean), NA,
                              pred_hyb_curr$r_mean  >= thr_hyb_maxSSS)

cell_area_km2 <- prod(res(bin_corr_curr_p10)) * 111^2
cat(sprintf("  Current range (P10):    corr=%.0f km2 | hyb=%.0f km2\n",
            global(bin_corr_curr_p10, "sum", na.rm = TRUE)[[1]] * cell_area_km2,
            global(bin_hyb_curr_p10,  "sum", na.rm = TRUE)[[1]] * cell_area_km2))
cat(sprintf("  Current range (maxSSS): corr=%.0f km2 | hyb=%.0f km2\n",
            global(bin_corr_curr_maxSSS, "sum", na.rm = TRUE)[[1]] * cell_area_km2,
            global(bin_hyb_curr_maxSSS,  "sum", na.rm = TRUE)[[1]] * cell_area_km2))

# ============================================================================ #
# 4. ITERATE OVER GCMs ----
# ============================================================================ #

# Storage for per-GCM results
all_metrics    <- list()
all_suit_corr  <- list()
all_suit_hyb   <- list()
gcms_processed <- character(0)

for (gcm in GCMS) {

  cat(sprintf("\n--- GCM: %s ---\n", gcm))

  # Create output dirs
  gcm_dir <- file.path(OUTDIR_SENS, gsub("[^A-Za-z0-9_-]", "_", gcm))
  dir.create(file.path(gcm_dir, "rasters"), showWarnings = FALSE, recursive = TRUE)

  # --- 3a. Download bioclimatic layers ---
  cat(sprintf("  [%s] Downloading bioclimatic layers...\n", gcm))
  bio_fut <- tryCatch(
    geodata::cmip6_world(
      model    = gcm,
      ssp      = SSP_CODE,
      time     = TIME_PERIOD,
      var      = "bioc",
      res      = 2.5,
      path     = "data",
      download = TRUE
    ),
    error = function(e) {
      warning(sprintf("Failed to download %s: %s", gcm, conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(bio_fut)) {
    cat(sprintf("  [%s] SKIPPED (download failed)\n", gcm))
    next
  }
  names(bio_fut) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))

  # --- 3b. Resample and mask ---
  cat(sprintf("  [%s] Resampling to study area...\n", gcm))
  bio_fut <- resample(bio_fut, meandb2_ref) %>% mask(meandb2_ref)

  # --- 3c. Load GCM/SSP-specific mechanistic layers (v2 change) ---
  cat(sprintf("  [%s] Loading %s-specific mechanistic layers...\n", gcm, SSP_LABEL))
  mech_fut <- load_mech_future(gcm = gcm, ssp = SSP_LABEL)
  meandb2_fut       <- mech_fut$meandb
  meanActivity2_fut <- mech_fut$meanActivity

  # --- 3c. Build predictor stacks ---
  cat(sprintf("  [%s] Building predictor stacks...\n", gcm))
  pred_corr_stack <- bio_fut[[vars_correlativo]]
  pred_hyb_stack  <- c(bio_fut[[vars_hybrido_climate]], meandb2_fut, meanActivity2_fut)

  pred_corr_stack <- mask(pred_corr_stack, land_mask, maskvalues = 0)
  pred_hyb_stack  <- mask(pred_hyb_stack,  land_mask, maskvalues = 0)

  # --- 3d. Project ensembles ---
  cat(sprintf("  [%s] Projecting correlative ensemble...\n", gcm))
  pred_corr <- predict_ensemble(mod_corr$models, pred_corr_stack)
  pred_corr$r_mean <- mask(pred_corr$r_mean, land_mask)
  pred_corr$r_sd   <- mask(pred_corr$r_sd,   land_mask)

  cat(sprintf("  [%s] Projecting hybrid ensemble...\n", gcm))
  pred_hyb <- predict_ensemble(mod_hyb$models, pred_hyb_stack)
  pred_hyb$r_mean <- mask(pred_hyb$r_mean, land_mask)
  pred_hyb$r_sd   <- mask(pred_hyb$r_sd,   land_mask)

  # Save continuous rasters
  gcm_tag <- gsub("[^A-Za-z0-9_-]", "_", gcm)
  writeRaster(pred_corr$r_mean,
              file.path(gcm_dir, "rasters", paste0("future_correlative_mean_", gcm_tag, ".tif")),
              overwrite = TRUE, datatype = "FLT4S")
  writeRaster(pred_hyb$r_mean,
              file.path(gcm_dir, "rasters", paste0("future_hybrid_mean_", gcm_tag, ".tif")),
              overwrite = TRUE, datatype = "FLT4S")

  # Store for summary panel
  all_suit_corr[[gcm]] <- pred_corr$r_mean
  all_suit_hyb[[gcm]]  <- pred_hyb$r_mean

  # --- 3e. Binarize ---
  cat(sprintf("  [%s] Binarizing...\n", gcm))

  bin_corr_p10    <- ifel(is.na(pred_corr$r_mean), NA, pred_corr$r_mean >= thr_corr_p10)
  bin_hyb_p10     <- ifel(is.na(pred_hyb$r_mean),  NA, pred_hyb$r_mean  >= thr_hyb_p10)
  bin_corr_maxSSS <- ifel(is.na(pred_corr$r_mean), NA, pred_corr$r_mean >= thr_corr_maxSSS)
  bin_hyb_maxSSS  <- ifel(is.na(pred_hyb$r_mean),  NA, pred_hyb$r_mean  >= thr_hyb_maxSSS)

  writeRaster(bin_corr_p10,
              file.path(gcm_dir, "rasters", paste0("binary_correlative_P10_", gcm_tag, ".tif")),
              overwrite = TRUE, datatype = "INT1U")
  writeRaster(bin_hyb_p10,
              file.path(gcm_dir, "rasters", paste0("binary_hybrid_P10_", gcm_tag, ".tif")),
              overwrite = TRUE, datatype = "INT1U")
  writeRaster(bin_corr_maxSSS,
              file.path(gcm_dir, "rasters", paste0("binary_correlative_maxSSS_", gcm_tag, ".tif")),
              overwrite = TRUE, datatype = "INT1U")
  writeRaster(bin_hyb_maxSSS,
              file.path(gcm_dir, "rasters", paste0("binary_hybrid_maxSSS_", gcm_tag, ".tif")),
              overwrite = TRUE, datatype = "INT1U")

  # --- 3f. Range change metrics ---
  cat(sprintf("  [%s] Computing range change metrics...\n", gcm))

  m_corr_p10    <- calculate_range_metrics(bin_corr_curr_p10,    bin_corr_p10,    "Correlative")
  m_hyb_p10     <- calculate_range_metrics(bin_hyb_curr_p10,     bin_hyb_p10,     "Hybrid")
  m_corr_maxSSS <- calculate_range_metrics(bin_corr_curr_maxSSS, bin_corr_maxSSS, "Correlative")
  m_hyb_maxSSS  <- calculate_range_metrics(bin_hyb_curr_maxSSS,  bin_hyb_maxSSS,  "Hybrid")

  metrics_p10    <- bind_rows(m_corr_p10, m_hyb_p10)    %>% mutate(GCM = gcm, Threshold = "P10")
  metrics_maxSSS <- bind_rows(m_corr_maxSSS, m_hyb_maxSSS) %>% mutate(GCM = gcm, Threshold = "maxSSS")

  gcm_metrics <- bind_rows(metrics_p10, metrics_maxSSS)
  readr::write_csv(gcm_metrics, file.path(gcm_dir, paste0("range_change_metrics_", gcm_tag, ".csv")))

  all_metrics[[gcm]] <- gcm_metrics
  gcms_processed <- c(gcms_processed, gcm)

  cat(sprintf("  [%s] Done\n", gcm))
}

cat(sprintf("\n  GCMs processed: %d / %d\n", length(gcms_processed), length(GCMS)))

# ============================================================================ #
# 5. COMPILE MULTI-GCM SUMMARY ----
# ============================================================================ #

cat("\n=== STEP 5: Compiling multi-GCM summary ===\n")

# All GCMs processed in the loop — no external CSV loading needed
summary_all <- bind_rows(all_metrics) %>%
  dplyr::select(GCM, Model, Threshold, Current_km2, Future_km2,
                Percent_change, Jaccard, everything())

readr::write_csv(summary_all,
                 file.path(OUTDIR_SENS, "summary_range_change_all_gcms.csv"))

cat("  Summary table saved\n")
cat("\nRange change summary (maxSSS threshold):\n")
summary_all %>%
  dplyr::filter(Threshold == "maxSSS") %>%
  dplyr::select(GCM, Model, Current_km2, Future_km2, Percent_change, Jaccard) %>%
  as.data.frame() %>% print(digits = 3)

# ============================================================================ #
# 6. FIGURE 1: Range change barplot across GCMs ----
# ============================================================================ #

cat("\n=== STEP 6: Multi-GCM range change barplot ===\n")

df_bar <- summary_all %>%
  dplyr::filter(Threshold == "maxSSS") %>%
  dplyr::select(GCM, Model, Percent_change)

# Compute mean and range across GCMs per model type
df_bar_summary <- df_bar %>%
  group_by(Model) %>%
  summarise(
    mean_change = mean(Percent_change),
    min_change  = min(Percent_change),
    max_change  = max(Percent_change),
    .groups = "drop"
  )

fig_bar <- ggplot(df_bar, aes(x = GCM, y = Percent_change, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
  scale_fill_manual(values = c("Correlative" = "#1F78B4", "Hybrid" = "#33A02C")) +
  labs(
    y     = "Range change (%)",
    x     = NULL,
    fill  = "Model type",
    title = "Range change across GCMs (SSP5-8.5, maxSSS threshold)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position     = "bottom",
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    plot.title          = element_text(face = "bold", hjust = 0.5, size = 10),
    axis.text.x         = element_text(angle = 30, hjust = 1, size = 8)
  )

ggsave(file.path(OUTDIR_SENS, "figures/range_change_across_gcms.pdf"),
       fig_bar, width = 180, height = 120, units = "mm", dpi = 300)
cat("  Barplot saved\n")

# ============================================================================ #
# 7. FIGURE 2: Suitability maps panel (5 GCMs x 2 model types) ----
# ============================================================================ #

cat("\n=== STEP 7: Multi-GCM suitability panel figure ===\n")

# Order: MIROC6 first, then remaining alphabetical
gcm_order <- c("MIROC6", sort(setdiff(gcms_processed, "MIROC6")))
gcm_order <- gcm_order[gcm_order %in% names(all_suit_corr)]

if (length(gcm_order) >= 2) {

  # Convert rasters to data frames for ggplot
  suit_dfs <- list()
  for (g in gcm_order) {
    df_c <- as.data.frame(all_suit_corr[[g]], xy = TRUE, na.rm = TRUE)
    names(df_c)[3] <- "suitability"
    df_c$GCM   <- g
    df_c$Model <- "Correlative"

    df_h <- as.data.frame(all_suit_hyb[[g]], xy = TRUE, na.rm = TRUE)
    names(df_h)[3] <- "suitability"
    df_h$GCM   <- g
    df_h$Model <- "Hybrid"

    suit_dfs[[g]] <- bind_rows(df_c, df_h)
  }

  df_panel <- bind_rows(suit_dfs) %>%
    mutate(
      GCM   = factor(GCM, levels = gcm_order),
      Model = factor(Model, levels = c("Correlative", "Hybrid"))
    )

  fig_panel <- ggplot(df_panel, aes(x = x, y = y, fill = suitability)) +
    geom_raster() +
    facet_grid(Model ~ GCM) +
    scale_fill_viridis_c(limits = c(0, 1), name = "Suitability", option = "viridis") +
    coord_equal(expand = FALSE) +
    labs(x = NULL, y = NULL,
         title = "Future suitability across GCMs (SSP5-8.5, 2041-2060)") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid       = element_blank(),
      strip.text       = element_text(size = 8, face = "bold"),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "bottom",
      legend.key.width = unit(20, "mm"),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 10)
    )

  n_gcms <- length(gcm_order)
  fig_width  <- max(180, n_gcms * 45)
  fig_height <- 120

  ggsave(file.path(OUTDIR_SENS, "figures/suitability_maps_all_gcms.pdf"),
         fig_panel, width = fig_width, height = fig_height, units = "mm", dpi = 300)
  cat("  Suitability panel saved\n")

} else {
  cat("  WARNING: Fewer than 2 GCMs available — panel figure skipped\n")
}

# ============================================================================ #
# DONE ----
# ============================================================================ #

cat("\n========================================================================\n")
cat(" MULTI-GCM SENSITIVITY ANALYSIS COMPLETE\n")
cat(sprintf(" GCMs processed: %s\n", paste(gcms_processed, collapse = ", ")))
cat(sprintf(" Outputs: %s\n", OUTDIR_SENS))
cat("========================================================================\n\n")

cat("Output summary:\n")
cat(sprintf("  Summary CSV:   %s\n",
            file.path(OUTDIR_SENS, "summary_range_change_all_gcms.csv")))
cat(sprintf("  Figures:       %d files\n",
            length(list.files(file.path(OUTDIR_SENS, "figures"), pattern = "\\.pdf$"))))
cat(sprintf("  Per-GCM dirs:  %d\n", length(gcms_processed)))

sessionInfo()
