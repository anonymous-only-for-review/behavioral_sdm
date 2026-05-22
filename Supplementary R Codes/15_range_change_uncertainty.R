################################################################################
## 15_range_change_uncertainty_v2.R — Propagate model uncertainty to range change
##
## V2 NOTE: Uses MIROC6/SSP585-specific mechanistic layers
## (meandb_MIROC6_ssp585, meanActivity_MIROC6_ssp585) instead of the constant
## "_warm_map" layers used in v1. Outputs go to range_change_uncertainty_v2/.
## Scope is unchanged (MIROC6 / SSP585 only).
##
## Purpose: Compute per-model range change from all 50 BRT fits (10 repeats x
## 5 folds) for both correlative and hybrid models. Report median, 95% CI,
## and IQR of range change estimates rather than point estimates from the
## ensemble mean.
##
## Also computes per-model binary maps for both thresholds (maxSSS, P10).
################################################################################

source("RCode_SDM_helpers.R")

# Optional smoke-test mode: set V2_SMOKE=1 to run a reduced pipeline
SMOKE <- isTRUE(Sys.getenv("V2_SMOKE") == "1")

# Output directory
UNC_DIR <- file.path(OUTDIR, if (SMOKE) "range_change_uncertainty_v2_smoke" else "range_change_uncertainty_v2")
dirs_unc <- c("", "figures")
invisible(sapply(file.path(UNC_DIR, dirs_unc), dir.create,
                 showWarnings = FALSE, recursive = TRUE))

# ============================================================================ #
# REBUILD RASTERS FROM SOURCE FILES ----
# ============================================================================ #

# NOTE: terra SpatRaster objects lose C++ pointers when loaded from RData.
# We rebuild rasters from source files, but load data frames and gbm models
# (which are pure R objects) from RData checkpoints.

cat("\n=== Rebuilding raster layers from source files ===\n")

# Mechanistic layers (reference grid)
meandb2 <- terra::rast("data/mechanistic layers/mechanistic_layers_v4/meandb2_map.grd")
names(meandb2) <- "deviation_mean"
meanActivity2 <- terra::rast("data/mechanistic layers/mechanistic_layers_v4/meanActivity2_map.grd")
names(meanActivity2) <- "Activity"

# Bioclimatic layers
bio_curr <- geodata::worldclim_global(var = "bio", res = 2.5, download = FALSE, path = "data")
names(bio_curr) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
bio_curr <- terra::resample(bio_curr, meandb2) %>% terra::mask(meandb2)

# Land mask
world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_v  <- terra::vect(world_sf) %>% terra::project(terra::crs(meandb2)) %>%
  terra::crop(terra::ext(meandb2))
land_mask <- terra::rasterize(world_v, meandb2[[1]], field = 1)

# Selected variable names
vars_correlativo <- readRDS(file.path(OUTDIR, "models/selected_variables_correlative.rds"))
vars_hibrido     <- readRDS(file.path(OUTDIR, "models/selected_variables_hybrid.rds"))
vars_hybrido_climate <- setdiff(vars_hibrido, c("deviation_mean", "Activity"))

# Current predictor stacks
pred_correlativo <- terra::mask(bio_curr[[vars_correlativo]], land_mask, maskvalues = 0)
pred_hibrido     <- c(bio_curr[[vars_hybrido_climate]], meandb2, meanActivity2)
pred_hibrido     <- terra::mask(pred_hibrido, land_mask, maskvalues = 0)

# Future climate
bio_fut <- geodata::cmip6_world(model = "MIROC6", ssp = "585", time = "2041-2060",
                                var = "bioc", res = 2.5, path = "data", download = FALSE)
names(bio_fut) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
bio_fut <- terra::resample(bio_fut, meandb2) %>% terra::mask(meandb2)

# V2: MIROC6 / SSP585-specific mechanistic layers
mech_fut <- load_mech_future(gcm = "MIROC6", ssp = "ssp585")
meandb2_fut       <- mech_fut$meandb
meanActivity2_fut <- mech_fut$meanActivity

pred_correlativo_fut <- bio_fut[[vars_correlativo]]
pred_hibrido_fut     <- c(bio_fut[[vars_hybrido_climate]], meandb2_fut, meanActivity2_fut)
pred_correlativo_fut <- terra::mask(pred_correlativo_fut, land_mask, maskvalues = 0)
pred_hibrido_fut     <- terra::mask(pred_hibrido_fut,     land_mask, maskvalues = 0)

# Load gbm models and data frames from checkpoints (pure R objects survive load)
cat("\n=== Loading model checkpoints ===\n")
load(file.path(OUTDIR, "models/03_model_fitting.RData"))
load(file.path(OUTDIR, "models/04_evaluation.RData"))

if (SMOKE) {
  n_smoke <- 3
  super_models_corr <- super_models_corr[seq_len(n_smoke)]
  super_models_hyb  <- super_models_hyb[seq_len(n_smoke)]
  cat(sprintf("[V2_SMOKE=1] Reduced to %d fits per model type\n", n_smoke))
}

# Thresholds
thresholds <- readRDS(file.path(OUTDIR, "results_summary/thresholds_maxSSS_ensemble.rds"))
thr_corr_maxSSS <- thresholds$correlativo_mean_ensemble
thr_hyb_maxSSS  <- thresholds$hibrido_mean_ensemble
thr_corr_p10    <- thr_results_corr_p10$mean_threshold
thr_hyb_p10     <- thr_results_hyb_p10$mean_threshold

cat(sprintf("  Thresholds — Correlative: maxSSS=%.4f, P10=%.4f\n", thr_corr_maxSSS, thr_corr_p10))
cat(sprintf("  Thresholds — Hybrid:      maxSSS=%.4f, P10=%.4f\n", thr_hyb_maxSSS, thr_hyb_p10))

# ============================================================================ #
# CELL AREA ----
# ============================================================================ #

# Approximate cell area using latitude-corrected grid
cell_area_km2 <- prod(res(pred_correlativo)) * 111.32^2 *
  cos(mean(ext(pred_correlativo)[3:4]) * pi / 180)
cat(sprintf("  Approx cell area: %.1f km2\n", cell_area_km2))

# ============================================================================ #
# PER-MODEL PREDICTIONS AND RANGE CHANGE ----
# ============================================================================ #

compute_per_model_range_change <- function(models, pred_curr, pred_fut,
                                           thr_maxSSS, thr_p10, model_label) {
  n_models <- length(models)
  results <- vector("list", n_models)

  cat(sprintf("\n=== Per-model predictions: %s (%d models) ===\n", model_label, n_models))

  for (i in seq_len(n_models)) {
    if (i %% 10 == 0) cat(sprintf("  Model %d/%d\n", i, n_models))
    m <- models[[i]]
    model_vars <- m$var.names

    # Predict current
    r_curr <- terra::predict(pred_curr[[model_vars]], m, n.trees = m$best_ntree, type = "response")

    # Predict future
    r_fut  <- terra::predict(pred_fut[[model_vars]], m, n.trees = m$best_ntree, type = "response")

    # Range area under each threshold
    for (thr_name in c("maxSSS", "P10")) {
      thr_val <- if (thr_name == "maxSSS") thr_maxSSS else thr_p10

      n_curr <- sum(values(r_curr, na.rm = TRUE) >= thr_val)
      n_fut  <- sum(values(r_fut,  na.rm = TRUE) >= thr_val)
      area_curr <- n_curr * cell_area_km2
      area_fut  <- n_fut  * cell_area_km2
      pct_change <- if (area_curr > 0) (area_fut - area_curr) / area_curr * 100 else NA_real_

      # Jaccard: overlap / union
      bin_c <- values(r_curr, na.rm = TRUE) >= thr_val
      bin_f <- values(r_fut,  na.rm = TRUE) >= thr_val
      ok <- complete.cases(bin_c, bin_f)
      intersection <- sum(bin_c[ok] & bin_f[ok])
      union_n      <- sum(bin_c[ok] | bin_f[ok])
      jaccard      <- if (union_n > 0) intersection / union_n else NA_real_

      results[[length(results) + 1]] <- data.frame(
        Model = model_label, Model_idx = i, Threshold = thr_name,
        Area_current_km2 = area_curr, Area_future_km2 = area_fut,
        Change_pct = pct_change, Jaccard = jaccard
      )
    }
  }
  dplyr::bind_rows(results)
}

# Run for both model types
rc_corr <- compute_per_model_range_change(
  super_models_corr, pred_correlativo, pred_correlativo_fut,
  thr_corr_maxSSS, thr_corr_p10, "Correlative"
)

rc_hyb <- compute_per_model_range_change(
  super_models_hyb, pred_hibrido, pred_hibrido_fut,
  thr_hyb_maxSSS, thr_hyb_p10, "Hybrid"
)

rc_all <- rbind(rc_corr, rc_hyb)
readr::write_csv(rc_all, file.path(UNC_DIR, "per_model_range_change.csv"))

# ============================================================================ #
# SUMMARY STATISTICS ----
# ============================================================================ #

cat("\n=== Summary statistics ===\n")

rc_summary <- rc_all %>%
  dplyr::group_by(Model, Threshold) %>%
  dplyr::summarise(
    n_models = dplyr::n(),
    Area_curr_median = median(Area_current_km2, na.rm = TRUE),
    Area_curr_q025   = quantile(Area_current_km2, 0.025, na.rm = TRUE),
    Area_curr_q975   = quantile(Area_current_km2, 0.975, na.rm = TRUE),
    Area_fut_median  = median(Area_future_km2, na.rm = TRUE),
    Area_fut_q025    = quantile(Area_future_km2, 0.025, na.rm = TRUE),
    Area_fut_q975    = quantile(Area_future_km2, 0.975, na.rm = TRUE),
    Change_pct_median = median(Change_pct, na.rm = TRUE),
    Change_pct_q025   = quantile(Change_pct, 0.025, na.rm = TRUE),
    Change_pct_q975   = quantile(Change_pct, 0.975, na.rm = TRUE),
    Change_pct_IQR    = IQR(Change_pct, na.rm = TRUE),
    Change_pct_mean   = mean(Change_pct, na.rm = TRUE),
    Change_pct_sd     = sd(Change_pct, na.rm = TRUE),
    Jaccard_median    = median(Jaccard, na.rm = TRUE),
    Jaccard_q025      = quantile(Jaccard, 0.025, na.rm = TRUE),
    Jaccard_q975      = quantile(Jaccard, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(rc_summary, file.path(UNC_DIR, "range_change_summary.csv"))

cat("\nRange change summary (median [95% CI]):\n")
for (i in seq_len(nrow(rc_summary))) {
  row <- rc_summary[i, ]
  cat(sprintf("  %s (%s): %.1f%% [%.1f%%, %.1f%%]  |  Jaccard: %.3f [%.3f, %.3f]\n",
              row$Model, row$Threshold,
              row$Change_pct_median, row$Change_pct_q025, row$Change_pct_q975,
              row$Jaccard_median, row$Jaccard_q025, row$Jaccard_q975))
}

# ============================================================================ #
# ENSEMBLE AGREEMENT MAP ----
# ============================================================================ #

cat("\n=== Computing ensemble agreement maps ===\n")

# For each model, predict binary (maxSSS). Sum across 50 models = proportion of models
# predicting presence. This gives a continuous "agreement" surface.

compute_agreement <- function(models, pred_stack, thr, label, period) {
  model_vars <- models[[1]]$var.names
  pred_stack <- pred_stack[[model_vars]]
  n <- length(models)
  r_sum <- NULL
  for (i in seq_len(n)) {
    if (i %% 10 == 0) cat(sprintf("  [%s %s] Agreement model %d/%d\n", label, period, i, n))
    m <- models[[i]]
    r_pred <- terra::predict(pred_stack, m, n.trees = m$best_ntree, type = "response")
    r_bin  <- ifel(r_pred >= thr, 1, 0)
    if (is.null(r_sum)) {
      r_sum <- r_bin
    } else {
      r_sum <- r_sum + r_bin
    }
  }
  r_sum / n  # proportion of models predicting presence
}

agree_corr_curr <- compute_agreement(super_models_corr, pred_correlativo, thr_corr_maxSSS, "Corr", "current")
agree_corr_fut  <- compute_agreement(super_models_corr, pred_correlativo_fut, thr_corr_maxSSS, "Corr", "future")
agree_hyb_curr  <- compute_agreement(super_models_hyb, pred_hibrido, thr_hyb_maxSSS, "Hyb", "current")
agree_hyb_fut   <- compute_agreement(super_models_hyb, pred_hibrido_fut, thr_hyb_maxSSS, "Hyb", "future")

agree_corr_curr <- mask(agree_corr_curr, land_mask)
agree_corr_fut  <- mask(agree_corr_fut,  land_mask)
agree_hyb_curr  <- mask(agree_hyb_curr,  land_mask)
agree_hyb_fut   <- mask(agree_hyb_fut,   land_mask)

writeRaster(agree_corr_curr, file.path(UNC_DIR, "agreement_corr_current.tif"), overwrite = TRUE, datatype = "FLT4S")
writeRaster(agree_corr_fut,  file.path(UNC_DIR, "agreement_corr_future.tif"),  overwrite = TRUE, datatype = "FLT4S")
writeRaster(agree_hyb_curr,  file.path(UNC_DIR, "agreement_hyb_current.tif"),  overwrite = TRUE, datatype = "FLT4S")
writeRaster(agree_hyb_fut,   file.path(UNC_DIR, "agreement_hyb_future.tif"),   overwrite = TRUE, datatype = "FLT4S")

# ============================================================================ #
# FIGURES ----
# ============================================================================ #

cat("\n=== Generating figures ===\n")

# Figure 1: Boxplot of per-model range change
p_box <- ggplot(rc_all, aes(x = interaction(Model, Threshold), y = Change_pct, fill = Model)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Correlative" = "#E69F00", "Hybrid" = "#009E73")) +
  labs(x = NULL, y = "Range change (%)",
       title = "Per-model range change estimates (50 fits per model type)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave(file.path(UNC_DIR, "figures/range_change_boxplot.png"), p_box,
       width = 8, height = 5, dpi = 300)

# Figure 2: Density plot
p_dens <- ggplot(rc_all %>% dplyr::filter(Threshold == "maxSSS"),
                 aes(x = Change_pct, fill = Model, color = Model)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Correlative" = "#E69F00", "Hybrid" = "#009E73")) +
  scale_color_manual(values = c("Correlative" = "#E69F00", "Hybrid" = "#009E73")) +
  labs(x = "Range change (%)", y = "Density",
       title = "Distribution of range change estimates (maxSSS threshold)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(UNC_DIR, "figures/range_change_density.png"), p_dens,
       width = 7, height = 4.5, dpi = 300)

# ============================================================================ #
# SAVE ----
# ============================================================================ #

save(rc_all, rc_summary,
     agree_corr_curr, agree_corr_fut, agree_hyb_curr, agree_hyb_fut,
     file = file.path(UNC_DIR, "range_change_uncertainty_results.RData"))

cat("\n  OK: 15_range_change_uncertainty.R complete\n")
cat(sprintf("  Outputs saved to: %s\n", UNC_DIR))

sessionInfo()
