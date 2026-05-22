################################################################################
## RCode_SDM_projections_future.R
##
## Project the trained correlative and hybrid SDM ensembles onto current and
## future (CMIP6 SSP2-4.5, MPI-ESM1-2-HR, 2041-2060) climate, and compute the
## downstream niche-overlap, range-change and MESS extrapolation metrics that
## drive Figs. 5-8 of:
##
##   Rubalcaba, Fandos & Diaz. "Behavioural thermoregulation and microclimate
##   reshape climate-driven range forecasts."
##
## Pipeline steps:
##   1. Load trained ensembles + selected variables + thresholds (from
##      RCode_SDM_correlative.R outputs)
##   2. Load current and future bioclim + mechanistic predictors
##   3. Project ensembles onto current and future climate, with per-cell SD
##   4. Binarize predictions with P10 and maxSSS thresholds
##   5. Compute MESS extrapolation rasters (point-of-reference = training data)
##   6. Compute niche overlap (Schoener's D, Pearson r, Spearman rho)
##   7. Compute range-change metrics (current vs future)
##   8. Build categorical agreement maps and continuous difference maps
##
## Notes:
##   - The chosen GCM/SSP for the main paper is MPI-ESM1-2-HR / SSP2-4.5,
##     selected via the diagnostic analyses in
##     "Supplementary R Codes/17_GCM_SSP_decision.R" and verified in
##     "Supplementary R Codes/18_verify_MPI_extension.R".
##   - GCM sensitivity (5 alternative GCMs under SSP585) and multi-GCM range-
##     change uncertainty are computed in
##     "Supplementary R Codes/13_sensitivity_GCM.R" and
##     "Supplementary R Codes/15_range_change_uncertainty.R".
##
## Inputs (relative to repo root):
##   - <OUTDIR>/models/{01_data_loading, 02_variable_selection,
##                      03_model_fitting, 04_evaluation}.RData
##   - <OUTDIR>/models/selected_variables_{correlative,hybrid}.rds
##   - <OUTDIR>/results_summary/thresholds_maxSSS_ensemble.rds
##   - Sources/meandb_map.grd, Sources/meanActivity_map.grd       (current)
##   - Sources/SDM/mechanistic_future/                            (future GCMs)
##   - data/cmip6/...                                             (auto-download)
##
## Outputs (under OUTDIR):
##   - rasters/{current,future}_{correlative,hybrid}_{mean,sd}.tif
##   - rasters/binary_{current,future}_{correlative,hybrid}_{p10,maxSSS}_PA.tif
##   - rasters/{categorical_present,categorical_future,land_mask}.tif
##   - rasters/{diff_continuous_present,diff_continuous_future}.tif
##   - rasters/mess_{correlative,hybrid}_future_pointsRef.tif
##   - niche_overlap_D_r.csv
##   - range_change_metrics{,_P10}.csv
##   - models/05_projections.RData, models/06_comparison.RData
##
## Author: Guillermo Fandos (gfandos@ucm.es)
################################################################################

source("RCode_SDM_helpers.R")
set.seed(MASTER_SEED)


# ============================================================================ #
# Configuration                                                                #
# ============================================================================ #

GCM_MODEL <- "MPI-ESM1-2-HR"
SSP_LABEL <- "ssp245"
SSP_GD    <- "245"   # geodata SSP code

# Paths to mechanistic layers (current = Juanvi's Sources/; future = SSP-specific)
MECH_CURR_MEANDB_PATH    <- file.path("Sources", "meandb_map.grd")
MECH_CURR_ACTIVITY_PATH  <- file.path("Sources", "meanActivity_map.grd")
# Future MPI/SSP245 mechanistic layers (loaded via load_mech_future from helpers)


# ============================================================================ #
# Step 1. Load trained ensembles, selected variables and thresholds            #
# ============================================================================ #

cat("\n=== STEP 1: Loading trained ensembles ===\n")

load(file.path(OUTDIR, "models/01_data_loading.RData"))
load(file.path(OUTDIR, "models/02_variable_selection.RData"))
load(file.path(OUTDIR, "models/03_model_fitting.RData"))
load(file.path(OUTDIR, "models/04_evaluation.RData"))

vars_correlativo <- readRDS(file.path(OUTDIR, "models/selected_variables_correlative.rds"))
vars_hibrido     <- readRDS(file.path(OUTDIR, "models/selected_variables_hybrid.rds"))

# Rebuild current rasters from disk (terra pointers do not survive .RData)
meandb2       <- terra::rast(MECH_CURR_MEANDB_PATH);   names(meandb2)       <- "deviation_mean"
meanActivity2 <- terra::rast(MECH_CURR_ACTIVITY_PATH); names(meanActivity2) <- "Activity"

bio_curr <- geodata::worldclim_global(var = "bio", res = 2.5,
                                      download = FALSE, path = "data")
names(bio_curr) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
bio_curr <- terra::resample(bio_curr, meandb2) |> terra::mask(meandb2)

world_sf  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_v   <- terra::vect(world_sf) |>
  terra::project(terra::crs(meandb2)) |>
  terra::crop(terra::ext(meandb2))
land_mask <- terra::rasterize(world_v, meandb2[[1]], field = 1)
terra::writeRaster(land_mask, file.path(OUTDIR, "rasters/land_mask.tif"),
                   overwrite = TRUE, datatype = "INT1U")

pred_correlativo <- terra::mask(bio_curr[[vars_correlativo]],     land_mask, maskvalues = 0)
pred_hibrido     <- c(bio_curr[[setdiff(vars_hibrido, c("deviation_mean", "Activity"))]],
                      meandb2, meanActivity2) |>
                    terra::mask(land_mask, maskvalues = 0)


# ============================================================================ #
# Step 2. Load future climate and mechanistic layers                           #
# ============================================================================ #

cat(sprintf("\n=== STEP 2: Loading %s / %s climate ===\n", GCM_MODEL, SSP_LABEL))

bio_fut <- geodata::cmip6_world(model = GCM_MODEL, ssp = SSP_GD,
                                time = "2041-2060", var = "bioc",
                                res = 2.5, path = "data/cmip6", download = TRUE)
names(bio_fut) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
bio_fut <- terra::resample(bio_fut, meandb2) |> terra::mask(meandb2)

mech_fut          <- load_mech_future(gcm = GCM_MODEL, ssp = SSP_LABEL)
meandb2_fut       <- mech_fut$meandb
meanActivity2_fut <- mech_fut$meanActivity

vars_hybrido_climate <- setdiff(vars_hibrido, c("deviation_mean", "Activity"))
pred_correlativo_fut <- bio_fut[[vars_correlativo]] |>
  terra::mask(land_mask, maskvalues = 0)
pred_hibrido_fut     <- c(bio_fut[[vars_hybrido_climate]], meandb2_fut, meanActivity2_fut) |>
  terra::mask(land_mask, maskvalues = 0)


# ============================================================================ #
# Step 3. Project ensembles onto current and future climate                    #
# ============================================================================ #

predict_ensemble <- function(models, pred_stack, return_sd = TRUE) {
  model_vars   <- models[[1]]$var.names
  missing_vars <- setdiff(model_vars, names(pred_stack))
  if (length(missing_vars) > 0)
    stop(sprintf("Missing vars: %s", paste(missing_vars, collapse = ", ")))
  pred_stack <- pred_stack[[model_vars]]
  preds <- lapply(seq_along(models), function(i) {
    m <- models[[i]]
    terra::predict(pred_stack, m, n.trees = m$best_ntree, type = "response")
  })
  pred_stack_all <- terra::rast(preds)
  r_mean <- terra::app(pred_stack_all, fun = mean, na.rm = TRUE)
  names(r_mean) <- "suitability_mean"
  if (!return_sd) return(r_mean)
  r_sd <- terra::app(pred_stack_all, fun = sd, na.rm = TRUE)
  names(r_sd) <- "suitability_sd"
  list(r_mean = r_mean, r_sd = r_sd, r_stack = pred_stack_all)
}

cat("\n=== STEP 3: Projecting ensembles ===\n")

pred_corr_curr <- predict_ensemble(super_models_corr, pred_correlativo)
pred_hyb_curr  <- predict_ensemble(super_models_hyb,  pred_hibrido)
pred_corr_fut  <- predict_ensemble(super_models_corr, pred_correlativo_fut)
pred_hyb_fut   <- predict_ensemble(super_models_hyb,  pred_hibrido_fut)

for (obj_name in c("pred_corr_curr", "pred_hyb_curr", "pred_corr_fut", "pred_hyb_fut")) {
  obj <- get(obj_name)
  obj$r_mean <- terra::mask(obj$r_mean, land_mask)
  obj$r_sd   <- terra::mask(obj$r_sd,   land_mask)
  assign(obj_name, obj)
}

terra::writeRaster(pred_corr_curr$r_mean,
                   file.path(OUTDIR, "rasters/current_correlative_mean.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(pred_corr_curr$r_sd,
                   file.path(OUTDIR, "rasters/current_correlative_sd.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(pred_hyb_curr$r_mean,
                   file.path(OUTDIR, "rasters/current_hybrid_mean.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(pred_hyb_curr$r_sd,
                   file.path(OUTDIR, "rasters/current_hybrid_sd.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(pred_corr_fut$r_mean,
                   file.path(OUTDIR, "rasters/future_correlative_mean.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(pred_corr_fut$r_sd,
                   file.path(OUTDIR, "rasters/future_correlative_sd.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(pred_hyb_fut$r_mean,
                   file.path(OUTDIR, "rasters/future_hybrid_mean.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(pred_hyb_fut$r_sd,
                   file.path(OUTDIR, "rasters/future_hybrid_sd.tif"),
                   overwrite = TRUE, datatype = "FLT4S")


# ============================================================================ #
# Step 4. Binarization (P10 and maxSSS)                                        #
# ============================================================================ #

cat("\n=== STEP 4: Binarising suitability rasters ===\n")

thr_corr_p10 <- thr_results_corr_p10$mean_threshold
thr_hyb_p10  <- thr_results_hyb_p10$mean_threshold

bin_corr_curr_p10    <- terra::ifel(is.na(pred_corr_curr$r_mean), NA, pred_corr_curr$r_mean >= thr_corr_p10)
bin_corr_fut_p10     <- terra::ifel(is.na(pred_corr_fut$r_mean),  NA, pred_corr_fut$r_mean  >= thr_corr_p10)
bin_hyb_curr_p10     <- terra::ifel(is.na(pred_hyb_curr$r_mean),  NA, pred_hyb_curr$r_mean  >= thr_hyb_p10)
bin_hyb_fut_p10      <- terra::ifel(is.na(pred_hyb_fut$r_mean),   NA, pred_hyb_fut$r_mean   >= thr_hyb_p10)

bin_corr_curr_maxSSS <- terra::ifel(is.na(pred_corr_curr$r_mean), NA, pred_corr_curr$r_mean >= thr_corr_maxSSS)
bin_corr_fut_maxSSS  <- terra::ifel(is.na(pred_corr_fut$r_mean),  NA, pred_corr_fut$r_mean  >= thr_corr_maxSSS)
bin_hyb_curr_maxSSS  <- terra::ifel(is.na(pred_hyb_curr$r_mean),  NA, pred_hyb_curr$r_mean  >= thr_hyb_maxSSS)
bin_hyb_fut_maxSSS   <- terra::ifel(is.na(pred_hyb_fut$r_mean),   NA, pred_hyb_fut$r_mean   >= thr_hyb_maxSSS)

terra::writeRaster(bin_corr_curr_p10,    file.path(OUTDIR, "rasters/binary_current_correlative_p10.tif"),     overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(bin_corr_fut_p10,     file.path(OUTDIR, "rasters/binary_future_correlative_p10.tif"),      overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(bin_hyb_curr_p10,     file.path(OUTDIR, "rasters/binary_current_hybrid_p10.tif"),          overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(bin_hyb_fut_p10,      file.path(OUTDIR, "rasters/binary_future_hybrid_p10.tif"),           overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(bin_corr_curr_maxSSS, file.path(OUTDIR, "rasters/binary_current_correlative_maxSSS_PA.tif"), overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(bin_corr_fut_maxSSS,  file.path(OUTDIR, "rasters/binary_future_correlative_maxSSS_PA.tif"),  overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(bin_hyb_curr_maxSSS,  file.path(OUTDIR, "rasters/binary_current_hybrid_maxSSS_PA.tif"),    overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(bin_hyb_fut_maxSSS,   file.path(OUTDIR, "rasters/binary_future_hybrid_maxSSS_PA.tif"),     overwrite = TRUE, datatype = "INT1U")


# ============================================================================ #
# Step 5. MESS extrapolation rasters (training-points reference)               #
# ============================================================================ #

cat("\n=== STEP 5: MESS extrapolation analysis ===\n")

# Reference data: training cells used in 03_model_fitting (already extracted)
ref_corr_df <- dat_corr[, vars_correlativo, drop = FALSE]
ref_hyb_df  <- dat_hyb[,  vars_hibrido,    drop = FALSE]

r_mess_corr <- mess_robusto(pred_correlativo_fut, ref_corr_df)
r_mess_hyb  <- mess_robusto(pred_hibrido_fut,    ref_hyb_df)
r_mess_corr <- terra::mask(r_mess_corr, land_mask)
r_mess_hyb  <- terra::mask(r_mess_hyb,  land_mask)

terra::writeRaster(r_mess_corr,
                   file.path(OUTDIR, "rasters/mess_correlative_future_pointsRef.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(r_mess_hyb,
                   file.path(OUTDIR, "rasters/mess_hybrid_future_pointsRef.tif"),
                   overwrite = TRUE, datatype = "FLT4S")


# ============================================================================ #
# Step 6. Categorical agreement maps                                           #
# ============================================================================ #

cat("\n=== STEP 6: Categorical agreement maps ===\n")

pres_cat <- cat_map_from(bin_corr_curr_maxSSS, bin_hyb_curr_maxSSS)
fut_cat  <- cat_map_from(bin_corr_fut_maxSSS,  bin_hyb_fut_maxSSS)
terra::writeRaster(pres_cat, file.path(OUTDIR, "rasters/categorical_present.tif"),
                   overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(fut_cat,  file.path(OUTDIR, "rasters/categorical_future.tif"),
                   overwrite = TRUE, datatype = "INT1U")


# ============================================================================ #
# Step 7. Niche overlap (current and future)                                   #
# ============================================================================ #

cat("\n=== STEP 7: Niche overlap (Schoener's D, Pearson r) ===\n")

D_current <- schoeners_d(pred_corr_curr$r_mean, pred_hyb_curr$r_mean)
D_future  <- schoeners_d(pred_corr_fut$r_mean,  pred_hyb_fut$r_mean)
r_current <- pearson_r(pred_corr_curr$r_mean, pred_hyb_curr$r_mean)
r_future  <- pearson_r(pred_corr_fut$r_mean,  pred_hyb_fut$r_mean)

cat(sprintf("  Schoener's D (current): %.4f\n", D_current))
cat(sprintf("  Schoener's D (future):  %.4f\n", D_future))
cat(sprintf("  Pearson r (current):    %.4f\n", r_current))
cat(sprintf("  Pearson r (future):     %.4f\n", r_future))

overlap_df <- data.frame(Scenario    = c("Current", "Future"),
                         Schoeners_D = c(D_current, D_future),
                         Pearson_r   = c(r_current, r_future))
readr::write_csv(overlap_df, file.path(OUTDIR, "niche_overlap_D_r.csv"))


# ============================================================================ #
# Step 8. Continuous difference maps                                           #
# ============================================================================ #

diff_present_cont <- pred_hyb_curr$r_mean - pred_corr_curr$r_mean
diff_future_cont  <- pred_hyb_fut$r_mean  - pred_corr_fut$r_mean
names(diff_present_cont) <- "diff_hybrid_minus_correlative_present"
names(diff_future_cont)  <- "diff_hybrid_minus_correlative_future"

terra::writeRaster(diff_present_cont,
                   file.path(OUTDIR, "rasters/diff_continuous_present.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(diff_future_cont,
                   file.path(OUTDIR, "rasters/diff_continuous_future.tif"),
                   overwrite = TRUE, datatype = "FLT4S")


# ============================================================================ #
# Step 9. Range-change metrics (P10 and maxSSS)                                #
# ============================================================================ #

cat("\n=== STEP 9: Range-change metrics ===\n")

calculate_range_metrics <- function(current, future, model_name) {
  cell_area_km2 <- prod(terra::res(current)) * 111 ^ 2
  stable <- current & future
  gained <- (!current) & future
  lost   <- current   & (!future)
  area_current <- terra::global(current, "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_future  <- terra::global(future,  "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_stable  <- terra::global(stable,  "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_gained  <- terra::global(gained,  "sum", na.rm = TRUE)[[1]] * cell_area_km2
  area_lost    <- terra::global(lost,    "sum", na.rm = TRUE)[[1]] * cell_area_km2
  jaccard      <- area_stable / (area_current + area_future - area_stable)
  data.frame(Model           = model_name,
             Current_km2     = area_current,
             Future_km2      = area_future,
             Stable_km2      = area_stable,
             Gained_km2      = area_gained,
             Lost_km2        = area_lost,
             Net_change_km2  = area_future - area_current,
             Percent_change  = 100 * (area_future - area_current) / area_current,
             Jaccard         = jaccard)
}

metrics_corr_p10 <- calculate_range_metrics(bin_corr_curr_p10, bin_corr_fut_p10, "Correlative_P10")
metrics_hyb_p10  <- calculate_range_metrics(bin_hyb_curr_p10,  bin_hyb_fut_p10,  "Hybrid_P10")
area_change_p10  <- dplyr::bind_rows(metrics_corr_p10, metrics_hyb_p10)
readr::write_csv(area_change_p10,
                 file.path(OUTDIR, "range_change_metrics_P10.csv"))

metrics_corr_maxSSS <- calculate_range_metrics(bin_corr_curr_maxSSS, bin_corr_fut_maxSSS, "Correlative")
metrics_hyb_maxSSS  <- calculate_range_metrics(bin_hyb_curr_maxSSS,  bin_hyb_fut_maxSSS,  "Hybrid")
area_change_maxSSS  <- dplyr::bind_rows(metrics_corr_maxSSS, metrics_hyb_maxSSS)
readr::write_csv(area_change_maxSSS,
                 file.path(OUTDIR, "range_change_metrics.csv"))


# ============================================================================ #
# Step 10. Save and finish                                                     #
# ============================================================================ #

save(pred_corr_curr, pred_hyb_curr, pred_corr_fut, pred_hyb_fut,
     bin_corr_curr_p10, bin_corr_fut_p10, bin_hyb_curr_p10, bin_hyb_fut_p10,
     bin_corr_curr_maxSSS, bin_corr_fut_maxSSS, bin_hyb_curr_maxSSS, bin_hyb_fut_maxSSS,
     pred_correlativo_fut, pred_hibrido_fut, meandb2_fut, meanActivity2_fut,
     file = file.path(OUTDIR, "models/05_projections.RData"))

save(D_current, D_future, r_current, r_future, overlap_df,
     diff_present_cont, diff_future_cont,
     area_change_p10, area_change_maxSSS,
     file = file.path(OUTDIR, "models/06_comparison.RData"))

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "session_info_projections_future.txt"))

cat("\n  OK: RCode_SDM_projections_future.R complete\n")
