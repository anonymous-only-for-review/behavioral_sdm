################################################################################
## RCode_SDM_correlative.R
##
## Species Distribution Models: correlative vs hybrid (correlative + behavioural
## mechanistic predictors) for Psammodromus algirus, fitted with Boosted
## Regression Trees, evaluated under repeated spatial cross-validation, and
## compared with a Wilcoxon signed-rank test.
##
## Pipeline steps:
##   1. Acquire climate data (WorldClim 2.1 baseline + CMIP6 future)
##   2. Load and spatially-thin GBIF occurrences
##   3. Load environmental and mechanistic predictors and align extents
##   4. Define accessible-area buffers (M_core, M_bg)
##   5. Generate bias-corrected background points (tempered_uniform)
##   6. Variable selection (Dormann's select07 + VIF)
##   7. Fit simple BRT ensembles (1 model per spatial fold)
##   8. Fit robust BRT ensembles (10 repeats x 5 spatial folds = 50 fits)
##   9. Compute thresholds (P10 and maxSSS) and held-out skill metrics
##  10. Wilcoxon signed-rank test on paired AUC and TSS across 50 fits
##
## Inputs (relative to repo root):
##   - SDM data/Psammodromus_GBIF_cleaned.RData    GBIF cleaned occurrences
##   - Sources/meandb_map.grd                      mechanistic deviation_mean
##   - Sources/meanActivity_map.grd                mechanistic Activity
##   - data/wc2.1_2.5m/...                         WorldClim 2.1 (auto-download)
##   - data/cmip6/...                              CMIP6 SSP2-4.5 (auto-download)
##
## Outputs (under OUTDIR, defined in RCode_SDM_helpers.R):
##   - models/01_data_loading.RData
##   - models/02_variable_selection.RData
##   - models/03_model_fitting.RData
##   - models/04_evaluation.RData
##   - models/selected_variables_correlative.rds
##   - models/selected_variables_hybrid.rds
##   - results_summary/thresholds_maxSSS_ensemble.rds
##   - model_performance_by_fold_50fits.csv
##   - model_performance_comparison_wilcoxon.csv
##
## Author: Guillermo Fandos (gfandos@ucm.es)
## Reproducibility: requires R >= 4.3, see renv.lock for exact package versions.
################################################################################

source("RCode_SDM_helpers.R")
set.seed(MASTER_SEED)


# ============================================================================ #
# Local input paths (override before sourcing if your layout differs)          #
# ============================================================================ #

GBIF_CLEANED_PATH <- "SDM data/Psammodromus_GBIF_cleaned.RData"

# Mechanistic layers from RCode_biophysical_model.R (Juanvi). For the
# correlative-vs-hybrid comparison we use the current-climate layers that
# Juanvi distributes in Sources/. If a derived/updated version is needed,
# override these paths.
MECH_CURR_MEANDB_PATH    <- file.path("Sources", "meandb_map.grd")
MECH_CURR_ACTIVITY_PATH  <- file.path("Sources", "meanActivity_map.grd")


# ============================================================================ #
# Step 1. Climate data acquisition                                             #
# ============================================================================ #
#
# WorldClim 2.1 (Fick & Hijmans 2017) and CMIP6 (Eyring et al. 2016) data are
# not redistributed here. Run this block once on a fresh checkout to populate
# data/. Subsequent runs read from disk.

cat("\n=== STEP 1: Climate data acquisition ===\n")

dir.create("data", showWarnings = FALSE, recursive = TRUE)

# Current bioclim, 2.5 arc-min
bio_curr <- geodata::worldclim_global(var  = "bio",
                                      res  = 2.5,
                                      path = "data")
names(bio_curr) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))

# Future bioclim, CMIP6 MPI-ESM1-2-HR, SSP2-4.5, 2041-2060, 2.5 arc-min
# (Decision rationale documented in RCode_SDM_projections_future.R)
bio_fut <- geodata::cmip6_world(model = "MPI-ESM1-2-HR",
                                ssp   = "245",
                                time  = "2041-2060",
                                var   = "bioc",
                                res   = 2.5,
                                path  = "data/cmip6")
names(bio_fut) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))


# ============================================================================ #
# Step 2. GBIF occurrences and spatial thinning                                #
# ============================================================================ #

cat("\n=== STEP 2: Loading and thinning GBIF occurrences ===\n")

load(GBIF_CLEANED_PATH)
stopifnot(exists("gbif_lizard_cleaned"))

cat(sprintf("  * Applying spatial thinning (min. distance = %d km)...\n",
            THIN_DIST_KM))
thinned <- spThin::thin(
  loc.data                 = gbif_lizard_cleaned,
  lat.col                  = "decimalLatitude",
  long.col                 = "decimalLongitude",
  spec.col                 = "species",
  thin.par                 = THIN_DIST_KM,
  reps                     = 1,
  locs.thinned.list.return = TRUE,
  write.files              = FALSE,
  verbose                  = FALSE
)

p_coords <- thinned[[1]] |>
  dplyr::rename(decimalLongitude = Longitude,
                decimalLatitude  = Latitude) |>
  dplyr::mutate(pr_ab = 1)

cat(sprintf("  * Retained %d/%d occurrences after thinning (%.1f%%)\n",
            nrow(p_coords), nrow(gbif_lizard_cleaned),
            100 * nrow(p_coords) / nrow(gbif_lizard_cleaned)))


# ============================================================================ #
# Step 3. Mechanistic predictors and extent alignment                          #
# ============================================================================ #

cat("\n=== STEP 3: Loading mechanistic predictors and aligning extents ===\n")

meandb2 <- terra::rast(MECH_CURR_MEANDB_PATH)
names(meandb2) <- "deviation_mean"

meanActivity2 <- terra::rast(MECH_CURR_ACTIVITY_PATH)
names(meanActivity2) <- "Activity"

# Resample bioclim to the mechanistic grid and mask to its extent
bio_curr <- terra::resample(bio_curr, meandb2) |> terra::mask(meandb2)

# Land mask from rnaturalearth
world_sf  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_v   <- terra::vect(world_sf) |>
  terra::project(terra::crs(meandb2)) |>
  terra::crop(terra::ext(meandb2))
land_mask <- terra::rasterize(world_v, meandb2[[1]], field = 1)

if (!terra::compareGeom(bio_curr[[1]], meandb2, stopOnError = FALSE)) {
  stop("Mechanistic layers extent doesn't match bioclim after resampling.")
}
cat("  OK: mechanistic layers aligned with bioclim grid\n")


# ============================================================================ #
# Step 4. Accessible-area buffers (M_core, M_bg)                               #
# ============================================================================ #

cat("\n=== STEP 4: Defining accessible-area buffers ===\n")

crs_r <- terra::crs(meandb2)
pres_sf        <- sf::st_as_sf(p_coords,
                               coords = c("decimalLongitude", "decimalLatitude"),
                               crs    = 4326)
world_sf_equal <- sf::st_transform(world_sf, CRS_EQUAL_AREA)
pres_equal     <- sf::st_transform(pres_sf,  CRS_EQUAL_AREA)

M_core_equal <- pres_equal |>
  sf::st_union() |>
  sf::st_buffer(dist = BUF_CORE_KM * 1000) |>
  sf::st_make_valid()

M_bg_equal <- pres_equal |>
  sf::st_union() |>
  sf::st_buffer(dist = (BUF_CORE_KM + BUF_BG_ADD_KM) * 1000) |>
  sf::st_make_valid()

land_equal <- sf::st_union(world_sf_equal) |> sf::st_make_valid()
M_core_land_equal <- suppressWarnings(sf::st_intersection(M_core_equal, land_equal)) |>
  sf::st_make_valid()
M_bg_land_equal   <- suppressWarnings(sf::st_intersection(M_bg_equal,   land_equal)) |>
  sf::st_make_valid()

M_core_wgs <- sf::st_transform(M_core_land_equal, crs = crs_r)
M_bg_wgs   <- sf::st_transform(M_bg_land_equal,   crs = crs_r)

r_ref    <- bio_curr[[1]]
M_core_r <- terra::rasterize(terra::vect(M_core_wgs), r_ref, field = 1, background = NA) |>
  terra::mask(land_mask)
M_bg_r   <- terra::rasterize(terra::vect(M_bg_wgs),   r_ref, field = 1, background = NA) |>
  terra::mask(land_mask)


# ============================================================================ #
# Step 5. Bias-corrected background sampling (tempered_uniform)                #
# ============================================================================ #

cat("\n=== STEP 5: Generating bias-corrected background points ===\n")

occ_sf   <- sf::st_as_sf(p_coords,
                         coords = c("decimalLongitude", "decimalLatitude"),
                         crs    = sf::st_crs(4326))
occ_rast <- terra::rasterize(terra::vect(occ_sf), bio_curr[[1]], fun = "length")
occ_rast[is.na(occ_rast)] <- 0
bias_surface <- terra::focal(occ_rast, w = 5, fun = "mean", na.rm = TRUE)

domain_mask <- switch(BG_DOMAIN,
                     "previous" = land_mask,
                     "M_bg"     = M_bg_r,
                     "M_core"   = M_core_r,
                     { warning("BG_DOMAIN not recognized; using 'previous'"); land_mask })

w <- terra::mask(bias_surface, domain_mask)
min_pos <- suppressWarnings(min(terra::values(w), na.rm = TRUE))
if (is.finite(min_pos)) {
  w[is.na(w) & !is.na(domain_mask)] <- min_pos
} else {
  w <- domain_mask
  w[!is.na(w)] <- 1
}

build_weights <- function(w, mode = "plain", gamma = 1, lambda = 0.0) {
  w_use <- w
  if (mode %in% c("tempered", "tempered_uniform", "stratified_rings")) {
    w_use <- w_use ^ gamma
  }
  if (mode %in% c("tempered_uniform", "stratified_rings")) {
    uniform <- domain_mask
    uniform[!is.na(uniform)] <- 1
    uniform <- uniform / terra::global(uniform, "sum", na.rm = TRUE)[[1]]
    w_use <- (1 - lambda) * w_use + lambda * uniform
  }
  w_use / terra::global(w_use, "sum", na.rm = TRUE)[[1]]
}

w_final <- build_weights(w, mode = BG_SAMPLING_MODE, gamma = GAMMA, lambda = LAMBDA)

set.seed(MASTER_SEED)
bg_cells <- terra::spatSample(w_final, size = N_BACKGROUND,
                              method   = "weights",
                              replace  = FALSE,
                              na.rm    = TRUE,
                              as.points = FALSE,
                              cells    = TRUE)
bg_coords <- terra::xyFromCell(bio_curr[[1]], bg_cells[, "cell"])
bg <- data.frame(decimalLongitude = bg_coords[, 1],
                 decimalLatitude  = bg_coords[, 2],
                 pr_ab            = 0)

presence_pa     <- dplyr::bind_rows(p_coords, bg)
presence_pa$id  <- seq_len(nrow(presence_pa))

cat(sprintf("  * BG_DOMAIN='%s' | MODE='%s' -> %d presences + %d background = %d points\n",
            BG_DOMAIN, BG_SAMPLING_MODE,
            sum(presence_pa$pr_ab == 1),
            sum(presence_pa$pr_ab == 0),
            nrow(presence_pa)))

# Drop points with NA in any predictor
pred_stack_check <- c(bio_curr, meandb2, meanActivity2)
vals_pa <- terra::extract(pred_stack_check,
                          presence_pa[, c("decimalLongitude", "decimalLatitude")])
ok_pa <- stats::complete.cases(vals_pa)
if (!all(ok_pa)) {
  n_drop      <- sum(!ok_pa)
  presence_pa <- presence_pa[ok_pa, , drop = FALSE]
  cat(sprintf("  * Dropped %d points with NA predictors (kept %d)\n",
              n_drop, nrow(presence_pa)))
}

# QC map of buffers + occurrences + background
M_core_plot  <- sf::st_as_sf(terra::as.polygons(M_core_r))
M_bg_plot    <- sf::st_as_sf(terra::as.polygons(M_bg_r))
sf::st_crs(M_core_plot) <- sf::st_crs(crs_r)
sf::st_crs(M_bg_plot)   <- sf::st_crs(crs_r)
pres_sf_plot <- sf::st_transform(pres_sf, crs_r)
bg_sf_plot   <- sf::st_as_sf(bg, coords = c("decimalLongitude", "decimalLatitude"), crs = crs_r)

ex <- terra::ext(r_ref)
qc_map <- ggplot() +
  geom_sf(data = sf::st_transform(world_sf, crs_r),
          fill = "grey97", color = "white", linewidth = 0.2) +
  geom_sf(data = M_core_plot, fill = NA, color = "#1E8449", linewidth = 0.6) +
  geom_sf(data = bg_sf_plot,   color = "#6E6E6E", alpha = 0.20, size = 0.25, show.legend = FALSE) +
  geom_sf(data = pres_sf_plot, color = "#F6A5A5", fill = "#F6A5A5",
          shape = 21, size = 1.2, stroke = 0.2) +
  coord_sf(xlim = c(ex[1], ex[2]), ylim = c(ex[3], ex[4]), expand = FALSE, clip = "on") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank())

ggsave(file.path(OUTDIR, "figures/qc_buffers_and_points.pdf"),
       qc_map, width = 160, height = 110, units = "mm", dpi = 300)

save(gbif_lizard_cleaned, p_coords, presence_pa,
     bio_curr, bio_fut, meandb2, meanActivity2,
     land_mask, M_core_r, M_bg_r, world_sf, crs_r,
     file = file.path(OUTDIR, "models/01_data_loading.RData"))


# ============================================================================ #
# Step 6. Variable selection (select07 + VIF)                                  #
# ============================================================================ #

cat("\n=== STEP 6: Variable selection (select07 + VIF) ===\n")

pred_full <- bio_curr

presence_pa_sf <- sf::st_as_sf(presence_pa,
                               coords = c("decimalLongitude", "decimalLatitude"),
                               crs    = terra::crs(pred_full))

data_full <- presence_pa_sf |>
  cbind(terra::extract(pred_full, _)) |>
  sf::st_drop_geometry() |>
  dplyr::select(-ID, -id) |>
  tidyr::drop_na()

sel <- select07(X         = data_full[, setdiff(names(data_full), "pr_ab")],
                y         = data_full$pr_ab,
                threshold = 0.7,
                method    = "spearman",
                univar    = "glm2")
vars_sel <- sel$pred_sel
cat("  * Variables after select07:\n")
print(vars_sel)

# Iterative VIF filter (drops worst until all VIFs <= thr)
vif_quick <- function(rst, thr = 5) {
  df   <- as.data.frame(rst, na.rm = TRUE)
  keep <- names(df)
  repeat {
    cm   <- tryCatch(cor(df, use = "pairwise.complete.obs"), error = function(e) NULL)
    if (is.null(cm)) break
    invR <- tryCatch(solve(cm), error = function(e) NULL)
    if (is.null(invR)) break
    v <- diag(invR)
    if (max(v, na.rm = TRUE) <= thr) break
    drop_var <- names(which.max(v))
    keep <- setdiff(keep, drop_var)
    df   <- df[, keep, drop = FALSE]
  }
  subset(rst, keep)
}

pred_correlativo0 <- subset(pred_full, vars_sel)
pred_correlativo  <- vif_quick(pred_correlativo0, thr = 5)
cat("  * Variables after VIF (correlative):\n")
print(names(pred_correlativo))

pred_hibrido0 <- c(pred_correlativo, meandb2, meanActivity2)
pred_hibrido  <- vif_quick(pred_hibrido0, thr = 5)
if (!"deviation_mean" %in% names(pred_hibrido)) pred_hibrido <- c(pred_hibrido, meandb2)
if (!"Activity"       %in% names(pred_hibrido)) pred_hibrido <- c(pred_hibrido, meanActivity2)
cat("  * Variables after VIF (hybrid):\n")
print(names(pred_hibrido))

vars_correlativo <- names(pred_correlativo)
vars_hibrido     <- names(pred_hibrido)

saveRDS(vars_correlativo, file.path(OUTDIR, "models/selected_variables_correlative.rds"))
saveRDS(vars_hibrido,     file.path(OUTDIR, "models/selected_variables_hybrid.rds"))

pred_correlativo <- terra::mask(pred_correlativo, land_mask, maskvalues = 0)
pred_hibrido     <- terra::mask(pred_hibrido,     land_mask, maskvalues = 0)

save(pred_correlativo, pred_hibrido, vars_correlativo, vars_hibrido,
     file = file.path(OUTDIR, "models/02_variable_selection.RData"))


# ============================================================================ #
# Step 7. Simple BRT ensembles (one model per spatial fold)                    #
# ============================================================================ #

cat("\n=== STEP 7: Fitting simple BRT ensembles (1 model per fold) ===\n")

mk_data <- function(pred) {
  pts_sf <- sf::st_as_sf(presence_pa,
                         coords = c("decimalLongitude", "decimalLatitude"),
                         crs    = terra::crs(pred))
  ext <- terra::extract(pred, pts_sf, cells = TRUE)
  out <- cbind(presence_pa, ext) |> tidyr::drop_na()
  out <- dplyr::distinct(out, cell, .keep_all = TRUE)
  out[, c("id", "pr_ab", "cell", names(pred))]
}

dat_corr <- mk_data(pred_correlativo)
dat_hyb  <- mk_data(pred_hibrido)
cat(sprintf("  * Correlative dataset: %d points\n", nrow(dat_corr)))
cat(sprintf("  * Hybrid dataset:      %d points\n", nrow(dat_hyb)))

coords_sf <- sf::st_as_sf(presence_pa,
                          coords = c("decimalLongitude", "decimalLatitude"),
                          crs    = terra::crs(pred_correlativo)) |>
  sf::st_transform(3857)
set.seed(MASTER_SEED)
sb        <- blockCV::cv_spatial(x         = coords_sf,
                                 column    = "pr_ab",
                                 k         = K_FOLDS,
                                 size      = 30 * 1000,
                                 selection = "systematic",
                                 progress  = FALSE)
fold_vec <- sb$folds_ids

dat_corr$fold <- fold_vec[match(dat_corr$id, presence_pa$id)]
dat_hyb$fold  <- fold_vec[match(dat_hyb$id,  presence_pa$id)]

fit_ensemble_gbm_simple <- function(data, predictors,
                                    response = "pr_ab", folds_col = "fold",
                                    n_trees           = N_TREES_MAX,
                                    interaction_depth = INTERACTION_DEPTH,
                                    shrinkage         = SHRINKAGE,
                                    bag_fraction      = BAG_FRACTION) {
  dat <- data |>
    dplyr::select(dplyr::all_of(c(predictors, response, folds_col))) |>
    tidyr::drop_na()
  fold_ids <- unique(dat[[folds_col]])
  models   <- vector("list", length(fold_ids))
  perf     <- data.frame()
  var_imp  <- list()

  for (i in seq_along(fold_ids)) {
    train_data <- dat |> dplyr::filter(.data[[folds_col]] != fold_ids[i])
    test_data  <- dat |> dplyr::filter(.data[[folds_col]] == fold_ids[i])
    m <- gbm::gbm(formula           = as.formula(paste(response, "~ .")),
                  data              = train_data |> dplyr::select(-!!rlang::sym(folds_col)),
                  distribution      = "bernoulli",
                  n.trees           = n_trees,
                  interaction.depth = interaction_depth,
                  shrinkage         = shrinkage,
                  bag.fraction      = bag_fraction,
                  cv.folds          = 0,
                  verbose           = FALSE,
                  n.cores           = 1)

    seq_nt   <- seq(100, n_trees, by = 100)
    pred_seq <- predict(m, newdata = test_data, n.trees = seq_nt, type = "response")
    auc_seq  <- apply(pred_seq, 2, function(p)
                      suppressMessages(pROC::auc(test_data[[response]], p)))
    best_nt  <- seq_nt[which.max(auc_seq)]

    models[[i]] <- m
    models[[i]]$best_ntree <- best_nt

    pred_test  <- predict(m, newdata = test_data,  n.trees = best_nt, type = "response")
    pred_train <- predict(m, newdata = train_data, n.trees = best_nt, type = "response")
    perf <- dplyr::bind_rows(perf,
      data.frame(fold = i, partition = "Test",
                 AUC = as.numeric(pROC::auc(test_data[[response]],  pred_test)),
                 n_trees = best_nt),
      data.frame(fold = i, partition = "Train",
                 AUC = as.numeric(pROC::auc(train_data[[response]], pred_train)),
                 n_trees = best_nt))
    vi          <- summary(m, n.trees = best_nt, plotit = FALSE)
    var_imp[[i]] <- as.data.frame(vi) |> dplyr::mutate(fold = i)
  }

  var_imp_combined <- dplyr::bind_rows(var_imp) |>
    dplyr::group_by(var) |>
    dplyr::summarise(RelInf_mean = mean(rel.inf),
                     RelInf_sd   = sd(rel.inf),
                     .groups     = "drop") |>
    dplyr::arrange(dplyr::desc(RelInf_mean)) |>
    dplyr::rename(Variable = var)

  list(models         = models,
       performance    = perf,
       var_importance = var_imp_combined,
       predictors     = predictors,
       response       = response)
}

mod_corr <- fit_ensemble_gbm_simple(dat_corr, predictors = names(pred_correlativo))
mod_hyb  <- fit_ensemble_gbm_simple(dat_hyb,  predictors = names(pred_hibrido))


# ============================================================================ #
# Step 8. Robust BRT: 10 repeats x 5 spatial folds = 50 fits                   #
# ============================================================================ #

cat("\n=== STEP 8: Repeated spatial CV BRT (50 fits) ===\n")

make_spatial_folds <- function(presence_pa, pred_ref,
                               k = K_FOLDS, block_km = 30, seed = 1L) {
  set.seed(seed)
  coords_sf <- sf::st_as_sf(presence_pa,
                            coords = c("decimalLongitude", "decimalLatitude"),
                            crs    = terra::crs(pred_ref)) |>
    sf::st_transform(3857)
  sb <- blockCV::cv_spatial(x         = coords_sf,
                            column    = "pr_ab",
                            k         = k,
                            size      = block_km * 1000,
                            selection = "random",
                            progress  = FALSE)
  sb$folds_ids
}

fit_ensemble_gbm_param <- function(data, predictors, folds_vec,
                                   shrinkage, interaction_depth,
                                   bag_fraction, n_trees,
                                   response = "pr_ab", folds_col = "fold") {
  dat <- data |>
    dplyr::select(dplyr::all_of(c(predictors, response, folds_col))) |>
    tidyr::drop_na()
  fold_ids <- sort(unique(folds_vec))
  models   <- vector("list", length(fold_ids))
  perf     <- vector("list", length(fold_ids))
  var_imp  <- vector("list", length(fold_ids))

  for (i in seq_along(fold_ids)) {
    te_id <- fold_ids[i]
    train_data <- dat |> dplyr::filter(.data[[folds_col]] != te_id)
    test_data  <- dat |> dplyr::filter(.data[[folds_col]] == te_id)
    w_train    <- make_class_weights(train_data[[response]])

    m <- gbm::gbm(formula           = stats::as.formula(paste(response, "~ .")),
                  data              = train_data |> dplyr::select(-!!rlang::sym(folds_col)),
                  weights           = w_train,
                  distribution      = "bernoulli",
                  n.trees           = n_trees,
                  interaction.depth = interaction_depth,
                  shrinkage         = shrinkage,
                  bag.fraction      = bag_fraction,
                  cv.folds          = 0,
                  keep.data         = FALSE,
                  verbose           = FALSE,
                  n.cores           = 1)

    seq_nt   <- seq(200, n_trees, by = 200)
    pred_seq <- predict(m, newdata = test_data, n.trees = seq_nt, type = "response")
    auc_seq  <- apply(pred_seq, 2, function(p)
                      suppressMessages(pROC::auc(test_data[[response]], p)))
    best_nt  <- seq_nt[which.max(auc_seq)]

    models[[i]] <- m
    models[[i]]$best_ntree <- best_nt

    pred_test  <- predict(m, newdata = test_data,  n.trees = best_nt, type = "response")
    pred_train <- predict(m, newdata = train_data, n.trees = best_nt, type = "response")

    obs_test <- test_data[[response]]
    thr_seq  <- seq(min(pred_test, na.rm = TRUE),
                    max(pred_test, na.rm = TRUE), length.out = 200)
    tss_vec  <- vapply(thr_seq, function(thr) {
      yhat <- as.numeric(pred_test >= thr)
      sens <- if (sum(obs_test == 1) > 0)
                sum(yhat == 1 & obs_test == 1) / sum(obs_test == 1) else NA_real_
      spec <- if (sum(obs_test == 0) > 0)
                sum(yhat == 0 & obs_test == 0) / sum(obs_test == 0) else NA_real_
      if (is.finite(sens) && is.finite(spec)) sens + spec - 1 else NA_real_
    }, numeric(1))
    tss_test <- max(tss_vec, na.rm = TRUE)

    perf[[i]] <- data.frame(fold     = i,
                            AUC_test = as.numeric(pROC::auc(obs_test, pred_test)),
                            TSS_test = tss_test,
                            AUC_train = as.numeric(pROC::auc(train_data[[response]], pred_train)),
                            n_trees   = best_nt)
    vi          <- summary(m, n.trees = best_nt, plotit = FALSE)
    var_imp[[i]] <- as.data.frame(vi) |> dplyr::mutate(fold = i)
  }

  perf_df <- dplyr::bind_rows(perf)
  var_imp_combined <- dplyr::bind_rows(var_imp) |>
    dplyr::group_by(var) |>
    dplyr::summarise(RelInf_mean = mean(rel.inf),
                     RelInf_sd   = sd(rel.inf),
                     .groups     = "drop") |>
    dplyr::arrange(dplyr::desc(RelInf_mean)) |>
    dplyr::rename(Variable = var)

  list(models         = models,
       performance    = perf_df,
       var_importance = var_imp_combined)
}

run_repeated_brt <- function(dat, predictors, label_model) {
  all_models    <- list()
  summary_table <- list()
  for (rep in 1:N_REPEATS) {
    cat(sprintf("\n[%s] Repeat %d/%d (seed=%d)\n",
                label_model, rep, N_REPEATS, MASTER_SEED + rep))
    folds_vec <- make_spatial_folds(presence_pa, pred_correlativo,
                                    k = K_FOLDS, block_km = 30,
                                    seed = MASTER_SEED + rep)
    dat$fold  <- folds_vec[match(dat$id, presence_pa$id)]
    fit       <- fit_ensemble_gbm_param(data              = dat,
                                        predictors        = predictors,
                                        folds_vec         = sort(unique(dat$fold)),
                                        shrinkage         = SHRINKAGE,
                                        interaction_depth = INTERACTION_DEPTH,
                                        bag_fraction      = BAG_FRACTION,
                                        n_trees           = N_TREES_MAX)
    all_models[[rep]]    <- fit$models
    summary_table[[rep]] <- cbind(Repeat = rep, Model = label_model, fit$performance)
  }
  list(best_models_list    = all_models,
       performance_by_fold = dplyr::bind_rows(summary_table))
}

rep_corr <- run_repeated_brt(dat_corr, names(pred_correlativo), "Correlative")
rep_hyb  <- run_repeated_brt(dat_hyb,  names(pred_hibrido),     "Hybrid")

readr::write_csv(rep_corr$performance_by_fold,
                 file.path(OUTDIR, "models/corr_perf_by_fold_repeats.csv"))
readr::write_csv(rep_hyb$performance_by_fold,
                 file.path(OUTDIR, "models/hyb_perf_by_fold_repeats.csv"))

# Flatten 10 repeats x 5 folds into a single super-ensemble of 50 models
flatten_models <- function(repeats_obj) {
  out <- list()
  for (i in seq_along(repeats_obj$best_models_list))
    for (j in seq_along(repeats_obj$best_models_list[[i]]))
      out[[length(out) + 1]] <- repeats_obj$best_models_list[[i]][[j]]
  out
}
super_models_corr <- flatten_models(rep_corr)
super_models_hyb  <- flatten_models(rep_hyb)

save(dat_corr, dat_hyb, mod_corr, mod_hyb, rep_corr, rep_hyb,
     super_models_corr, super_models_hyb, fold_vec,
     file = file.path(OUTDIR, "models/03_model_fitting.RData"))


# ============================================================================ #
# Step 9. Threshold selection and held-out skill metrics                       #
# ============================================================================ #

cat("\n=== STEP 9: Threshold selection (P10 + maxSSS) ===\n")

# maxTSS == maxSSS for binary y: TSS = sensitivity + specificity - 1.
apply_threshold_and_evaluate <- function(models, data,
                                         threshold_method = "p10",
                                         folds_col        = "fold") {
  fold_ids <- sort(unique(data[[folds_col]]))
  results <- lapply(seq_along(models), function(i) {
    m         <- models[[i]]
    fold_id   <- fold_ids[i]
    test_data <- data |> dplyr::filter(.data[[folds_col]] == fold_id)
    pred_test <- predict(m, newdata = test_data, n.trees = m$best_ntree, type = "response")

    if (threshold_method == "p10") {
      pres_preds <- pred_test[test_data$pr_ab == 1]
      threshold  <- stats::quantile(pres_preds, probs = 0.1, na.rm = TRUE)
    } else if (threshold_method %in% c("maxTSS", "maxSSS")) {
      obs     <- test_data$pr_ab
      seq_thr <- seq(min(pred_test, na.rm = TRUE),
                     max(pred_test, na.rm = TRUE), length.out = 100)
      tssv <- vapply(seq_thr, function(thr) {
        yhat <- as.numeric(pred_test >= thr)
        sens <- ifelse(sum(obs == 1) > 0, sum(yhat == 1 & obs == 1) / sum(obs == 1), NA_real_)
        spec <- ifelse(sum(obs == 0) > 0, sum(yhat == 0 & obs == 0) / sum(obs == 0), NA_real_)
        sens + spec - 1
      }, numeric(1))
      threshold <- seq_thr[which.max(tssv)]
    }

    cm   <- table(Observed  = test_data$pr_ab,
                  Predicted = as.numeric(pred_test >= threshold))
    sens <- if (all(dim(cm) == c(2, 2))) cm[2, 2] / sum(cm[2, ]) else NA_real_
    spec <- if (all(dim(cm) == c(2, 2))) cm[1, 1] / sum(cm[1, ]) else NA_real_
    tss  <- if (is.finite(sens) && is.finite(spec)) sens + spec - 1 else NA_real_
    auc  <- as.numeric(pROC::auc(test_data$pr_ab, pred_test))
    data.frame(fold        = i,
               fold_id     = fold_id,
               threshold   = threshold,
               AUC         = auc,
               TSS         = tss,
               Sensitivity = sens,
               Specificity = spec,
               Partition   = "Test")
  })
  eval_summary <- dplyr::bind_rows(results)
  list(mean_threshold = mean(eval_summary$threshold, na.rm = TRUE),
       eval_summary   = eval_summary)
}

thr_results_corr_p10 <- apply_threshold_and_evaluate(mod_corr$models, dat_corr, "p10")
thr_results_hyb_p10  <- apply_threshold_and_evaluate(mod_hyb$models,  dat_hyb,  "p10")

# maxSSS threshold via PresenceAbsence with method fallback
collect_test_preds <- function(models, data,
                               folds_col = "fold", response = "pr_ab") {
  fold_ids <- sort(unique(data[[folds_col]]))
  out <- vector("list", length(fold_ids))
  for (i in seq_along(fold_ids)) {
    m  <- models[[i]]
    te <- data[data[[folds_col]] == fold_ids[i], , drop = FALSE]
    pr <- as.numeric(predict(m, newdata = te, n.trees = m$best_ntree, type = "response"))
    out[[i]] <- data.frame(ID        = paste0("f", i, "_", seq_len(nrow(te))),
                           observed  = te[[response]],
                           predicted = pr)
  }
  do.call(rbind, out)
}

compute_maxSSS_threshold <- function(DATA) {
  try_methods <- c("MaxSens+Spec", "max.sensitivity+specificity", "Sens=Spec")
  for (mm in try_methods) {
    th <- try(PresenceAbsence::optimal.thresholds(DATA = DATA, opt.methods = mm),
              silent = TRUE)
    if (!inherits(th, "try-error") && "predicted" %in% names(th) && nrow(th) >= 1)
      return(as.numeric(th$predicted[1]))
  }
  ord <- order(DATA$predicted)
  p   <- DATA$predicted[ord]
  y   <- DATA$observed[ord]
  thr_seq   <- unique(p)
  sens_spec <- vapply(thr_seq, function(t) {
    yhat <- as.numeric(p >= t)
    sens <- ifelse(sum(y == 1) > 0, sum(yhat == 1 & y == 1) / sum(y == 1), NA_real_)
    spec <- ifelse(sum(y == 0) > 0, sum(yhat == 0 & y == 0) / sum(y == 0), NA_real_)
    sens + spec
  }, numeric(1))
  thr_seq[which.max(sens_spec)]
}

DATA_corr       <- collect_test_preds(mod_corr$models, dat_corr)
DATA_hyb        <- collect_test_preds(mod_hyb$models,  dat_hyb)
thr_corr_maxSSS <- compute_maxSSS_threshold(DATA_corr)
thr_hyb_maxSSS  <- compute_maxSSS_threshold(DATA_hyb)
cat(sprintf("  maxSSS threshold: corr=%.4f | hyb=%.4f\n",
            thr_corr_maxSSS, thr_hyb_maxSSS))


# ============================================================================ #
# Step 10. Wilcoxon signed-rank test (paired AUC and TSS, n = 50 fits)         #
# ============================================================================ #

cat("\n=== STEP 10: Wilcoxon signed-rank test (50 paired fits) ===\n")

perf_corr_50 <- rep_corr$performance_by_fold
perf_hyb_50  <- rep_hyb$performance_by_fold

comparison_50 <- dplyr::inner_join(
  perf_corr_50 |> dplyr::select(Repeat, fold, AUC_corr = AUC_test, TSS_corr = TSS_test),
  perf_hyb_50  |> dplyr::select(Repeat, fold, AUC_hyb  = AUC_test, TSS_hyb  = TSS_test),
  by = c("Repeat", "fold")
)

wilcox_auc <- wilcox.test(comparison_50$AUC_hyb, comparison_50$AUC_corr, paired = TRUE)
wilcox_tss <- wilcox.test(comparison_50$TSS_hyb, comparison_50$TSS_corr, paired = TRUE)

comparison_stats <- data.frame(
  Metric           = c("AUC", "TSS"),
  Correlative_mean = c(mean(comparison_50$AUC_corr, na.rm = TRUE),
                       mean(comparison_50$TSS_corr, na.rm = TRUE)),
  Correlative_sd   = c(sd(comparison_50$AUC_corr, na.rm = TRUE),
                       sd(comparison_50$TSS_corr, na.rm = TRUE)),
  Hybrid_mean      = c(mean(comparison_50$AUC_hyb, na.rm = TRUE),
                       mean(comparison_50$TSS_hyb, na.rm = TRUE)),
  Hybrid_sd        = c(sd(comparison_50$AUC_hyb, na.rm = TRUE),
                       sd(comparison_50$TSS_hyb, na.rm = TRUE)),
  Median_delta     = c(median(comparison_50$AUC_hyb - comparison_50$AUC_corr, na.rm = TRUE),
                       median(comparison_50$TSS_hyb - comparison_50$TSS_corr, na.rm = TRUE)),
  W_statistic      = c(wilcox_auc$statistic, wilcox_tss$statistic),
  p_value          = c(wilcox_auc$p.value,   wilcox_tss$p.value),
  n_fits           = c(nrow(comparison_50),  nrow(comparison_50))
)

readr::write_csv(comparison_50,
                 file.path(OUTDIR, "model_performance_by_fold_50fits.csv"))
readr::write_csv(comparison_stats,
                 file.path(OUTDIR, "model_performance_comparison_wilcoxon.csv"))
print(comparison_stats, digits = 4)

save(thr_results_corr_p10, thr_results_hyb_p10,
     thr_corr_maxSSS, thr_hyb_maxSSS,
     comparison_50, comparison_stats,
     file = file.path(OUTDIR, "models/04_evaluation.RData"))

# Persist maxSSS thresholds for downstream RCodes (integrated, projections)
saveRDS(list(correlativo_mean_ensemble = thr_corr_maxSSS,
             hibrido_mean_ensemble     = thr_hyb_maxSSS),
        file.path(OUTDIR, "results_summary/thresholds_maxSSS_ensemble.rds"))


# ============================================================================ #
# Reproducibility footer                                                       #
# ============================================================================ #

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "session_info_correlative.txt"))

cat("\n  OK: RCode_SDM_correlative.R complete\n")
