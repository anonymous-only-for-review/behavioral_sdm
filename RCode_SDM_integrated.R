################################################################################
## RCode_SDM_integrated.R
##
## Robustness analyses for the integrated (correlative + behavioural mechanistic)
## SDM in:
##
##   Rubalcaba, Fandos & Diaz. "Behavioural thermoregulation and microclimate
##   reshape climate-driven range forecasts."
##
## Three independent robustness analyses (run after RCode_SDM_correlative.R and
## RCode_SDM_projections_future.R have populated OUTDIR):
##
##   Section A. Prediction uncertainty (CV maps + MESS + manuscript values)
##              -- reads pre-computed rasters/CSVs from OUTDIR.
##
##   Section B. Null model test (reviewer DA-1 / R1-W1)
##              Demonstrates that hybrid-vs-correlative divergence is driven by
##              the mechanistic content, not by simply adding more predictors.
##              Three null-hybrid models with 2 spatially-structured but
##              non-mechanistic extra predictors:
##                - Null-A: topographic (TRI + slope)
##                - Null-B: quadratic climate (bio01^2 + bio12^2)
##                - Null-C: Gaussian random fields with same autocorrelation as
##                          deviation_mean (Krige-fitted spherical variogram).
##              Reduced 5 repeats x 5 folds = 25 fits per null hybrid.
##
##   Section C. Controlled comparison (reviewer R3)
##              Isolates the mechanistic contribution by fitting a matched-
##              correlative model that uses ONLY the bioclimatic variables
##              present in the hybrid. If divergence persists under matched
##              predictor sets, the mechanistic transformation is doing real
##              work; otherwise divergence was an artifact of predictor sets.
##              Full 10 repeats x 5 folds = 50 fits.
##
## Sections can be enabled/disabled at the top via RUN_A / RUN_B / RUN_C.
##
## Author: Guillermo Fandos (gfandos@ucm.es)
################################################################################

source("RCode_SDM_helpers.R")
set.seed(MASTER_SEED)

RUN_A <- TRUE   # Uncertainty + MESS + manuscript_values
RUN_B <- TRUE   # Null model test (5 reps)
RUN_C <- TRUE   # Controlled comparison (10 reps)

# Reduced repeats for the null model test (supplementary analysis)
N_REPEATS_NULL <- 5

# Local input paths (must match RCode_SDM_correlative.R)
GBIF_CLEANED_PATH        <- "SDM data/Psammodromus_GBIF_cleaned.RData"
MECH_CURR_MEANDB_PATH    <- file.path("Sources", "meandb_map.grd")
MECH_CURR_ACTIVITY_PATH  <- file.path("Sources", "meanActivity_map.grd")
MECH_FUT_WARM_MEANDB     <- file.path("Sources", "meandb_warm_map.grd")
MECH_FUT_WARM_ACTIVITY   <- file.path("Sources", "meanActivity_warm_map.grd")


# ============================================================================ #
# Shared utilities for Sections B and C                                        #
# ============================================================================ #

#' Rebuild presence/absence and predictor stacks from source files
#'
#' terra SpatRaster objects lose their C++ pointers when loaded from .RData,
#' so we always rebuild rasters from disk instead of relying on saved objects.
rebuild_pipeline_data <- function() {
  # Mechanistic layers (current)
  meandb2       <- terra::rast(MECH_CURR_MEANDB_PATH)
  names(meandb2) <- "deviation_mean"
  meanActivity2 <- terra::rast(MECH_CURR_ACTIVITY_PATH)
  names(meanActivity2) <- "Activity"

  # Bioclim (cached on disk under data/)
  bio_curr <- geodata::worldclim_global(var = "bio", res = 2.5,
                                        download = FALSE, path = "data")
  names(bio_curr) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
  bio_curr <- terra::resample(bio_curr, meandb2) |> terra::mask(meandb2)

  # Land mask
  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world_v  <- terra::vect(world_sf) |>
    terra::project(terra::crs(meandb2)) |>
    terra::crop(terra::ext(meandb2))
  land_mask <- terra::rasterize(world_v, meandb2[[1]], field = 1)

  # Occurrences (re-thinned with the same parameters)
  load(GBIF_CLEANED_PATH)
  thinned <- spThin::thin(loc.data                 = gbif_lizard_cleaned,
                          lat.col                  = "decimalLatitude",
                          long.col                 = "decimalLongitude",
                          spec.col                 = "species",
                          thin.par                 = THIN_DIST_KM,
                          reps                     = 1,
                          locs.thinned.list.return = TRUE,
                          write.files              = FALSE,
                          verbose                  = FALSE)
  p_coords <- thinned[[1]] |>
    dplyr::rename(decimalLongitude = Longitude,
                  decimalLatitude  = Latitude) |>
    dplyr::mutate(pr_ab = 1)

  # Buffers
  crs_r          <- terra::crs(meandb2)
  pres_sf        <- sf::st_as_sf(p_coords,
                                 coords = c("decimalLongitude", "decimalLatitude"),
                                 crs    = 4326)
  world_sf_equal <- sf::st_transform(world_sf, CRS_EQUAL_AREA)
  pres_equal     <- sf::st_transform(pres_sf,  CRS_EQUAL_AREA)
  M_core_equal   <- pres_equal |>
    sf::st_union() |>
    sf::st_buffer(dist = BUF_CORE_KM * 1000) |>
    sf::st_make_valid()
  land_equal        <- sf::st_union(world_sf_equal) |> sf::st_make_valid()
  M_core_land_equal <- suppressWarnings(sf::st_intersection(M_core_equal, land_equal)) |>
    sf::st_make_valid()
  M_core_wgs <- sf::st_transform(M_core_land_equal, crs = crs_r)
  M_core_r   <- terra::rasterize(terra::vect(M_core_wgs), bio_curr[[1]],
                                 field = 1, background = NA) |>
    terra::mask(land_mask)

  # Bias-corrected background
  occ_sf      <- sf::st_as_sf(p_coords,
                              coords = c("decimalLongitude", "decimalLatitude"),
                              crs    = 4326)
  occ_rast    <- terra::rasterize(terra::vect(occ_sf), bio_curr[[1]], fun = "length")
  occ_rast[is.na(occ_rast)] <- 0
  bias_surface <- terra::focal(occ_rast, w = 5, fun = "mean", na.rm = TRUE)
  domain_mask  <- M_core_r
  w <- terra::mask(bias_surface, domain_mask)
  min_pos <- suppressWarnings(min(terra::values(w), na.rm = TRUE))
  if (is.finite(min_pos)) {
    w[is.na(w) & !is.na(domain_mask)] <- min_pos
  } else {
    w <- domain_mask; w[!is.na(w)] <- 1
  }

  build_weights <- function(w, mode = "plain", gamma = 1, lambda = 0.0) {
    w_use <- w
    if (mode %in% c("tempered", "tempered_uniform")) w_use <- w_use ^ gamma
    if (mode == "tempered_uniform") {
      uniform <- domain_mask; uniform[!is.na(uniform)] <- 1
      uniform <- uniform / terra::global(uniform, "sum", na.rm = TRUE)[[1]]
      w_use <- (1 - lambda) * w_use + lambda * uniform
    }
    w_use / terra::global(w_use, "sum", na.rm = TRUE)[[1]]
  }
  w_final <- build_weights(w, mode = BG_SAMPLING_MODE, gamma = GAMMA, lambda = LAMBDA)

  set.seed(MASTER_SEED)
  bg_cells  <- terra::spatSample(w_final, size = N_BACKGROUND, method = "weights",
                                 replace = FALSE, na.rm = TRUE,
                                 as.points = FALSE, cells = TRUE)
  bg_coords <- terra::xyFromCell(bio_curr[[1]], bg_cells[, "cell"])
  bg <- data.frame(decimalLongitude = bg_coords[, 1],
                   decimalLatitude  = bg_coords[, 2],
                   pr_ab            = 0)
  presence_pa     <- dplyr::bind_rows(p_coords, bg)
  presence_pa$id  <- seq_len(nrow(presence_pa))
  pred_stack_check <- c(bio_curr, meandb2, meanActivity2)
  vals_pa <- terra::extract(pred_stack_check,
                            presence_pa[, c("decimalLongitude", "decimalLatitude")])
  ok_pa <- stats::complete.cases(vals_pa)
  if (!all(ok_pa)) presence_pa <- presence_pa[ok_pa, , drop = FALSE]

  list(bio_curr = bio_curr, meandb2 = meandb2, meanActivity2 = meanActivity2,
       land_mask = land_mask, M_core_r = M_core_r, presence_pa = presence_pa,
       p_coords = p_coords, crs_r = crs_r, world_sf = world_sf)
}

#' Build a matrix of (presence/absence, predictors) from a SpatRaster stack
mk_data <- function(presence_pa, pred) {
  pts_sf <- sf::st_as_sf(presence_pa,
                         coords = c("decimalLongitude", "decimalLatitude"),
                         crs    = terra::crs(pred))
  ext <- terra::extract(pred, pts_sf, cells = TRUE)
  out <- cbind(presence_pa, ext) |> tidyr::drop_na()
  out <- dplyr::distinct(out, cell, .keep_all = TRUE)
  out[, c("id", "pr_ab", "cell", names(pred))]
}

#' Spatial blockCV folds (random selection)
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

#' Fit a single ensemble (one BRT per fold) with class-balanced weights
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

    perf[[i]] <- data.frame(fold      = i,
                            AUC_test  = as.numeric(pROC::auc(obs_test, pred_test)),
                            TSS_test  = tss_test,
                            AUC_train = as.numeric(pROC::auc(train_data[[response]], pred_train)),
                            n_trees   = best_nt)
  }
  list(models = models, performance = dplyr::bind_rows(perf))
}

#' Repeated-CV BRT runner
run_repeated_brt <- function(dat, predictors, label_model, presence_pa, pred_correlativo,
                             n_repeats = N_REPEATS) {
  all_models    <- list()
  summary_table <- list()
  for (rep in 1:n_repeats) {
    cat(sprintf("\n[%s] Repeat %d/%d (seed=%d)\n",
                label_model, rep, n_repeats, MASTER_SEED + rep))
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

flatten_models <- function(repeats_obj) {
  out <- list()
  for (i in seq_along(repeats_obj$best_models_list))
    for (j in seq_along(repeats_obj$best_models_list[[i]]))
      out[[length(out) + 1]] <- repeats_obj$best_models_list[[i]][[j]]
  out
}

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

spearman_r <- function(r1, r2) {
  v1 <- terra::values(r1, mat = FALSE)
  v2 <- terra::values(r2, mat = FALSE)
  ok <- is.finite(v1) & is.finite(v2)
  stats::cor(v1[ok], v2[ok], method = "spearman")
}


# ============================================================================ #
# SECTION A. Prediction uncertainty + MESS + manuscript_values                 #
# ============================================================================ #

if (RUN_A) {
  cat("\n############################################################\n")
  cat("# SECTION A. Uncertainty + MESS + manuscript values\n")
  cat("############################################################\n\n")

  # CV-of-prediction maps
  r_corr_curr_mean <- terra::rast(file.path(OUTDIR, "rasters/current_correlative_mean.tif"))
  r_corr_curr_sd   <- terra::rast(file.path(OUTDIR, "rasters/current_correlative_sd.tif"))
  r_corr_fut_mean  <- terra::rast(file.path(OUTDIR, "rasters/future_correlative_mean.tif"))
  r_corr_fut_sd    <- terra::rast(file.path(OUTDIR, "rasters/future_correlative_sd.tif"))
  r_hyb_curr_mean  <- terra::rast(file.path(OUTDIR, "rasters/current_hybrid_mean.tif"))
  r_hyb_curr_sd    <- terra::rast(file.path(OUTDIR, "rasters/current_hybrid_sd.tif"))
  r_hyb_fut_mean   <- terra::rast(file.path(OUTDIR, "rasters/future_hybrid_mean.tif"))
  r_hyb_fut_sd     <- terra::rast(file.path(OUTDIR, "rasters/future_hybrid_sd.tif"))

  cv_corr_curr <- r_corr_curr_sd / r_corr_curr_mean
  cv_corr_fut  <- r_corr_fut_sd  / r_corr_fut_mean
  cv_hyb_curr  <- r_hyb_curr_sd  / r_hyb_curr_mean
  cv_hyb_fut   <- r_hyb_fut_sd   / r_hyb_fut_mean

  terra::writeRaster(cv_corr_curr, file.path(OUTDIR, "rasters/cv_current_correlative.tif"),
                     overwrite = TRUE, datatype = "FLT4S")
  terra::writeRaster(cv_corr_fut,  file.path(OUTDIR, "rasters/cv_future_correlative.tif"),
                     overwrite = TRUE, datatype = "FLT4S")
  terra::writeRaster(cv_hyb_curr,  file.path(OUTDIR, "rasters/cv_current_hybrid.tif"),
                     overwrite = TRUE, datatype = "FLT4S")
  terra::writeRaster(cv_hyb_fut,   file.path(OUTDIR, "rasters/cv_future_hybrid.tif"),
                     overwrite = TRUE, datatype = "FLT4S")

  summarize_cv <- function(cv_raster, desc) {
    v <- terra::values(cv_raster, mat = FALSE)
    v <- v[is.finite(v) & v >= 0]
    data.frame(Desc      = desc,
               N_layers  = NA,
               CV_median = median(v, na.rm = TRUE),
               CV_p75    = quantile(v, 0.75, na.rm = TRUE),
               CV_p90    = quantile(v, 0.90, na.rm = TRUE),
               CV_p95    = quantile(v, 0.95, na.rm = TRUE))
  }
  cv_summary <- rbind(
    summarize_cv(cv_corr_curr, "Correlative_current"),
    summarize_cv(cv_corr_fut,  "Correlative_future"),
    summarize_cv(cv_hyb_curr,  "Hybrid_current"),
    summarize_cv(cv_hyb_fut,   "Hybrid_future")
  )
  readr::write_csv(cv_summary,
                   file.path(OUTDIR, "results_summary/cv_spatial_summary.csv"))

  # MESS analysis (uses pre-computed rasters from projections script)
  r_mess_corr <- terra::rast(file.path(OUTDIR, "rasters/mess_correlative_future_pointsRef.tif"))
  r_mess_hyb  <- terra::rast(file.path(OUTDIR, "rasters/mess_hybrid_future_pointsRef.tif"))

  mess_vals_corr <- terra::values(r_mess_corr, mat = FALSE)
  mess_vals_hyb  <- terra::values(r_mess_hyb,  mat = FALSE)
  df_mess <- dplyr::bind_rows(
    data.frame(Model = "Correlative", value = mess_vals_corr),
    data.frame(Model = "Hybrid",      value = mess_vals_hyb)
  ) |> tidyr::drop_na()
  perc_extrap_all <- df_mess |>
    dplyr::mutate(Extrap = value < 0) |>
    dplyr::group_by(Model) |>
    dplyr::summarise(Extrap_pct = 100 * mean(Extrap), .groups = "drop")

  bin_corr_fut_maxSSS <- terra::rast(file.path(OUTDIR, "rasters/binary_future_correlative_maxSSS_PA.tif"))
  bin_hyb_fut_maxSSS  <- terra::rast(file.path(OUTDIR, "rasters/binary_future_hybrid_maxSSS_PA.tif"))
  bin_corr_curr_p10   <- terra::rast(file.path(OUTDIR, "rasters/binary_current_correlative_p10.tif"))
  bin_hyb_curr_p10    <- terra::rast(file.path(OUTDIR, "rasters/binary_current_hybrid_p10.tif"))

  compute_mess_in_range <- function(r_mess, r_bin_fut, r_bin_curr, model_label) {
    mess_fut <- terra::mask(r_mess, r_bin_fut, maskvalues = 0)
    vf <- terra::values(mess_fut, mat = FALSE); vf <- vf[is.finite(vf)]
    fut_lt0  <- if (length(vf) > 0) 100 * mean(vf < 0)   else NA
    fut_lt10 <- if (length(vf) > 0) 100 * mean(vf < -10) else NA
    mess_curr <- terra::mask(r_mess, r_bin_curr, maskvalues = 0)
    vc <- terra::values(mess_curr, mat = FALSE); vc <- vc[is.finite(vc)]
    curr_lt0  <- if (length(vc) > 0) 100 * mean(vc < 0) else NA
    data.frame(Model            = model_label,
               N_fut_cells       = length(vf),
               Fut_MESS_lt0_pct  = fut_lt0,
               Fut_MESS_lt10_pct = fut_lt10,
               N_curr_cells      = length(vc),
               Curr_MESS_lt0_pct = curr_lt0)
  }
  mess_range_df <- rbind(
    compute_mess_in_range(r_mess_corr, bin_corr_fut_maxSSS, bin_corr_curr_p10, "Correlative"),
    compute_mess_in_range(r_mess_hyb,  bin_hyb_fut_maxSSS,  bin_hyb_curr_p10,  "Hybrid")
  )
  readr::write_csv(perc_extrap_all,
                   file.path(OUTDIR, "mess_extrapolation_percentages_all_cells.csv"))
  readr::write_csv(mess_range_df,
                   file.path(OUTDIR, "mess_extrapolation_percentages.csv"))

  fig_mess <- plot_mess(r_mess_corr, "Correlative: MESS") +
    plot_mess(r_mess_hyb, "Hybrid: MESS") +
    patchwork::plot_layout(ncol = 2)
  ggsave(file.path(OUTDIR, "figures/fig8_mess_correlative_hybrid_future.pdf"),
         fig_mess, width = 180, height = 90, units = "mm", dpi = 300)

  # Categorical agreement, delta suitability summary, performance by repeat,
  # and consolidated manuscript_values.txt are produced by RCode_SDM_make_figures.R
  # (which consumes the same disk artefacts).
}


# ============================================================================ #
# SECTION B. Null model test                                                   #
# ============================================================================ #

if (RUN_B) {
  cat("\n############################################################\n")
  cat("# SECTION B. Null model test\n")
  cat("############################################################\n\n")

  if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat")
  library(gstat)

  NULL_DIR <- file.path(OUTDIR, "null_model_test")
  invisible(sapply(file.path(NULL_DIR, c("", "null_predictors", "models",
                                          "rasters", "figures")),
                   dir.create, showWarnings = FALSE, recursive = TRUE))

  pd <- rebuild_pipeline_data()
  bio_curr     <- pd$bio_curr
  meandb2      <- pd$meandb2
  meanActivity2 <- pd$meanActivity2
  land_mask    <- pd$land_mask
  presence_pa  <- pd$presence_pa

  # Reproduce the variable selection that the main pipeline uses
  pred_full       <- bio_curr
  presence_pa_sf  <- sf::st_as_sf(presence_pa,
                                  coords = c("decimalLongitude", "decimalLatitude"),
                                  crs    = terra::crs(pred_full))
  data_full <- presence_pa_sf |>
    cbind(terra::extract(pred_full, _)) |>
    sf::st_drop_geometry() |>
    dplyr::select(-ID, -id) |>
    tidyr::drop_na()

  sel       <- select07(X      = data_full[, setdiff(names(data_full), "pr_ab")],
                        y      = data_full$pr_ab,
                        threshold = 0.7,
                        method = "spearman",
                        univar = "glm2")
  vars_sel  <- sel$pred_sel

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

  pred_correlativo <- vif_quick(subset(pred_full, vars_sel), thr = 5)
  pred_hibrido0    <- c(pred_correlativo, meandb2, meanActivity2)
  pred_hibrido     <- vif_quick(pred_hibrido0, thr = 5)
  if (!"deviation_mean" %in% names(pred_hibrido)) pred_hibrido <- c(pred_hibrido, meandb2)
  if (!"Activity"       %in% names(pred_hibrido)) pred_hibrido <- c(pred_hibrido, meanActivity2)

  vars_correlativo <- names(pred_correlativo)
  vars_hibrido     <- names(pred_hibrido)

  pred_correlativo <- terra::mask(pred_correlativo, land_mask, maskvalues = 0)
  pred_hibrido     <- terra::mask(pred_hibrido,     land_mask, maskvalues = 0)

  # --- Reference correlative + hybrid ensembles (reduced repeats) --------------
  dat_corr <- mk_data(presence_pa, pred_correlativo)
  dat_hyb  <- mk_data(presence_pa, pred_hibrido)

  rep_corr_ref <- run_repeated_brt(dat_corr, vars_correlativo, "Correlative",
                                   presence_pa, pred_correlativo,
                                   n_repeats = N_REPEATS_NULL)
  rep_hyb_ref  <- run_repeated_brt(dat_hyb,  vars_hibrido,     "Hybrid",
                                   presence_pa, pred_correlativo,
                                   n_repeats = N_REPEATS_NULL)
  super_models_corr <- flatten_models(rep_corr_ref)
  super_models_hyb  <- flatten_models(rep_hyb_ref)

  pred_corr_curr <- predict_ensemble(super_models_corr, pred_correlativo)
  pred_hyb_curr  <- predict_ensemble(super_models_hyb,  pred_hibrido)
  pred_corr_curr$r_mean <- terra::mask(pred_corr_curr$r_mean, land_mask)
  pred_hyb_curr$r_mean  <- terra::mask(pred_hyb_curr$r_mean,  land_mask)

  # Future climate (preserved as in original 11_null_model_test: MIROC6/SSP585)
  bio_fut <- geodata::cmip6_world(model = "MIROC6", ssp = "585",
                                  time = "2041-2060", var = "bioc",
                                  res = 2.5, path = "data", download = FALSE)
  names(bio_fut) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
  bio_fut <- terra::resample(bio_fut, meandb2) |> terra::mask(meandb2)

  meandb2_fut       <- terra::rast(MECH_FUT_WARM_MEANDB)
  names(meandb2_fut) <- "deviation_mean"
  meanActivity2_fut <- terra::rast(MECH_FUT_WARM_ACTIVITY)
  names(meanActivity2_fut) <- "Activity"

  pred_correlativo_fut <- bio_fut[[vars_correlativo]] |>
    terra::mask(land_mask, maskvalues = 0)
  vars_hyb_clim    <- setdiff(vars_hibrido, c("deviation_mean", "Activity"))
  pred_hibrido_fut <- c(bio_fut[[vars_hyb_clim]], meandb2_fut, meanActivity2_fut) |>
    terra::mask(land_mask, maskvalues = 0)

  pred_corr_fut <- predict_ensemble(super_models_corr, pred_correlativo_fut)
  pred_hyb_fut  <- predict_ensemble(super_models_hyb,  pred_hibrido_fut)
  pred_corr_fut$r_mean <- terra::mask(pred_corr_fut$r_mean, land_mask)
  pred_hyb_fut$r_mean  <- terra::mask(pred_hyb_fut$r_mean,  land_mask)

  # --- Generate null predictors ------------------------------------------------
  rescale01 <- function(r) {
    v <- terra::values(r, mat = FALSE)
    mn <- min(v, na.rm = TRUE); mx <- max(v, na.rm = TRUE)
    if (mx == mn) return(r * 0)
    (r - mn) / (mx - mn)
  }

  # Null-A: topographic
  elev    <- geodata::elevation_global(res = 2.5, path = "data", download = TRUE)
  elev    <- terra::resample(elev, meandb2) |> terra::mask(meandb2)
  tri_r   <- rescale01(terra::terrain(elev, v = "TRI"));   names(tri_r)   <- "null_TRI"
  slope_r <- rescale01(terra::terrain(elev, v = "slope")); names(slope_r) <- "null_slope"
  terra::writeRaster(tri_r,   file.path(NULL_DIR, "null_predictors/null_TRI.tif"),   overwrite = TRUE)
  terra::writeRaster(slope_r, file.path(NULL_DIR, "null_predictors/null_slope.tif"), overwrite = TRUE)

  # Null-B: quadratic climate
  bio01_sq <- rescale01(bio_curr[["bio01"]] ^ 2); names(bio01_sq) <- "null_bio01sq"
  bio12_sq <- rescale01(bio_curr[["bio12"]] ^ 2); names(bio12_sq) <- "null_bio12sq"
  terra::writeRaster(bio01_sq, file.path(NULL_DIR, "null_predictors/null_bio01sq.tif"), overwrite = TRUE)
  terra::writeRaster(bio12_sq, file.path(NULL_DIR, "null_predictors/null_bio12sq.tif"), overwrite = TRUE)

  # Null-C: Gaussian random fields fitted to deviation_mean autocorrelation
  dev_df <- as.data.frame(meandb2, xy = TRUE, na.rm = TRUE); names(dev_df)[3] <- "z"
  set.seed(MASTER_SEED)
  idx_sub <- sample(nrow(dev_df), min(5000, nrow(dev_df)))
  dev_sub <- dev_df[idx_sub, ]
  sp::coordinates(dev_sub) <- ~ x + y
  vg_emp  <- gstat::variogram(z ~ 1, data = dev_sub, cutoff = 15, width = 0.5)
  vg_fit  <- gstat::fit.variogram(vg_emp,
                                  gstat::vgm(psill = var(dev_df$z, na.rm = TRUE),
                                              model = "Sph", range = 5, nugget = 0.01))

  grid_df <- as.data.frame(meandb2, xy = TRUE, na.rm = TRUE)
  grid_sp <- grid_df[, c("x", "y")]
  sp::coordinates(grid_sp) <- ~ x + y
  sp::gridded(grid_sp)     <- TRUE

  set.seed(MASTER_SEED + 999)
  g_null     <- gstat::gstat(formula = z ~ 1, dummy = TRUE, beta = 0,
                             model = vg_fit, nmax = 50)
  sim_fields <- predict(g_null, newdata = grid_sp, nsim = 2)
  null_rf1   <- terra::rasterize(terra::vect(sim_fields, geom = c("x", "y")),
                                 meandb2, field = "sim1") |> rescale01()
  null_rf2   <- terra::rasterize(terra::vect(sim_fields, geom = c("x", "y")),
                                 meandb2, field = "sim2") |> rescale01()
  names(null_rf1) <- "null_RF1"
  names(null_rf2) <- "null_RF2"
  terra::writeRaster(null_rf1, file.path(NULL_DIR, "null_predictors/null_RF1.tif"), overwrite = TRUE)
  terra::writeRaster(null_rf2, file.path(NULL_DIR, "null_predictors/null_RF2.tif"), overwrite = TRUE)

  # --- Null predictor stacks (with VIF filter, no forced retention) -----------
  pred_nullA <- vif_quick(c(pred_correlativo, tri_r,    slope_r),  thr = 5) |>
    terra::mask(land_mask, maskvalues = 0)
  pred_nullB <- vif_quick(c(pred_correlativo, bio01_sq, bio12_sq), thr = 5) |>
    terra::mask(land_mask, maskvalues = 0)
  pred_nullC <- vif_quick(c(pred_correlativo, null_rf1, null_rf2), thr = 5) |>
    terra::mask(land_mask, maskvalues = 0)

  null_survival <- data.frame(
    Model        = c("Null-A", "Null-A", "Null-B", "Null-B", "Null-C", "Null-C"),
    Predictor    = c("null_TRI", "null_slope", "null_bio01sq", "null_bio12sq",
                     "null_RF1", "null_RF2"),
    Survived_VIF = c("null_TRI"     %in% names(pred_nullA),
                     "null_slope"   %in% names(pred_nullA),
                     "null_bio01sq" %in% names(pred_nullB),
                     "null_bio12sq" %in% names(pred_nullB),
                     "null_RF1"     %in% names(pred_nullC),
                     "null_RF2"     %in% names(pred_nullC))
  )
  readr::write_csv(null_survival, file.path(NULL_DIR, "null_predictor_survival.csv"))

  # --- Fit the three null hybrids ---------------------------------------------
  dat_nullA <- mk_data(presence_pa, pred_nullA)
  dat_nullB <- mk_data(presence_pa, pred_nullB)
  dat_nullC <- mk_data(presence_pa, pred_nullC)

  rep_nullA <- run_repeated_brt(dat_nullA, names(pred_nullA), "Null-A_Topo",
                                presence_pa, pred_correlativo, n_repeats = N_REPEATS_NULL)
  rep_nullB <- run_repeated_brt(dat_nullB, names(pred_nullB), "Null-B_QuadClim",
                                presence_pa, pred_correlativo, n_repeats = N_REPEATS_NULL)
  rep_nullC <- run_repeated_brt(dat_nullC, names(pred_nullC), "Null-C_RandomRF",
                                presence_pa, pred_correlativo, n_repeats = N_REPEATS_NULL)
  super_nullA <- flatten_models(rep_nullA)
  super_nullB <- flatten_models(rep_nullB)
  super_nullC <- flatten_models(rep_nullC)

  perf_null <- dplyr::bind_rows(rep_nullA$performance_by_fold,
                                rep_nullB$performance_by_fold,
                                rep_nullC$performance_by_fold)
  readr::write_csv(perf_null, file.path(NULL_DIR, "models/null_model_performance.csv"))

  # --- Project null hybrids (current + future) --------------------------------
  bio01_sq_fut <- rescale01(bio_fut[["bio01"]] ^ 2); names(bio01_sq_fut) <- "null_bio01sq"
  bio12_sq_fut <- rescale01(bio_fut[["bio12"]] ^ 2); names(bio12_sq_fut) <- "null_bio12sq"

  build_fut_stack <- function(pred_curr, bio_fut_all, null_curr_layers, null_fut_layers = NULL) {
    clim_names <- intersect(names(pred_curr), names(bio_fut_all))
    null_names <- setdiff(names(pred_curr), names(bio_fut_all))
    fut_stack  <- bio_fut_all[[clim_names]]
    if (!is.null(null_fut_layers)) {
      fut_stack <- c(fut_stack, null_fut_layers[[null_names]])
    } else {
      fut_stack <- c(fut_stack, null_curr_layers[[null_names]])
    }
    terra::mask(fut_stack, land_mask, maskvalues = 0)
  }

  pred_nullA_fut <- build_fut_stack(pred_nullA, bio_fut, c(tri_r,    slope_r))
  pred_nullB_fut <- build_fut_stack(pred_nullB, bio_fut, NULL,
                                    null_fut_layers = c(bio01_sq_fut, bio12_sq_fut))
  pred_nullC_fut <- build_fut_stack(pred_nullC, bio_fut, c(null_rf1, null_rf2))

  pred_nullA_curr_ens <- predict_ensemble(super_nullA, pred_nullA)
  pred_nullA_fut_ens  <- predict_ensemble(super_nullA, pred_nullA_fut)
  pred_nullB_curr_ens <- predict_ensemble(super_nullB, pred_nullB)
  pred_nullB_fut_ens  <- predict_ensemble(super_nullB, pred_nullB_fut)
  pred_nullC_curr_ens <- predict_ensemble(super_nullC, pred_nullC)
  pred_nullC_fut_ens  <- predict_ensemble(super_nullC, pred_nullC_fut)

  for (obj_name in c("pred_nullA_curr_ens", "pred_nullA_fut_ens",
                     "pred_nullB_curr_ens", "pred_nullB_fut_ens",
                     "pred_nullC_curr_ens", "pred_nullC_fut_ens")) {
    obj <- get(obj_name)
    obj$r_mean <- terra::mask(obj$r_mean, land_mask)
    obj$r_sd   <- terra::mask(obj$r_sd,   land_mask)
    assign(obj_name, obj)
  }
  terra::writeRaster(pred_nullA_fut_ens$r_mean,
                     file.path(NULL_DIR, "rasters/nullA_future_mean.tif"), overwrite = TRUE)
  terra::writeRaster(pred_nullB_fut_ens$r_mean,
                     file.path(NULL_DIR, "rasters/nullB_future_mean.tif"), overwrite = TRUE)
  terra::writeRaster(pred_nullC_fut_ens$r_mean,
                     file.path(NULL_DIR, "rasters/nullC_future_mean.tif"), overwrite = TRUE)

  # --- Compare divergence -----------------------------------------------------
  compute_divergence <- function(null_curr, null_fut, corr_curr, corr_fut, label) {
    data.frame(Comparison = label,
               D_current   = schoeners_d(null_curr$r_mean, corr_curr$r_mean),
               D_future    = schoeners_d(null_fut$r_mean,  corr_fut$r_mean),
               r_current   = pearson_r(null_curr$r_mean, corr_curr$r_mean),
               r_future    = pearson_r(null_fut$r_mean,  corr_fut$r_mean),
               rho_current = spearman_r(null_curr$r_mean, corr_curr$r_mean),
               rho_future  = spearman_r(null_fut$r_mean,  corr_fut$r_mean))
  }

  div_real  <- compute_divergence(pred_hyb_curr,        pred_hyb_fut,
                                  pred_corr_curr, pred_corr_fut, "Real_Hybrid")
  div_nullA <- compute_divergence(pred_nullA_curr_ens, pred_nullA_fut_ens,
                                  pred_corr_curr, pred_corr_fut, "Null-A_Topo")
  div_nullB <- compute_divergence(pred_nullB_curr_ens, pred_nullB_fut_ens,
                                  pred_corr_curr, pred_corr_fut, "Null-B_QuadClim")
  div_nullC <- compute_divergence(pred_nullC_curr_ens, pred_nullC_fut_ens,
                                  pred_corr_curr, pred_corr_fut, "Null-C_RandomRF")

  diff_real  <- pred_hyb_fut$r_mean       - pred_corr_fut$r_mean
  diff_nullA <- pred_nullA_fut_ens$r_mean - pred_corr_fut$r_mean
  diff_nullB <- pred_nullB_fut_ens$r_mean - pred_corr_fut$r_mean
  diff_nullC <- pred_nullC_fut_ens$r_mean - pred_corr_fut$r_mean

  div_real$spatial_corr_with_real  <- 1.0
  div_nullA$spatial_corr_with_real <- pearson_r(diff_nullA, diff_real)
  div_nullB$spatial_corr_with_real <- pearson_r(diff_nullB, diff_real)
  div_nullC$spatial_corr_with_real <- pearson_r(diff_nullC, diff_real)

  div_real$corr_delta_vs_deviation  <- pearson_r(diff_real,  meandb2_fut)
  div_nullA$corr_delta_vs_deviation <- pearson_r(diff_nullA, meandb2_fut)
  div_nullB$corr_delta_vs_deviation <- pearson_r(diff_nullB, meandb2_fut)
  div_nullC$corr_delta_vs_deviation <- pearson_r(diff_nullC, meandb2_fut)

  null_summary <- dplyr::bind_rows(div_real, div_nullA, div_nullB, div_nullC)
  readr::write_csv(null_summary, file.path(NULL_DIR, "null_model_summary.csv"))
  print(null_summary, digits = 3)

  # --- Diagnostic figures -----------------------------------------------------
  diff_stack <- c(diff_real, diff_nullA, diff_nullB, diff_nullC)
  names(diff_stack) <- c("Real Hybrid", "Null-A (Topo)",
                         "Null-B (Quad Clim)", "Null-C (Random)")
  df_diff <- as.data.frame(diff_stack, xy = TRUE, na.rm = TRUE) |>
    tidyr::pivot_longer(cols = -c(x, y), names_to = "Model", values_to = "delta_suit")

  fig_delta_maps <- ggplot(df_diff, aes(x = x, y = y, fill = delta_suit)) +
    geom_raster() +
    facet_wrap(~ Model, ncol = 2) +
    scale_fill_gradient2(low = "#D73027", mid = "white", high = "#4575B4",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         oob = scales::squish,
                         name = expression(Delta * " suitability")) +
    coord_equal(expand = FALSE) +
    labs(title    = "Divergence from correlative model (future)",
         subtitle = "Null models vs. real hybrid") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(NULL_DIR, "figures/fig_delta_maps_comparison.pdf"),
         fig_delta_maps, width = 200, height = 180, units = "mm", dpi = 300)

  metrics_long <- null_summary |>
    dplyr::select(Comparison, D_future, rho_future,
                  spatial_corr_with_real, corr_delta_vs_deviation) |>
    tidyr::pivot_longer(-Comparison, names_to = "Metric", values_to = "Value") |>
    dplyr::mutate(Metric = dplyr::recode(Metric,
      D_future                = "Schoener's D (future)",
      rho_future              = "Spearman rho (future)",
      spatial_corr_with_real  = "Corr(delta-map, real hybrid delta-map)",
      corr_delta_vs_deviation = "Corr(delta-map, thermoreg. inaccuracy)"))

  fig_metrics <- ggplot(metrics_long, aes(x = Comparison, y = Value, fill = Comparison)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Real_Hybrid"     = "#1F78B4",
                                  "Null-A_Topo"    = "#A6CEE3",
                                  "Null-B_QuadClim" = "#B2DF8A",
                                  "Null-C_RandomRF" = "#FB9A99")) +
    labs(title = "Null model test: divergence metrics",
         y = "Value", x = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text  = element_text(face = "bold"))
  ggsave(file.path(NULL_DIR, "figures/fig_divergence_metrics.pdf"),
         fig_metrics, width = 200, height = 160, units = "mm", dpi = 300)

  save(null_summary, null_survival,
       diff_real, diff_nullA, diff_nullB, diff_nullC,
       pred_nullA_curr_ens, pred_nullA_fut_ens,
       pred_nullB_curr_ens, pred_nullB_fut_ens,
       pred_nullC_curr_ens, pred_nullC_fut_ens,
       file = file.path(NULL_DIR, "models/null_model_results.RData"))
}


# ============================================================================ #
# SECTION C. Controlled comparison (matched-correlative vs hybrid)             #
# ============================================================================ #

if (RUN_C) {
  cat("\n############################################################\n")
  cat("# SECTION C. Controlled comparison\n")
  cat("############################################################\n\n")

  CTRL_DIR <- file.path(OUTDIR, "controlled_comparison")
  invisible(sapply(file.path(CTRL_DIR, c("", "models", "rasters", "figures")),
                   dir.create, showWarnings = FALSE, recursive = TRUE))

  # Reuse rebuilt rasters from Section B if available; otherwise rebuild.
  if (!exists("bio_curr") || !exists("presence_pa")) {
    pd <- rebuild_pipeline_data()
    bio_curr      <- pd$bio_curr
    meandb2       <- pd$meandb2
    meanActivity2 <- pd$meanActivity2
    land_mask     <- pd$land_mask
    presence_pa   <- pd$presence_pa
  }

  vars_correlativo     <- readRDS(file.path(OUTDIR, "models/selected_variables_correlative.rds"))
  vars_hibrido         <- readRDS(file.path(OUTDIR, "models/selected_variables_hybrid.rds"))
  vars_hybrido_climate <- setdiff(vars_hibrido, c("deviation_mean", "Activity"))

  pred_correlativo <- terra::mask(bio_curr[[vars_correlativo]],     land_mask, maskvalues = 0)
  pred_hibrido     <- c(bio_curr[[vars_hybrido_climate]], meandb2, meanActivity2) |>
    terra::mask(land_mask, maskvalues = 0)
  pred_matched     <- terra::mask(bio_curr[[vars_hybrido_climate]], land_mask, maskvalues = 0)

  # Load gbm models from the main pipeline (pure R objects, no terra pointers)
  load(file.path(OUTDIR, "models/03_model_fitting.RData"))
  load(file.path(OUTDIR, "models/04_evaluation.RData"))

  dat_matched <- mk_data(presence_pa, pred_matched)

  rep_matched <- run_repeated_brt(dat_matched, names(pred_matched),
                                  "Matched_Correlative",
                                  presence_pa, pred_correlativo,
                                  n_repeats = N_REPEATS)
  super_models_matched <- flatten_models(rep_matched)
  readr::write_csv(rep_matched$performance_by_fold,
                   file.path(CTRL_DIR, "matched_perf_by_fold.csv"))

  # Performance comparison (3 models)
  perf_corr_pf <- rep_corr$performance_by_fold
  perf_hyb_pf  <- rep_hyb$performance_by_fold
  perf_mat_pf  <- rep_matched$performance_by_fold

  perf_summary <- data.frame(
    Model           = c("Correlative", "Matched_Correlative", "Hybrid"),
    n_predictors    = c(length(vars_correlativo),
                        length(vars_hybrido_climate),
                        length(vars_hibrido)),
    AUC_mean        = c(mean(perf_corr_pf$AUC_test),
                        mean(perf_mat_pf$AUC_test),
                        mean(perf_hyb_pf$AUC_test)),
    AUC_sd          = c(sd(perf_corr_pf$AUC_test),
                        sd(perf_mat_pf$AUC_test),
                        sd(perf_hyb_pf$AUC_test)),
    TSS_mean        = c(mean(perf_corr_pf$TSS_test),
                        mean(perf_mat_pf$TSS_test),
                        mean(perf_hyb_pf$TSS_test)),
    TSS_sd          = c(sd(perf_corr_pf$TSS_test),
                        sd(perf_mat_pf$TSS_test),
                        sd(perf_hyb_pf$TSS_test)),
    has_mechanistic = c(FALSE, FALSE, TRUE)
  )

  paired_mh <- dplyr::inner_join(
    perf_mat_pf |> dplyr::select(Repeat, fold, AUC_matched = AUC_test, TSS_matched = TSS_test),
    perf_hyb_pf |> dplyr::select(Repeat, fold, AUC_hybrid  = AUC_test, TSS_hybrid  = TSS_test),
    by = c("Repeat", "fold"))
  paired_mc <- dplyr::inner_join(
    perf_mat_pf  |> dplyr::select(Repeat, fold, AUC_matched = AUC_test, TSS_matched = TSS_test),
    perf_corr_pf |> dplyr::select(Repeat, fold, AUC_corr    = AUC_test, TSS_corr    = TSS_test),
    by = c("Repeat", "fold"))

  wilcox_mh_auc <- wilcox.test(paired_mh$AUC_hybrid, paired_mh$AUC_matched, paired = TRUE)
  wilcox_mh_tss <- wilcox.test(paired_mh$TSS_hybrid, paired_mh$TSS_matched, paired = TRUE)
  wilcox_mc_auc <- wilcox.test(paired_mc$AUC_matched, paired_mc$AUC_corr,    paired = TRUE)
  wilcox_mc_tss <- wilcox.test(paired_mc$TSS_matched, paired_mc$TSS_corr,    paired = TRUE)

  wilcox_summary <- data.frame(
    Comparison = c("Hybrid_vs_Matched", "Hybrid_vs_Matched",
                   "Matched_vs_Correlative", "Matched_vs_Correlative"),
    Metric  = c("AUC", "TSS", "AUC", "TSS"),
    W       = c(wilcox_mh_auc$statistic, wilcox_mh_tss$statistic,
                wilcox_mc_auc$statistic, wilcox_mc_tss$statistic),
    p_value = c(wilcox_mh_auc$p.value, wilcox_mh_tss$p.value,
                wilcox_mc_auc$p.value, wilcox_mc_tss$p.value),
    median_delta = c(median(paired_mh$AUC_hybrid - paired_mh$AUC_matched),
                     median(paired_mh$TSS_hybrid - paired_mh$TSS_matched),
                     median(paired_mc$AUC_matched - paired_mc$AUC_corr),
                     median(paired_mc$TSS_matched - paired_mc$TSS_corr)))

  readr::write_csv(perf_summary,   file.path(CTRL_DIR, "performance_comparison_3models.csv"))
  readr::write_csv(wilcox_summary, file.path(CTRL_DIR, "wilcoxon_tests.csv"))
  print(perf_summary,   digits = 4)
  print(wilcox_summary, digits = 4)

  # Project all three models, current + future
  pred_corr_curr_ctl <- predict_ensemble(super_models_corr,    pred_correlativo)
  pred_hyb_curr_ctl  <- predict_ensemble(super_models_hyb,     pred_hibrido)
  pred_matched_curr  <- predict_ensemble(super_models_matched, pred_matched)
  for (obj in c("pred_corr_curr_ctl", "pred_hyb_curr_ctl", "pred_matched_curr")) {
    o <- get(obj); o$r_mean <- terra::mask(o$r_mean, land_mask); o$r_sd <- terra::mask(o$r_sd, land_mask)
    assign(obj, o)
  }

  bio_fut <- geodata::cmip6_world(model = "MIROC6", ssp = "585",
                                  time = "2041-2060", var = "bioc",
                                  res = 2.5, path = "data", download = FALSE)
  names(bio_fut) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
  bio_fut <- terra::resample(bio_fut, meandb2) |> terra::mask(meandb2)

  meandb2_fut       <- terra::rast(MECH_FUT_WARM_MEANDB);   names(meandb2_fut)       <- "deviation_mean"
  meanActivity2_fut <- terra::rast(MECH_FUT_WARM_ACTIVITY); names(meanActivity2_fut) <- "Activity"

  pred_correlativo_fut   <- bio_fut[[vars_correlativo]]      |>
    terra::mask(land_mask, maskvalues = 0)
  pred_matched_fut_stack <- bio_fut[[vars_hybrido_climate]]  |>
    terra::mask(land_mask, maskvalues = 0)
  pred_hibrido_fut       <- c(bio_fut[[vars_hybrido_climate]], meandb2_fut, meanActivity2_fut) |>
    terra::mask(land_mask, maskvalues = 0)

  pred_corr_fut_ctl <- predict_ensemble(super_models_corr,    pred_correlativo_fut)
  pred_hyb_fut_ctl  <- predict_ensemble(super_models_hyb,     pred_hibrido_fut)
  pred_matched_fut  <- predict_ensemble(super_models_matched, pred_matched_fut_stack)
  for (obj in c("pred_corr_fut_ctl", "pred_hyb_fut_ctl", "pred_matched_fut")) {
    o <- get(obj); o$r_mean <- terra::mask(o$r_mean, land_mask); o$r_sd <- terra::mask(o$r_sd, land_mask)
    assign(obj, o)
  }

  terra::writeRaster(pred_matched_curr$r_mean,
                     file.path(CTRL_DIR, "rasters/current_matched_mean.tif"),
                     overwrite = TRUE, datatype = "FLT4S")
  terra::writeRaster(pred_matched_fut$r_mean,
                     file.path(CTRL_DIR, "rasters/future_matched_mean.tif"),
                     overwrite = TRUE, datatype = "FLT4S")

  # --- Niche overlap (the key controlled test) --------------------------------
  D_corr_hyb_curr   <- schoeners_d(pred_corr_curr_ctl$r_mean, pred_hyb_curr_ctl$r_mean)
  D_match_hyb_curr  <- schoeners_d(pred_matched_curr$r_mean,  pred_hyb_curr_ctl$r_mean)
  D_corr_match_curr <- schoeners_d(pred_corr_curr_ctl$r_mean, pred_matched_curr$r_mean)
  rho_corr_hyb_curr  <- spearman_r(pred_corr_curr_ctl$r_mean, pred_hyb_curr_ctl$r_mean)
  rho_match_hyb_curr <- spearman_r(pred_matched_curr$r_mean,  pred_hyb_curr_ctl$r_mean)

  D_corr_hyb_fut    <- schoeners_d(pred_corr_fut_ctl$r_mean, pred_hyb_fut_ctl$r_mean)
  D_match_hyb_fut   <- schoeners_d(pred_matched_fut$r_mean,  pred_hyb_fut_ctl$r_mean)
  rho_corr_hyb_fut  <- spearman_r(pred_corr_fut_ctl$r_mean, pred_hyb_fut_ctl$r_mean)
  rho_match_hyb_fut <- spearman_r(pred_matched_fut$r_mean,  pred_hyb_fut_ctl$r_mean)

  delta_corr  <- pred_corr_fut_ctl$r_mean - pred_corr_curr_ctl$r_mean
  delta_match <- pred_matched_fut$r_mean  - pred_matched_curr$r_mean
  delta_hyb   <- pred_hyb_fut_ctl$r_mean  - pred_hyb_curr_ctl$r_mean

  vals_ok <- complete.cases(terra::values(delta_corr,  na.rm = FALSE),
                            terra::values(delta_match, na.rm = FALSE),
                            terra::values(delta_hyb,   na.rm = FALSE))
  corr_delta_corr_hyb  <- cor(terra::values(delta_corr)[vals_ok],
                              terra::values(delta_hyb)[vals_ok], method = "spearman")
  corr_delta_match_hyb <- cor(terra::values(delta_match)[vals_ok],
                              terra::values(delta_hyb)[vals_ok], method = "spearman")

  overlap_results <- data.frame(
    Comparison         = c("Correlative_vs_Hybrid", "Matched_vs_Hybrid", "Correlative_vs_Matched"),
    D_current          = c(D_corr_hyb_curr,  D_match_hyb_curr, D_corr_match_curr),
    D_future           = c(D_corr_hyb_fut,   D_match_hyb_fut,  NA),
    rho_current        = c(rho_corr_hyb_curr, rho_match_hyb_curr, NA),
    rho_future         = c(rho_corr_hyb_fut,  rho_match_hyb_fut,  NA),
    delta_spatial_corr = c(corr_delta_corr_hyb, corr_delta_match_hyb, NA))
  readr::write_csv(overlap_results, file.path(CTRL_DIR, "overlap_results.csv"))
  print(overlap_results, digits = 4)

  # --- Range change comparison ------------------------------------------------
  thresholds <- readRDS(file.path(OUTDIR, "results_summary/thresholds_maxSSS_ensemble.rds"))
  thr_corr   <- thresholds$correlativo_mean_ensemble
  thr_hyb    <- thresholds$hibrido_mean_ensemble

  # maxSSS threshold for the matched-correlative ensemble
  pres_pts <- presence_pa[presence_pa$pr_ab == 1, ]
  pres_sf  <- sf::st_as_sf(pres_pts,
                           coords = c("decimalLongitude", "decimalLatitude"),
                           crs    = terra::crs(pred_matched))
  pres_locs  <- terra::extract(pred_matched_curr$r_mean, pres_sf)[, 2]
  ens_vals   <- terra::values(pred_matched_curr$r_mean, na.rm = TRUE)
  bg_vals    <- sample(ens_vals, min(5000, length(ens_vals)))
  obs_combo  <- c(rep(1, length(pres_locs)), rep(0, length(bg_vals)))
  pred_combo <- c(pres_locs, bg_vals)
  ok <- complete.cases(obs_combo, pred_combo)
  thr_seq    <- seq(min(pred_combo[ok]), max(pred_combo[ok]), length.out = 200)
  sss_vec    <- vapply(thr_seq, function(thr) {
    yhat <- as.numeric(pred_combo[ok] >= thr)
    sens <- sum(yhat == 1 & obs_combo[ok] == 1) / max(sum(obs_combo[ok] == 1), 1)
    spec <- sum(yhat == 0 & obs_combo[ok] == 0) / max(sum(obs_combo[ok] == 0), 1)
    sens + spec - 1
  }, numeric(1))
  thr_matched_maxSSS <- thr_seq[which.max(sss_vec)]

  cell_area_km2 <- prod(terra::res(pred_corr_curr_ctl$r_mean)) * 111.32 ^ 2 *
                   cos(mean(terra::ext(pred_corr_curr_ctl$r_mean)[3:4]) * pi / 180)

  compute_range_change <- function(r_curr, r_fut, thr, label) {
    bin_curr <- terra::ifel(is.na(r_curr), NA, r_curr >= thr)
    bin_fut  <- terra::ifel(is.na(r_fut),  NA, r_fut  >= thr)
    n_curr <- sum(terra::values(bin_curr, na.rm = TRUE) == 1)
    n_fut  <- sum(terra::values(bin_fut,  na.rm = TRUE) == 1)
    area_curr  <- n_curr * cell_area_km2
    area_fut   <- n_fut  * cell_area_km2
    pct_change <- (area_fut - area_curr) / area_curr * 100
    data.frame(Model            = label,
               Threshold        = "maxSSS",
               Area_current_km2 = area_curr,
               Area_future_km2  = area_fut,
               Change_pct       = pct_change)
  }

  range_change <- rbind(
    compute_range_change(pred_corr_curr_ctl$r_mean, pred_corr_fut_ctl$r_mean,
                          thr_corr,           "Correlative"),
    compute_range_change(pred_matched_curr$r_mean,  pred_matched_fut$r_mean,
                          thr_matched_maxSSS, "Matched_Correlative"),
    compute_range_change(pred_hyb_curr_ctl$r_mean,  pred_hyb_fut_ctl$r_mean,
                          thr_hyb,            "Hybrid"))
  readr::write_csv(range_change, file.path(CTRL_DIR, "range_change_comparison.csv"))
  print(range_change, digits = 4)

  save(rep_matched, super_models_matched, pred_matched,
       pred_matched_curr, pred_matched_fut,
       overlap_results, range_change, wilcox_summary, perf_summary,
       thr_matched_maxSSS, delta_corr, delta_match, delta_hyb,
       file = file.path(CTRL_DIR, "models/controlled_comparison_results.RData"))
}


# ============================================================================ #
# Reproducibility footer                                                       #
# ============================================================================ #

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "session_info_integrated.txt"))

cat("\n  OK: RCode_SDM_integrated.R complete\n")
