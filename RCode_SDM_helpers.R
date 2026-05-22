################################################################################
## RCode_SDM_helpers.R
##
## Shared library for the species distribution modelling (SDM) pipeline of:
##
##   Rubalcaba, Fandos & Diaz. "Behavioural thermoregulation and microclimate
##   reshape climate-driven range forecasts."
##
## All RCode_SDM_*.R analysis scripts source this file at startup. It provides:
##   1. Package management (auto-install + library loading)
##   2. Project-wide constants and seeds (BRT hyperparameters, paths, CRS)
##   3. Plotting palettes and ggplot themes
##   4. Statistical SDM utilities (select07, predictSDM, evalSDM, crossvalSDM,
##      TSS, expl_deviance, partial-response curves)
##   5. Raster utilities (mask alignment, MESS with driver fallback, niche
##      overlap metrics, class-balanced weights, categorical/MESS map plots)
##   6. Loader for GCM/SSP-specific future mechanistic layers
##
## Authorship of the statistical helpers (select07, evalSDM, crossvalSDM,
## predictSDM, TSS, expl_deviance, partial_response, inflated_response):
## adapted from the "mecofun" teaching package (Zurell et al., MEcoFun working
## group), with minor edits for compatibility with this pipeline.
##
## Required packages: terra, sf, dplyr, tidyr, ggplot2, patchwork, blockCV,
##   dismo, gbm, ROCR, PresenceAbsence, readr, rasterVis, purrr, spThin, pdp,
##   geodata, pROC, rnaturalearth, rnaturalearthdata, predicts, mgcv, glm2,
##   lhs, maxnet, randomForest, tibble.
################################################################################

# ---------------------------------------------------------------------------- #
# 0. Package management
# ---------------------------------------------------------------------------- #

required_packages <- c(
  "terra", "sf", "dplyr", "tidyr", "ggplot2", "patchwork",
  "blockCV", "dismo", "gbm", "ROCR", "PresenceAbsence",
  "readr", "rasterVis", "purrr", "spThin", "pdp",
  "geodata", "pROC", "rnaturalearth", "rnaturalearthdata",
  "predicts", "mgcv"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))


# ---------------------------------------------------------------------------- #
# 1. Project-wide constants and seeds
# ---------------------------------------------------------------------------- #

MASTER_SEED       <- 123
set.seed(MASTER_SEED)

# Default output directory; override before sourcing if needed.
if (!exists("OUTDIR")) OUTDIR <- "output"

# Sampling and partitioning
N_BACKGROUND      <- 5000
K_FOLDS           <- 5
THIN_DIST_KM      <- 3

# BRT hyperparameters (fixed; no grid search)
N_REPEATS         <- 10
SHRINKAGE         <- 0.001
INTERACTION_DEPTH <- 5
BAG_FRACTION      <- 0.75
N_TREES_MAX       <- 10000

# Background sampling (bias-corrected, tempered_uniform)
BG_DOMAIN         <- "M_core"
BG_SAMPLING_MODE  <- "tempered_uniform"
GAMMA             <- 0.6
LAMBDA            <- 0.30
OUTER_QUOTA       <- 0.35

# Buffer distances and equal-area CRS
BUF_CORE_KM       <- 50
BUF_BG_ADD_KM     <- 100
CRS_EQUAL_AREA    <- "EPSG:6933"

# Initialise output subdirectories
.outdir_subdirs <- c("", "rasters", "models", "figures",
                     "threshold_sensitivity", "results_summary",
                     "performance_metrics")
invisible(sapply(file.path(OUTDIR, .outdir_subdirs),
                 dir.create, showWarnings = FALSE, recursive = TRUE))


# ---------------------------------------------------------------------------- #
# 2. Plotting palettes and themes
# ---------------------------------------------------------------------------- #

pub_colors <- c("1" = "#A6611A", "2" = "#1F78B4", "3" = "#B2DF8A")
col_diff   <- c(low = "#D73027", mid = "white", high = "#4575B4")

theme_set(theme_void(base_size = 10))
theme_minimal_nc <- function() {
  theme_void(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title    = element_text(size = 9, face = "bold"),
          legend.text     = element_text(size = 8),
          plot.margin     = margin(5, 5, 5, 5))
}

# ---------------------------------------------------------------------------- #
# 3. Variable selection (Dormann et al. 2013 select07)
# ---------------------------------------------------------------------------- #

#' Pairwise correlation-based predictor selection (select07)
#'
#' For each pair of predictors with abs(rho) >= threshold, drops the one with
#' lower univariate explanatory power (AIC of glm/gam against the response).
#' See Dormann et al. (2013, Ecography) for rationale.
select07 <- function(X, y,
                     family    = "binomial",
                     univar    = "glm2",
                     threshold = 0.7,
                     method    = "spearman",
                     sequence  = NULL,
                     weights   = NULL) {

  var.imp <- function(variable, response, univar, family, weights) {
    m1 <- switch(univar,
                 glm1 = glm(response ~ variable, family = family, weights = weights),
                 glm2 = glm(response ~ poly(variable, 2), family = family, weights = weights),
                 gam  = mgcv::gam(response ~ s(variable, k = 4), family = family, weights = weights))
    AIC(m1)
  }

  cm <- cor(X, method = method)

  if (is.null(sequence)) {
    a <- try(var.imp(X[, 1], y, univar = univar, family = family, weights = weights))
    if (!is.numeric(a)) stop("invalid univar method")
    imp <- apply(X, 2, var.imp, response = y, family = family,
                 univar = univar, weights = weights)
    sort.imp <- names(sort(imp))
  } else {
    sort.imp <- sequence
  }

  pairs <- which(abs(cm) >= threshold, arr.ind = TRUE)
  index <- which(pairs[, 1] == pairs[, 2])
  pairs <- pairs[-index, ]

  exclude <- NULL
  for (i in seq_along(sort.imp)) {
    if ((sort.imp[i] %in% rownames(pairs)) && !(sort.imp[i] %in% exclude)) {
      cv <- cm[setdiff(rownames(cm), exclude), sort.imp[i]]
      cv <- cv[setdiff(names(cv), sort.imp[1:i])]
      exclude <- c(exclude, names(which(abs(cv) >= threshold)))
    }
  }

  pred_sel <- sort.imp[!(sort.imp %in% unique(exclude)), drop = FALSE]
  list(AIC = sort(imp), cor_mat = cm, pred_sel = pred_sel)
}


# ---------------------------------------------------------------------------- #
# 4. Skill metrics
# ---------------------------------------------------------------------------- #

#' True Skill Statistic from a confusion matrix
TSS <- function(cmx) {
  PresenceAbsence::sensitivity(cmx, st.dev = FALSE) +
    PresenceAbsence::specificity(cmx, st.dev = FALSE) - 1
}

#' Explained deviance for a binary SDM
expl_deviance <- function(obs, pred,
                          family  = "binomial",
                          weights = rep(1, length(obs))) {
  if (family == "binomial") {
    pred <- pmin(pmax(pred, 1e-05), 0.9999)
  }
  null_pred <- rep(mean(obs), length(obs))
  1 - (dismo::calc.deviance(obs, pred,      family = family, weights = weights) /
       dismo::calc.deviance(obs, null_pred, family = family, weights = weights))
}


# ---------------------------------------------------------------------------- #
# 5. Prediction and evaluation
# ---------------------------------------------------------------------------- #

#' Unified prediction wrapper for heterogeneous SDM model classes
#'
#' Returns probability/suitability on the response scale for each supported
#' model class (Bioclim, Domain, glm/Gam/gam, negbin, rpart, randomForest,
#' gbm, maxnet).
predictSDM <- function(model, newdata) {
  switch(class(model)[1],
         Bioclim       = predict(model, newdata),
         Domain        = predict(model, newdata),
         glm           = predict(model, newdata, type = "response"),
         Gam           = predict(model, newdata, type = "response"),
         gam           = predict(model, newdata, type = "response"),
         negbin        = predict(model, newdata, type = "response"),
         rpart         = predict(model, newdata),
         randomForest.formula = switch(model$type,
                                       regression     = predict(model, newdata, type = "response"),
                                       classification = predict(model, newdata, type = "prob")[, 2]),
         randomForest         = switch(model$type,
                                       regression     = predict(model, newdata, type = "response"),
                                       classification = predict(model, newdata, type = "prob")[, 2]),
         gbm           = switch(ifelse(is.null(model$gbm.call), "GBM", "GBM.STEP"),
                                GBM.STEP = predict.gbm(model, newdata,
                                                       n.trees = model$gbm.call$best.trees,
                                                       type    = "response"),
                                GBM      = predict.gbm(model, newdata,
                                                       n.trees = model$n.trees,
                                                       type    = "response")),
         maxnet        = predict(model, newdata, type = "logistic"))
}

#' Performance metrics on a single hold-out
#'
#' Returns AUC, TSS, Kappa, Sens, Spec, PCC and D2 for a vector of predictions
#' against observed presence/absence. The threshold (for TSS/Kappa/Sens/Spec/
#' PCC) is either supplied or selected via PresenceAbsence::optimal.thresholds
#' using `thresh.method` (default "MaxSens+Spec", equivalent to maxTSS / maxSSS
#' for the binary case used here).
evalSDM <- function(observation, predictions,
                    thresh        = NULL,
                    thresh.method = "MaxSens+Spec",
                    req.sens      = 0.85,
                    req.spec      = 0.85,
                    FPC           = 1,
                    FNC           = 1,
                    weigths       = rep(1, length(observation))) {

  thresh.dat <- data.frame(ID  = seq_along(observation),
                           obs = observation,
                           pred = predictions)

  if (is.null(thresh)) {
    thresh.mat <- PresenceAbsence::optimal.thresholds(
      DATA = thresh.dat,
      req.sens = req.sens, req.spec = req.spec,
      FPC = FPC, FNC = FNC)
    thresh <- thresh.mat[thresh.mat$Method == thresh.method, 2]
  }

  cmx.opt <- PresenceAbsence::cmx(DATA = thresh.dat, threshold = thresh)

  data.frame(
    AUC    = PresenceAbsence::auc(thresh.dat, st.dev = FALSE),
    TSS    = TSS(cmx.opt),
    Kappa  = PresenceAbsence::Kappa(cmx.opt, st.dev = FALSE),
    Sens   = PresenceAbsence::sensitivity(cmx.opt, st.dev = FALSE),
    Spec   = PresenceAbsence::specificity(cmx.opt, st.dev = FALSE),
    PCC    = PresenceAbsence::pcc(cmx.opt, st.dev = FALSE),
    D2     = expl_deviance(observation, predictions, weights = weigths),
    thresh = thresh)
}


# ---------------------------------------------------------------------------- #
# 6. Cross-validation
# ---------------------------------------------------------------------------- #

#' k-fold cross-validation with model-class dispatch
#'
#' Refits the supplied model on each training fold (using the appropriate
#' constructor for the model class) and returns a vector of out-of-fold
#' predictions aligned with `traindat`. Folds can be supplied via `kfold`
#' as either an integer (random stratified by species) or a numeric vector
#' of pre-assigned fold IDs (used for spatial blocking).
crossvalSDM <- function(model,
                        kfold           = 5,
                        traindat,
                        colname_species,
                        colname_pred,
                        env_r           = NULL,
                        colname_coord   = NULL,
                        weights         = NULL) {

  weights.full <- weights

  if (length(kfold) == 1) {
    ks <- dismo::kfold(traindat, k = kfold, by = traindat[, colname_species])
  } else {
    ks <- kfold
    kfold <- length(unique(kfold))
  }

  cross_val_preds <- numeric(length = nrow(traindat))

  for (i in seq_len(kfold)) {
    cv_train <- traindat[ks != i, ]
    cv_test  <- traindat[ks == i, ]

    if (class(model)[1] == "gbm") {
      cv_train_gbm <- cv_train
      if (!is.null(model$gbm.call)) {
        names(cv_train_gbm)[names(cv_train_gbm) == colname_species] <- model$response.name
      }
    }
    if (!is.null(weights)) weights <- weights.full[ks != i]

    modtmp <- switch(class(model)[1],
      Bioclim = dismo::bioclim(env_r[[colname_pred]],
                               cv_train[cv_train[, colname_species] == 1, colname_coord]),
      Domain  = dismo::domain(env_r[[colname_pred]],
                              cv_train[cv_train[, colname_species] == 1, colname_coord]),
      glm     = update(model, data = cv_train),
      Gam     = update(model, data = cv_train, weights = weights),
      gam     = update(model, data = cv_train, weights = weights),
      negbin  = update(model, data = cv_train),
      rpart   = update(model, data = cv_train),
      randomForest         = update(model, data = cv_train),
      randomForest.formula = update(model, data = cv_train),
      gbm     = switch(ifelse(is.null(model$gbm.call), "GBM", "GBM.STEP"),
        GBM.STEP = gbm::gbm(model$call, "bernoulli",
                            data              = cv_train_gbm[, c(colname_pred, model$response.name)],
                            n.trees           = model$gbm.call$best.trees,
                            shrinkage         = model$gbm.call$learning.rate,
                            bag.fraction      = model$gbm.call$bag.fraction,
                            interaction.depth = model$gbm.call$tree.complexity),
        GBM      = gbm::gbm(model$call, "bernoulli",
                            data              = cv_train_gbm[, c(colname_pred, model$response.name)],
                            n.trees           = model$n.trees,
                            shrinkage         = model$shrinkage,
                            bag.fraction      = model$bag.fraction,
                            interaction.depth = model$interaction.depth)),
      maxnet  = maxnet::maxnet(p = cv_train[, colname_species],
                               data = cv_train[, colname_pred, drop = FALSE]))

    cross_val_preds[ks == i] <- predictSDM(modtmp, cv_test[, colname_pred, drop = FALSE])
  }

  cross_val_preds
}


# ---------------------------------------------------------------------------- #
# 7. Partial-response curves
# ---------------------------------------------------------------------------- #

#' Partial response curve (predictor i held free, others at column means)
partial_response <- function(object, predictors,
                             select.columns = NULL,
                             label          = NULL,
                             len            = 50,
                             col            = "black",
                             ylab           = NULL,
                             ...) {
  inflated_response(object, predictors, select.columns, label, len,
                    method = "mean", col.curves = col, ylab = ylab, ...)
}

#' Inflated response curves with optional Latin Hypercube background
#'
#' Plots response curves over predictor `i` while sampling the remaining
#' predictors using one of three strategies: `stat3` (min/median/max),
#' `stat6` (full quantiles), or `mean` (mean of each remaining predictor).
#' When the number of combinations exceeds `lhsample`, a Latin Hypercube is
#' drawn instead.
inflated_response <- function(object, predictors,
                              select.columns = NULL,
                              label          = NULL,
                              len            = 50,
                              lhsample       = 100,
                              lwd            = 1,
                              ylab           = NULL,
                              method         = "stat3",
                              disp           = "all",
                              overlay.mean   = TRUE,
                              col.curves     = "grey",
                              col.novel      = "grey",
                              col.mean       = "black",
                              lwd.known      = 2,
                              lwd.mean       = 2,
                              ylim           = c(0, 1),
                              ...) {
  if (is.null(select.columns)) select.columns <- seq_len(ncol(predictors))

  for (i in select.columns) {
    summaries <- data.frame(matrix(0, 6, ncol(predictors)))
    for (iz in 1:ncol(predictors)) summaries[, iz] <- summary(predictors[, iz])

    if (method == "stat3") {
      summaries.j <- as.matrix(summaries[c(1, 4, 6), -i], ncol = (ncol(predictors) - 1))
      comb <- min(lhsample, 3 ^ (ncol(predictors) - 1)); nc <- 3
    } else if (method == "stat6") {
      summaries.j <- as.matrix(summaries[, -i],          ncol = (ncol(predictors) - 1))
      comb <- min(lhsample, 6 ^ (ncol(predictors) - 1)); nc <- 6
    } else if (method == "mean") {
      summaries.j <- as.matrix(summaries[4, -i],         ncol = (ncol(predictors) - 1))
      comb <- 1; nc <- 1; overlay.mean <- FALSE
    }

    dummy.j <- as.matrix(predictors[1:len, -i], ncol = (ncol(predictors) - 1))

    if (comb < lhsample) {
      mat <- vector("list", ncol(dummy.j))
      for (m in 1:ncol(dummy.j)) mat[[m]] <- 1:nc
      mat <- expand.grid(mat)
    } else {
      mat <- round(qunif(lhs::randomLHS(lhsample, ncol(dummy.j)),
                         1, nrow(summaries.j)), 0)
    }

    if (is.null(label)) label <- names(predictors)

    for (r in 1:nrow(mat)) {
      for (j in 1:ncol(dummy.j)) {
        dummy.j[, j] <- as.vector(rep(summaries.j[mat[r, j], j], len))
      }
      dummy <- data.frame(seq(min(predictors[, i]), max(predictors[, i]), length = len), dummy.j)
      names(dummy)[-1] <- names(predictors)[-i]
      names(dummy)[1]  <- names(predictors)[i]
      curves <- predictSDM(object, dummy)

      if (disp == "all") {
        if (r == 1) {
          if (i == 1)
            plot(dummy[, names(predictors)[i]], curves, type = "l", ylim = ylim,
                 xlab = label[i], lwd = lwd, col = col.curves, ylab = ylab, ...)
          else
            plot(dummy[, names(predictors)[i]], curves, type = "l", ylim = ylim,
                 xlab = label[i], lwd = lwd, col = col.curves, ylab = "", ...)
        } else {
          lines(dummy[, names(predictors)[i]], curves, lwd = lwd, col = col.curves, ...)
        }
      }

      if (disp == "eo.mask") {
        novel <- eo.mask(predictors, dummy)
        curves.known <- curves; curves.known[novel == 1] <- NA
        curves.novel <- curves; curves.novel[novel == 0] <- NA

        if (r == 1) {
          if (i == 1) {
            plot(dummy[, names(predictors)[i]], curves.known, type = "l", ylim = ylim,
                 xlab = label[i], lwd = lwd.known, col = col.curves, ylab = ylab, ...)
            lines(dummy[, names(predictors)[i]], curves.novel,
                  lwd = lwd, col = col.novel, lty = "dotted", ...)
          } else {
            plot(dummy[, names(predictors)[i]], curves.known, type = "l", ylim = ylim,
                 xlab = label[i], lwd = lwd.known, col = col.curves, ylab = "", ...)
            lines(dummy[, names(predictors)[i]], curves.novel,
                  lwd = lwd, col = col.novel, lty = "dotted", ...)
          }
        } else {
          lines(dummy[, names(predictors)[i]], curves.known, lwd = lwd.known, col = col.curves, ...)
          lines(dummy[, names(predictors)[i]], curves.novel,
                lwd = lwd, col = col.novel, lty = "dotted", ...)
        }
      }
    }

    if (isTRUE(overlay.mean)) {
      dummy <- predictors[1:len, ]
      dummy[, i] <- seq(min(predictors[, i]), max(predictors[, i]), length = len)
      for (j in 1:ncol(predictors)) {
        if (j != i) dummy[, j] <- rep(mean(predictors[, j]), len)
      }
      curves <- predictSDM(object, dummy)
      lines(dummy[, names(predictors)[i]], curves, lwd = lwd.mean, col = col.mean, ...)
    }
  }
}


# ---------------------------------------------------------------------------- #
# 8. Raster and map utilities
# ---------------------------------------------------------------------------- #

#' Resample/project a mask raster onto a target raster
#'
#' Tries direct resample first; falls back to project-then-resample if the
#' source CRS differs from the target.
align_mask_to <- function(mask_src, target) {
  m <- try(terra::resample(mask_src, target, method = "near"), silent = TRUE)
  if (inherits(m, "try-error")) {
    m <- terra::project(mask_src, terra::crs(target), method = "near")
    m <- terra::resample(m, target, method = "near")
  }
  terra::crop(m, terra::ext(target))
}

#' 4-class agreement map from two binary (correlative / hybrid) rasters
#'
#' Returns a SpatRaster with values:
#'   1 = correlative only, 2 = hybrid only, 3 = both, NA = neither.
cat_map_from <- function(r_bin_corr, r_bin_hyb) {
  both      <- r_bin_corr & r_bin_hyb
  only_corr <- r_bin_corr & !r_bin_hyb
  only_hyb  <- !r_bin_corr & r_bin_hyb
  absent    <- !r_bin_corr & !r_bin_hyb
  out <- absent * 0 + only_corr * 1 + only_hyb * 2 + both * 3
  out[out == 0] <- NA
  out
}

#' MESS computation with driver/package fallback chain
#'
#' Tries `predicts::mess` (terra-native), then `predicts::mess` on a raster
#' stack, then `dismo::mess`. Returns a SpatRaster.
mess_robusto <- function(proj_r, ref_df) {
  ref_df <- stats::na.omit(ref_df)
  ok1 <- try(predicts::mess(proj_r, ref_df), silent = TRUE)
  if (!inherits(ok1, "try-error")) return(ok1)
  proj_r_raster <- raster::stack(raster::raster(proj_r[[1]]))
  if (terra::nlyr(proj_r) > 1) {
    for (i in 2:terra::nlyr(proj_r))
      proj_r_raster <- raster::addLayer(proj_r_raster, raster::raster(proj_r[[i]]))
  }
  ok2 <- try(predicts::mess(proj_r_raster, ref_df), silent = TRUE)
  if (!inherits(ok2, "try-error")) return(terra::rast(ok2))
  ok3 <- try(dismo::mess(proj_r_raster, ref_df), silent = TRUE)
  if (!inherits(ok3, "try-error")) return(terra::rast(ok3))
  stop("Could not compute MESS with predicts or dismo.")
}

#' Class-balanced observation weights for binary y
make_class_weights <- function(y) {
  p <- mean(y == 1)
  w <- ifelse(y == 1, 0.5 / p, 0.5 / (1 - p))
  as.numeric(w)
}

#' Schoener's D between two suitability rasters
schoeners_d <- function(r1, r2) {
  v1 <- terra::values(r1, mat = FALSE)
  v2 <- terra::values(r2, mat = FALSE)
  ok <- is.finite(v1) & is.finite(v2) & (v1 >= 0) & (v2 >= 0)
  p1 <- v1[ok] / sum(v1[ok])
  p2 <- v2[ok] / sum(v2[ok])
  1 - 0.5 * sum(abs(p1 - p2))
}

#' Pearson correlation between two suitability rasters
pearson_r <- function(r1, r2) {
  v1 <- terra::values(r1, mat = FALSE)
  v2 <- terra::values(r2, mat = FALSE)
  ok <- is.finite(v1) & is.finite(v2)
  stats::cor(v1[ok], v2[ok], method = "pearson")
}

#' Categorical agreement map (ggplot)
plot_categorical <- function(r_cat, land_mask, title = NULL,
                             ocean_col   = "#B3DDF2",
                             land_na_col = "#EDEDED") {
  land_mask_aligned <- align_mask_to(land_mask, r_cat)
  ex <- terra::ext(r_cat)
  df_cat  <- as.data.frame(r_cat,             xy = TRUE, na.rm = FALSE)
  df_land <- as.data.frame(land_mask_aligned, xy = TRUE, na.rm = FALSE)
  stopifnot(nrow(df_cat) == nrow(df_land))
  names(df_cat)[3]  <- "value"
  names(df_land)[3] <- "land"
  df <- dplyr::bind_cols(df_cat, land = df_land$land)
  is_land <- !is.na(df$land)

  ggplot() +
    geom_raster(data = df[!is_land, ],                  aes(x = x, y = y), fill = ocean_col) +
    geom_raster(data = df[is_land & is.na(df$value), ], aes(x = x, y = y), fill = land_na_col) +
    geom_raster(data = df[!is.na(df$value), ],          aes(x = x, y = y, fill = factor(value))) +
    scale_fill_manual(values = c("1" = "#E69F00", "2" = "#009E73", "3" = "#756BB1"),
                      labels = c("1" = "Correlative only", "2" = "Hybrid only", "3" = "Both"),
                      drop = FALSE, name = "Model agreement") +
    coord_equal(xlim = c(ex[1], ex[2]), ylim = c(ex[3], ex[4]), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid       = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
          axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          legend.position  = "bottom")
}

#' MESS surface map (ggplot)
plot_mess <- function(r_mess, title = NULL) {
  if (inherits(r_mess, "SpatRaster")) {
    df <- terra::as.data.frame(r_mess, xy = TRUE, na.rm = FALSE)
    names(df) <- c("x", "y", "value")
  } else {
    df <- rasterVis::gplot(r_mess)$data
  }
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = "#67001f", mid = "white", high = "#053061",
                         midpoint = 0, na.value = "grey92", name = "MESS") +
    coord_equal(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid       = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title       = element_text(face = "bold", size = 10, hjust = 0.5))
}


# ---------------------------------------------------------------------------- #
# 9. Future mechanistic layer loader
# ---------------------------------------------------------------------------- #

# Default location for GCM/SSP-specific future mechanistic layers
# (deviation_mean, Activity), generated by the biophysical model in
# RCode_biophysical_model.R. These layers are NOT shipped in the GitHub
# repository: they will be archived on Zenodo together with the rest of
# intermediate model outputs upon acceptance, or can be regenerated from
# the biophysical model script. Place the layers in the directory below
# (or override MECH_FUT_DIR before sourcing this script).
#
# File naming convention:
#   meandb_<GCM>_ssp<SSP>_<period>.grd
#   meanActivity_<GCM>_ssp<SSP>_<period>.grd
#
# Availability (as of 2026-04-16):
#   ssp585 : MIROC6, MPI-ESM1-2-HR, CNRM-CM6-1, EC-Earth3-Veg, UKESM1-0-LL
#   ssp245 : MIROC6
if (!exists("MECH_FUT_DIR")) {
  MECH_FUT_DIR <- file.path("Sources", "SDM", "mechanistic_future")
}

#' Load GCM/SSP-specific future mechanistic layers
#'
#' @param gcm    Character. GCM identifier (e.g. "MIROC6", "UKESM1-0-LL").
#' @param ssp    Character. SSP code: "ssp585" or "ssp245".
#' @param period Character. Time window suffix. Default "2041-2060".
#' @return List with two SpatRasters renamed to match the predictors expected
#'   by the trained ensembles ("deviation_mean" and "Activity").
load_mech_future <- function(gcm, ssp, period = "2041-2060") {
  stopifnot(is.character(gcm), length(gcm) == 1)
  stopifnot(ssp %in% c("ssp585", "ssp245"))

  meandb_path  <- file.path(MECH_FUT_DIR,
                            sprintf("meandb_%s_%s_%s.grd",       gcm, ssp, period))
  meanact_path <- file.path(MECH_FUT_DIR,
                            sprintf("meanActivity_%s_%s_%s.grd", gcm, ssp, period))

  missing <- c(meandb_path, meanact_path)[!file.exists(c(meandb_path, meanact_path))]
  if (length(missing) > 0) {
    stop(sprintf("Missing mechanistic layer file(s) for GCM=%s, SSP=%s:\n  %s",
                 gcm, ssp, paste(missing, collapse = "\n  ")))
  }

  meandb_r  <- terra::rast(meandb_path)
  meanact_r <- terra::rast(meanact_path)
  names(meandb_r)  <- "deviation_mean"
  names(meanact_r) <- "Activity"

  list(meandb = meandb_r, meanActivity = meanact_r)
}
