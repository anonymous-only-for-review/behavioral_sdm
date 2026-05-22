################################################################################
## RCode_SDM_make_figures.R
##
## Generates publication-ready figures and raster exports for the SDM section
## of:
##
##   Rubalcaba, Fandos & Diaz. "Behavioural thermoregulation and microclimate
##   reshape climate-driven range forecasts."
##
## Three concatenated sections (run in order; all sections read from disk
## artefacts produced by RCode_SDM_correlative.R and RCode_SDM_projections_future.R):
##
##   1. Main and Extended-Data figures  (Paul Tol CVD-safe palette, MESS overlay,
##      soft agreement maps, latitudinal profiles, niche-divergence panels).
##   2. Supplementary figures           (Fig S2.2 PDP response curves,
##                                       Fig S2.3 model-performance boxplots).
##   3. Raster export bundle            (GeoTIFF + ASCII with manifest, README
##                                       and per-class area statistics).
##
## Author: Guillermo Fandos (gfandos@ucm.es)
################################################################################

source("RCode_SDM_helpers.R")
set.seed(MASTER_SEED)


# ============================================================================ #
# Section 1: Main and Extended-Data figures (formerly 07_figures.R)            #
# ============================================================================ #

# ============================================================================ #
# ADDITIONAL PACKAGES FOR ADVANCED FIGURES
# ============================================================================ #

for (pkg in c("ggspatial", "colorspace", "scales", "scico")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

# ============================================================================ #
# LOAD ALL UPSTREAM DATA — from rasters / RDS / CSV on disk
# ============================================================================ #

cat("\n=== Loading data from disk ===\n")

# -- Continuous prediction rasters --
r_corr_curr_mean <- terra::rast(file.path(OUTDIR, "rasters/current_correlative_mean.tif"))
r_corr_curr_sd   <- terra::rast(file.path(OUTDIR, "rasters/current_correlative_sd.tif"))
r_hyb_curr_mean  <- terra::rast(file.path(OUTDIR, "rasters/current_hybrid_mean.tif"))
r_hyb_curr_sd    <- terra::rast(file.path(OUTDIR, "rasters/current_hybrid_sd.tif"))
r_corr_fut_mean  <- terra::rast(file.path(OUTDIR, "rasters/future_correlative_mean.tif"))
r_corr_fut_sd    <- terra::rast(file.path(OUTDIR, "rasters/future_correlative_sd.tif"))
r_hyb_fut_mean   <- terra::rast(file.path(OUTDIR, "rasters/future_hybrid_mean.tif"))
r_hyb_fut_sd     <- terra::rast(file.path(OUTDIR, "rasters/future_hybrid_sd.tif"))

# -- Binary maps (maxSSS threshold) --
bin_corr_curr_maxSSS <- terra::rast(file.path(OUTDIR, "rasters/binary_current_correlative_maxSSS_PA.tif"))
bin_hyb_curr_maxSSS  <- terra::rast(file.path(OUTDIR, "rasters/binary_current_hybrid_maxSSS_PA.tif"))
bin_corr_fut_maxSSS  <- terra::rast(file.path(OUTDIR, "rasters/binary_future_correlative_maxSSS_PA.tif"))
bin_hyb_fut_maxSSS   <- terra::rast(file.path(OUTDIR, "rasters/binary_future_hybrid_maxSSS_PA.tif"))

# -- Pre-computed difference and categorical rasters --
diff_present_cont <- terra::rast(file.path(OUTDIR, "rasters/diff_continuous_present.tif"))
diff_future_cont  <- terra::rast(file.path(OUTDIR, "rasters/diff_continuous_future.tif"))
pres_cat_maxSSS   <- terra::rast(file.path(OUTDIR, "rasters/categorical_present.tif"))
fut_cat_maxSSS    <- terra::rast(file.path(OUTDIR, "rasters/categorical_future.tif"))

# -- Land mask --
land_mask <- terra::rast(file.path(OUTDIR, "rasters/land_mask.tif"))

# -- Thresholds --
thr_data <- readRDS(file.path(OUTDIR, "performance_metrics/thresholds_maxSSS.rds"))
thr_corr_maxSSS <- thr_data$correlativo_mean
thr_hyb_maxSSS  <- thr_data$hibrido_mean

# -- Performance by fold (for Fig 2) --
perf_by_fold <- read.csv(file.path(OUTDIR, "performance_metrics/performance_by_fold.csv"),
                         stringsAsFactors = FALSE)

# -- Ensemble model objects (for Fig 1 var importance, Fig 3 PDP) --
ens_corr_path <- file.path(OUTDIR, "models/ensemble_correlative.rds")
ens_hyb_path  <- file.path(OUTDIR, "models/ensemble_hybrid.rds")
has_ensembles <- file.exists(ens_corr_path) && file.exists(ens_hyb_path)

if (has_ensembles) {
  mod_corr <- readRDS(ens_corr_path)
  mod_hyb  <- readRDS(ens_hyb_path)
  cat("  Ensemble models loaded (var_importance + GBM objects available).\n")
} else {
  mod_corr <- mod_hyb <- NULL
  cat("  Ensemble models NOT found — Fig 1 (var importance) and Fig 3 (PDP) will be skipped.\n")
}

# -- Training data for PDP (Fig 3) --
dat_hyb_path <- file.path(OUTDIR, "models/dat_hyb_compact.rds")
if (file.exists(dat_hyb_path)) {
  dat_hyb <- readRDS(dat_hyb_path)
} else {
  dat_hyb <- NULL
}

# -- Range change metrics (for Fig 4) --
range_change_path <- file.path(OUTDIR, "range_change_metrics.csv")
if (file.exists(range_change_path)) {
  area_change <- read.csv(range_change_path, stringsAsFactors = FALSE)
  cat("  Range change metrics loaded.\n")
} else {
  area_change <- NULL
  cat("  range_change_metrics.csv NOT found — Fig 4 (range shift) will be skipped.\n")
}

cat("  Data loaded successfully.\n")

cat("\n=== Generating publication-ready figures ===\n")

# ============================================================================ #
# CONSTANTS
# ============================================================================ #

MESS_THR <- -10
W2C  <- 180; HLOW <- 90; HMID <- 120; HHI <- 140

# ============================================================================ #
# CATEGORICAL PALETTE: Paul Tol "Bright" (CVD-safe, maximum contrast)
# Reference: Tol, P. (2021). Colour Schemes. Technical Note SRON/EPS/TN/09-002
# ============================================================================ #

cols_cat <- c(
  "1" = "#0077BB",   # Correlative only: Blue
  "2" = "#EE7733",   # Hybrid only: Orange
  "3" = "#882255"    # Both: Dark magenta
)

# Background colors
OCEAN_COL   <- "#DFF1FB"
LAND_NA_COL <- "#FFFFFF"
EXTRAP_COL  <- "#E5E5E5"

# Model colors for legends/histograms
COL_CORR <- "#0072B2"
COL_HYB  <- "#D55E00"

# ============================================================================ #
# THEME: Nature Communications
# ============================================================================ #

theme_nc <- function() {
  theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text        = element_text(size = 8),
      axis.title       = element_text(size = 9)
    )
}

# X-axis labels without decimals
axis_x_no_decimals <- scale_x_continuous(labels = scales::label_number(accuracy = 1))

# ============================================================================ #
# HELPER: save_figure() — exports both PDF (cairo_pdf) and PNG
# ============================================================================ #

save_figure <- function(plot_obj, filename_stub, width_mm, height_mm, dpi = 300) {
  dir.create(file.path(OUTDIR, "figures"), showWarnings = FALSE, recursive = TRUE)

  ggsave(
    file.path(OUTDIR, "figures", paste0(filename_stub, ".pdf")),
    plot_obj, width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf
  )

  ggsave(
    file.path(OUTDIR, "figures", paste0(filename_stub, ".png")),
    plot_obj, width = width_mm, height = height_mm, units = "mm", dpi = dpi
  )

  message("  Saved: ", filename_stub)
}

# ============================================================================ #
# HELPER: get_coast_cropped() — coastline, CRS transform, spatial crop
# ============================================================================ #

get_coast_cropped <- function(reference_raster) {
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    message("  Package 'rnaturalearth' not available.")
    return(NULL)
  }

  coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
  crs_r <- sf::st_crs(terra::crs(reference_raster))
  coast_tr <- try(suppressWarnings(sf::st_transform(coast, crs = crs_r)), silent = TRUE)

  if (inherits(coast_tr, "try-error")) return(NULL)

  # Spatial crop to raster extent
  ext <- terra::ext(reference_raster)
  bbox <- sf::st_bbox(c(xmin = ext[1], xmax = ext[2],
                        ymin = ext[3], ymax = ext[4]), crs = crs_r)

  coast_cropped <- try(sf::st_crop(coast_tr, bbox), silent = TRUE)

  if (inherits(coast_cropped, "try-error")) return(coast_tr)

  coast_cropped
}

# ============================================================================ #
# LOAD COASTLINE
# ============================================================================ #

message("Processing coastline...")
# Use the first available continuous prediction as reference raster
ref_raster <- r_corr_curr_mean
coast_sf   <- get_coast_cropped(ref_raster)

# ============================================================================ #
# LOAD MESS RASTERS (if available)
# ============================================================================ #

mess_corr_path <- file.path(OUTDIR, "rasters/mess_correlative_future_pointsRef.tif")
mess_hyb_path  <- file.path(OUTDIR, "rasters/mess_hybrid_future_pointsRef.tif")
has_mess <- file.exists(mess_corr_path) && file.exists(mess_hyb_path)

if (has_mess) {
  message("Loading MESS layers...")
  r_mess_corr <- terra::rast(mess_corr_path)
  r_mess_hyb  <- terra::rast(mess_hyb_path)
} else {
  message("  MESS layers not found; overlay/masking will be skipped.")
}

# ============================================================================ #
# ADVANCED PLOTTERS (override basic versions from 00_config.R)
# ============================================================================ #

# --- plot_categorical_nc() -------------------------------------------------
# Full implementation with coastline, MESS overlay, scale bar, Paul Tol palette

plot_categorical_nc <- function(
    r_cat, land_mask, title = NULL,
    ocean_col = OCEAN_COL, land_na_col = LAND_NA_COL,
    coast = NULL, coast_color = "grey30", coast_width = 0.3,
    mess = NULL, mess_thr = MESS_THR, extrap_col = EXTRAP_COL,
    add_scalebar = TRUE
) {
  land_mask_aligned <- align_mask_to(land_mask, r_cat)
  if (!is.null(mess)) mess <- align_mask_to(mess, r_cat)

  ex <- terra::ext(r_cat)

  df_cat  <- terra::as.data.frame(r_cat,             xy = TRUE, na.rm = FALSE)
  df_land <- terra::as.data.frame(land_mask_aligned, xy = TRUE, na.rm = FALSE)

  stopifnot(nrow(df_cat) == nrow(df_land))

  names(df_cat)[3]  <- "value"
  names(df_land)[3] <- "land"
  df <- dplyr::bind_cols(df_cat, land = df_land$land)

  is_land <- !is.na(df$land)

  # Categorical variable
  df$cat <- NA_character_
  df$cat[df$value == 3 & !is.na(df$value)] <- "Both"
  df$cat[df$value == 2 & !is.na(df$value)] <- "Hybrid only"
  df$cat[df$value == 1 & !is.na(df$value)] <- "Correlative only"

  lvls <- c("Both", "Hybrid only", "Correlative only")
  df$cat <- factor(df$cat, levels = lvls)

  # Palette: Both (magenta), Hybrid only (orange), Correlative only (blue)
  pal <- unname(c(cols_cat["3"], cols_cat["2"], cols_cat["1"]))

  p <- ggplot() +
    geom_raster(data = df[!is_land, ], aes(x = x, y = y), fill = ocean_col) +
    geom_raster(data = df[is_land & is.na(df$cat), ], aes(x = x, y = y), fill = land_na_col) +
    geom_raster(data = df[!is.na(df$cat), ], aes(x = x, y = y, fill = cat)) +
    scale_fill_manual(
      values = pal,
      limits = lvls,
      breaks = lvls,
      drop   = FALSE,
      name   = "Model agreement",
      guide  = guide_legend(override.aes = list(alpha = 1))
    )

  # MESS mask
  if (!is.null(mess)) {
    df_m <- terra::as.data.frame(mess, xy = TRUE, na.rm = FALSE)
    names(df_m) <- c("x", "y", "mess")

    p <- p + geom_raster(
      data = df_m[df_m$mess < mess_thr & !is.na(df_m$mess), , drop = FALSE],
      aes(x = x, y = y),
      fill = extrap_col,
      inherit.aes = FALSE
    )
  }

  # Coastline
  if (!is.null(coast) && inherits(coast, "sf")) {
    p <- p + geom_sf(
      data = coast,
      inherit.aes = FALSE,
      color = coast_color,
      linewidth = coast_width
    )
  }

  # Scale bar
  if (add_scalebar && !is.null(coast)) {
    p <- p + ggspatial::annotation_scale(
      location   = "br",
      width_hint = 0.2,
      style      = "ticks",
      line_width = 0.5,
      height     = unit(0.15, "cm")
    )
  }

  p <- p +
    coord_sf(
      xlim   = c(ex[1], ex[2]),
      ylim   = c(ex[3], ex[4]),
      expand = FALSE,
      crs    = sf::st_crs(terra::crs(r_cat))
    ) +
    labs(title = title, x = NULL, y = NULL) +
    axis_x_no_decimals +
    theme_nc()

  p
}

# --- plot_diff_continuous() -------------------------------------------------
# Continuous difference map with scico "vik" palette, coastline, MESS overlay

plot_diff_continuous <- function(
    r_diff, land_mask, title = NULL, lim = NULL,
    coast = NULL, coast_color = "grey30", coast_width = 0.3,
    mess = NULL, mess_thr = MESS_THR, extrap_col = EXTRAP_COL,
    add_scalebar = TRUE
) {
  r_diff <- terra::mask(r_diff, land_mask)
  if (!is.null(mess)) mess <- align_mask_to(mess, r_diff)

  df <- terra::as.data.frame(r_diff, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", "value")
  df$value[!is.finite(df$value)] <- NA

  if (is.null(lim)) lim <- max(abs(range(df$value, na.rm = TRUE)))

  ex <- terra::ext(r_diff)

  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster()

  if (requireNamespace("scico", quietly = TRUE)) {
    p <- p + scico::scale_fill_scico(
      palette = "vik",
      limits  = c(-lim, lim),
      oob     = scales::squish,
      name    = expression(Delta * " suitability\n(Hybrid - Correlative)")
    )
  } else {
    p <- p + scale_fill_gradientn(
      colours = colorspace::divergingx_hcl(11, palette = "Blue-Red 3"),
      limits  = c(-lim, lim),
      oob     = scales::squish,
      name    = expression(Delta * " suitability\n(Hybrid - Correlative)")
    )
  }

  if (!is.null(mess)) {
    dfm <- terra::as.data.frame(mess, xy = TRUE, na.rm = FALSE)
    names(dfm) <- c("x", "y", "mess")

    p <- p + geom_raster(
      data = dfm[dfm$mess < mess_thr & !is.na(dfm$mess), , drop = FALSE],
      aes(x = x, y = y),
      inherit.aes = FALSE,
      fill = extrap_col
    )
  }

  if (!is.null(coast) && inherits(coast, "sf")) {
    p <- p + geom_sf(
      data = coast,
      inherit.aes = FALSE,
      color = coast_color,
      linewidth = coast_width
    )
  }

  if (add_scalebar && !is.null(coast)) {
    p <- p + ggspatial::annotation_scale(
      location   = "br",
      width_hint = 0.2,
      style      = "ticks",
      line_width = 0.5,
      height     = unit(0.15, "cm")
    )
  }

  p <- p +
    coord_sf(
      xlim   = c(ex[1], ex[2]),
      ylim   = c(ex[3], ex[4]),
      expand = FALSE,
      crs    = sf::st_crs(terra::crs(r_diff))
    ) +
    labs(title = title, x = NULL, y = NULL) +
    axis_x_no_decimals +
    theme_nc()

  p
}

# --- plot_diff_soft_agreement() ---------------------------------------------
# Green/gray tint overlay with contour lines for model agreement

plot_diff_soft_agreement <- function(
    r_diff, suit_hyb, suit_corr, land_mask, title = NULL,
    global_lim = NULL,
    coast = NULL, coast_color = "grey30", coast_width = 0.3,
    mess = NULL, mess_thr = MESS_THR, extrap_col = EXTRAP_COL,
    green_tint = "#009E73", gray_tint = "#7F8C8D",
    max_tint_alpha_good = 0.20, max_tint_alpha_bad = 0.16,
    draw_contour_good = TRUE, contour_thr = 0.6,
    contour_col = "#006D2C", contour_alpha = 0.9, contour_lwd = 0.3,
    add_scalebar = TRUE
) {
  suit_hyb  <- align_mask_to(suit_hyb,  r_diff)
  suit_corr <- align_mask_to(suit_corr, r_diff)
  land_mask <- align_mask_to(land_mask, r_diff)
  if (!is.null(mess)) mess <- align_mask_to(mess, r_diff)

  r_diff    <- terra::mask(r_diff, land_mask)
  suit_hyb  <- terra::mask(suit_hyb,  land_mask)
  suit_corr <- terra::mask(suit_corr, land_mask)

  df_diff <- terra::as.data.frame(r_diff,    xy = TRUE, na.rm = FALSE)
  df_hyb  <- terra::as.data.frame(suit_hyb,  xy = TRUE, na.rm = FALSE)
  df_corr <- terra::as.data.frame(suit_corr, xy = TRUE, na.rm = FALSE)

  names(df_diff) <- c("x", "y", "diff")
  names(df_hyb)  <- c("x", "y", "hyb")
  names(df_corr) <- c("x", "y", "corr")

  df <- df_diff %>%
    dplyr::left_join(df_hyb,  by = c("x", "y")) %>%
    dplyr::left_join(df_corr, by = c("x", "y"))

  df$hyb  <- scales::squish(df$hyb,  c(0, 1))
  df$corr <- scales::squish(df$corr, c(0, 1))
  df$soft_good <- pmin(df$hyb, df$corr)
  df$soft_bad  <- pmin(1 - df$hyb, 1 - df$corr)
  df$alpha_good <- max_tint_alpha_good * df$soft_good
  df$alpha_bad  <- max_tint_alpha_bad  * df$soft_bad

  lim <- if (is.null(global_lim)) {
    max(abs(range(df$diff, na.rm = TRUE)))
  } else {
    global_lim
  }

  rs <- terra::res(r_diff)
  wx <- rs[1]; wy <- rs[2]
  ex <- terra::ext(r_diff)

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_raster(aes(fill = diff))

  if (requireNamespace("scico", quietly = TRUE)) {
    p <- p + scico::scale_fill_scico(
      palette = "vik",
      limits  = c(-lim, lim),
      oob     = scales::squish,
      name    = expression(Delta * " suitability\n(Hybrid - Correlative)")
    )
  } else {
    p <- p + scale_fill_gradientn(
      colours = colorspace::divergingx_hcl(11, palette = "Blue-Red 3"),
      limits  = c(-lim, lim),
      oob     = scales::squish,
      name    = expression(Delta * " suitability\n(Hybrid - Correlative)")
    )
  }

  p <- p + geom_raster(
    data = df,
    aes(x = x, y = y),
    fill = green_tint,
    alpha = df$alpha_good,
    inherit.aes = FALSE
  )

  p <- p + geom_raster(
    data = df,
    aes(x = x, y = y),
    fill = gray_tint,
    alpha = df$alpha_bad,
    inherit.aes = FALSE
  )

  if (isTRUE(draw_contour_good)) {
    df$both_high <- as.numeric(df$hyb >= contour_thr & df$corr >= contour_thr)
    p <- p + geom_contour(
      data = df,
      aes(z = both_high),
      breaks    = 0.5,
      colour    = contour_col,
      linewidth = contour_lwd,
      alpha     = contour_alpha
    )
  }

  if (!is.null(mess)) {
    dfm <- terra::as.data.frame(mess, xy = TRUE, na.rm = FALSE)
    names(dfm) <- c("x", "y", "mess")

    p <- p + geom_raster(
      data = dfm[dfm$mess < mess_thr & !is.na(dfm$mess), , drop = FALSE],
      aes(x = x, y = y),
      fill = extrap_col,
      inherit.aes = FALSE
    )
  }

  if (!is.null(coast) && inherits(coast, "sf")) {
    p <- p + geom_sf(
      data = coast,
      inherit.aes = FALSE,
      color = coast_color,
      linewidth = coast_width
    )
  }

  if (add_scalebar && !is.null(coast)) {
    p <- p + ggspatial::annotation_scale(
      location   = "br",
      width_hint = 0.2,
      style      = "ticks",
      line_width = 0.5,
      height     = unit(0.15, "cm")
    )
  }

  p <- p +
    coord_sf(
      xlim   = c(ex[1], ex[2]),
      ylim   = c(ex[3], ex[4]),
      expand = FALSE,
      crs    = sf::st_crs(terra::crs(r_diff))
    ) +
    labs(title = title, x = NULL, y = NULL) +
    axis_x_no_decimals +
    theme_nc()

  p
}

# --- lat_profile() ----------------------------------------------------------
# Latitudinal binning of delta-suitability for profiles

lat_profile <- function(r_diff, land_mask,
                        lat_min = 36, lat_max = 43, bin_width = 0.25,
                        mess = NULL, mess_thr = -10, apply_mess = FALSE,
                        weight_by_area = FALSE) {
  # Align and mask
  r_diff <- terra::mask(r_diff, land_mask)
  if (apply_mess && !is.null(mess)) {
    mess_al <- align_mask_to(mess, r_diff)
    # Keep only cells with MESS >= threshold (no severe extrapolation)
    r_diff <- terra::mask(r_diff, mess_al >= mess_thr)
  }
  # Ensure lon/lat (WGS84); reproject if needed
  if (!grepl("longlat", terra::crs(r_diff), ignore.case = TRUE)) {
    r_diff <- terra::project(r_diff, "EPSG:4326", method = "bilinear")
  }
  # Extract data
  df <- terra::as.data.frame(r_diff, xy = TRUE, na.rm = FALSE)
  names(df) <- c("lon", "lat", "delta")
  df <- df[is.finite(df$delta) & !is.na(df$lat) & !is.na(df$lon), , drop = FALSE]

  # Filter latitude range
  df <- dplyr::filter(df, lat >= lat_min, lat <= lat_max)
  if (!nrow(df)) return(tibble::tibble(lat_mid = numeric(), mean = numeric(), se = numeric(), n = integer()))

  # Latitudinal bins
  breaks <- seq(lat_min, lat_max, by = bin_width)
  if (tail(breaks, 1) < lat_max) breaks <- c(breaks, lat_max)
  df$lat_bin <- cut(df$lat, breaks = breaks, include.lowest = TRUE, right = FALSE)

  # Bin centers for X axis
  lat_mids <- sapply(levels(df$lat_bin), function(lb) {
    bb <- strsplit(gsub("\\[|\\)|\\]", "", lb), ",")[[1]]
    bb <- as.numeric(bb)
    mean(bb)
  })
  df_mid <- dplyr::mutate(df, lat_mid = lat_mids[lat_bin])

  # Optional area weighting by cos(lat)
  if (isTRUE(weight_by_area)) {
    df_mid$w <- cos(df_mid$lat * pi / 180)
    df_mid$w[!is.finite(df_mid$w) | df_mid$w < 0] <- 0
  } else {
    df_mid$w <- 1
  }

  # Summary by bin
  prof <- df_mid %>%
    dplyr::group_by(lat_bin, lat_mid) %>%
    dplyr::summarise(
      mean = weighted.mean(delta, w, na.rm = TRUE),
      n    = dplyr::n(),
      sd   = stats::sd(delta, na.rm = TRUE),
      se   = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    dplyr::arrange(lat_mid)

  prof
}

# ============================================================================ #
#                           FIGURE GENERATION
# ============================================================================ #

# ---- Fig 1: Variable importance --------------------------------------------

if (!is.null(mod_corr) && !is.null(mod_hyb) &&
    !is.null(mod_corr$var_importance) && !is.null(mod_hyb$var_importance)) {

  message("\n=== Figure 1: Variable importance ===")

  var_imp_combined <- dplyr::bind_rows(
    mod_corr$var_importance %>% dplyr::mutate(Model = "Correlative"),
    mod_hyb$var_importance  %>% dplyr::mutate(Model = "Hybrid")
  ) %>%
    dplyr::mutate(
      Model = factor(Model, levels = c("Correlative", "Hybrid")),
      Type  = ifelse(Variable %in% c("deviation_mean", "Activity"), "Mechanistic", "Climatic")
    )

  fig1 <- ggplot(var_imp_combined, aes(x = RelInf_mean, y = reorder(Variable, RelInf_mean))) +
    geom_col(aes(fill = Type), width = 0.7) +
    geom_errorbar(aes(xmin = RelInf_mean - RelInf_sd, xmax = RelInf_mean + RelInf_sd),
                  width = 0.2, linewidth = 0.3) +
    facet_wrap(~Model, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Mechanistic" = COL_HYB, "Climatic" = COL_CORR)) +
    labs(x = "Relative influence (%)", y = NULL, fill = NULL) +
    theme_nc() +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 11),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())

  save_figure(fig1, "fig1_variable_importance", W2C, HMID)

} else {
  message("  [SKIP] Fig 1 -- Ensemble models not available on disk")
}

# ---- Fig 2: Model performance boxplots ------------------------------------
# ADAPTED: uses performance_metrics/performance_by_fold.csv instead of
# thr_results_corr_p10$eval_summary / thr_results_hyb_p10$eval_summary

message("\n=== Figure 2: Model performance boxplots ===")

perf_combined <- perf_by_fold %>%
  dplyr::mutate(Model = factor(Model, levels = c("Correlative", "Hybrid")))

# Check that required columns exist
if (all(c("AUC", "TSS", "Model") %in% names(perf_combined))) {
  fig2 <- perf_combined %>%
    tidyr::pivot_longer(cols = c(AUC, TSS), names_to = "Metric", values_to = "Value") %>%
    ggplot(aes(x = Model, y = Value, fill = Model)) +
    geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
    facet_wrap(~Metric, scales = "free_y") +
    scale_fill_manual(values = c("Correlative" = COL_CORR, "Hybrid" = COL_HYB)) +
    labs(y = NULL, x = NULL) +
    theme_nc() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 11),
          panel.grid.major.x = element_blank())

  save_figure(fig2, "fig2_model_performance", HMID, HLOW)
} else {
  message("  [SKIP] Fig 2 -- performance_by_fold.csv missing AUC/TSS/Model columns")
}

# ---- Fig 3: Response curves PDP -------------------------------------------

if (!is.null(mod_hyb) && !is.null(mod_hyb$models) && !is.null(dat_hyb)) {

  # Check that all required predictors are present in dat_hyb
  missing_vars <- setdiff(mod_hyb$predictors, names(dat_hyb))

  if (length(missing_vars) == 0) {

    message("\n=== Figure 3: Response curves (PDP) ===")

    compute_pdp_with_uncertainty <- function(models, variable, data, n_points = 50) {
      pdp_list <- lapply(seq_along(models), function(i) {
        m <- models[[i]]
        nt <- if (!is.null(m$gbm.call$best.trees)) m$gbm.call$best.trees else m$n.trees
        pd <- pdp::partial(m, pred.var = variable, train = data,
                           n.trees = nt, type = "classification",
                           grid.resolution = n_points)
        data.frame(model_id = i, var_value = pd[[variable]], yhat = pd$yhat)
      })
      dplyr::bind_rows(pdp_list) %>%
        dplyr::group_by(var_value) %>%
        dplyr::summarise(
          yhat_mean  = mean(yhat, na.rm = TRUE),
          yhat_sd    = sd(yhat, na.rm = TRUE),
          yhat_lower = quantile(yhat, 0.025, na.rm = TRUE),
          yhat_upper = quantile(yhat, 0.975, na.rm = TRUE),
          .groups    = "drop"
        )
    }

    pdp_deviation <- compute_pdp_with_uncertainty(mod_hyb$models, "deviation_mean", dat_hyb)
    pdp_activity  <- compute_pdp_with_uncertainty(mod_hyb$models, "Activity",       dat_hyb)

    fig3a <- ggplot(pdp_deviation, aes(x = var_value, y = yhat_mean)) +
      geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper), fill = COL_HYB, alpha = 0.2) +
      geom_line(linewidth = 1.2, color = COL_HYB) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.5) +
      labs(x = "Thermoregulatory inaccuracy (\u00B0C)", y = "Predicted suitability") +
      theme_nc() + theme(panel.grid.minor = element_blank())

    fig3b <- ggplot(pdp_activity, aes(x = var_value, y = yhat_mean)) +
      geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper), fill = COL_CORR, alpha = 0.2) +
      geom_line(linewidth = 1.2, color = COL_CORR) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.5) +
      labs(x = "Activity window (h/day)", y = "Predicted suitability") +
      theme_nc() + theme(panel.grid.minor = element_blank())

    fig3 <- fig3a + fig3b
    save_figure(fig3, "fig3_response_curves_with_uncertainty", W2C, 70)

  } else {
    message("  [SKIP] Fig 3 -- dat_hyb_compact.rds missing predictors: ", paste(missing_vars, collapse = ", "))
    message("         To enable: download bio15 raster and rebuild dat_hyb_compact.rds")
  }

} else {
  message("  [SKIP] Fig 3 -- Response curves (PDP) require ensemble model objects and training data")
}

# ---- Fig 4: Range shift bars ----------------------------------------------

if (!is.null(area_change) && all(c("Model", "Gained_km2", "Lost_km2") %in% names(area_change))) {

  message("\n=== Figure 4: Range shift magnitude ===")

  area_change_long <- area_change %>%
    tidyr::pivot_longer(cols = c(Gained_km2, Lost_km2),
                        names_to = "Change_type", values_to = "Area_km2") %>%
    dplyr::mutate(
      Area_km2    = ifelse(Change_type == "Lost_km2", -Area_km2, Area_km2),
      Change_type = factor(Change_type,
                           levels = c("Gained_km2", "Lost_km2"),
                           labels = c("Expansion", "Contraction"))
    )

  fig4 <- ggplot(area_change_long, aes(x = Model, y = Area_km2, fill = Change_type)) +
    geom_col(position = "identity", alpha = 0.8, width = 0.6) +
    scale_fill_manual(values = c("Expansion" = COL_CORR, "Contraction" = COL_HYB)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    labs(y = expression("Area change (km"^2*")"), x = NULL, fill = NULL) +
    scale_y_continuous(labels = scales::comma) +
    theme_nc() +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

  save_figure(fig4, "fig4_range_shift_magnitude", 100, 100)

} else {
  message("  [SKIP] Fig 4 -- range_change_metrics.csv not found or missing columns")
}

# ============================================================================ #
# Fig 5: CATEGORICAL MAPS (overlay + masked)
# ============================================================================ #

message("\n=== Figure 5: Model comparison (categorical) ===")

# Use pre-computed categorical rasters loaded from disk
# (pres_cat_maxSSS and fut_cat_maxSSS already loaded above)

# Fig 5A: Overlay
fig5_cat_overlay <- patchwork::wrap_plots(
  plot_categorical_nc(
    pres_cat_maxSSS, land_mask,
    title = "Present: Hybrid vs Correlative (maxSSS)",
    coast = coast_sf
  ),
  plot_categorical_nc(
    fut_cat_maxSSS, land_mask,
    title = "Future: Hybrid vs Correlative (MESS overlay)",
    coast = coast_sf,
    mess     = if (has_mess) r_mess_hyb else NULL,
    mess_thr = MESS_THR
  ),
  ncol   = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Figure 5: Model agreement under present and future climate",
    theme = theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
  )

save_figure(fig5_cat_overlay, "fig5_categorical_maxSSS_OVERLAY", W2C, HLOW)

# Fig 5B: Masked
if (has_mess) {
  message("  Generating masked version...")

  fut_mask_corr <- terra::mask(r_corr_fut_mean, r_mess_corr >= MESS_THR)
  fut_mask_hyb  <- terra::mask(r_hyb_fut_mean,  r_mess_hyb  >= MESS_THR)

  bin_corr_fut_mask <- terra::ifel(is.na(fut_mask_corr), NA,
                                   fut_mask_corr >= thr_corr_maxSSS)
  bin_hyb_fut_mask  <- terra::ifel(is.na(fut_mask_hyb), NA,
                                   fut_mask_hyb >= thr_hyb_maxSSS)

  fut_cat_maxSSS_MASKED <- cat_map_from(bin_corr_fut_mask, bin_hyb_fut_mask)

  fig5_cat_masked <- patchwork::wrap_plots(
    plot_categorical_nc(
      pres_cat_maxSSS, land_mask,
      title = "Present: Hybrid vs Correlative (maxSSS)",
      coast = coast_sf
    ),
    plot_categorical_nc(
      fut_cat_maxSSS_MASKED, land_mask,
      title = "Future: Hybrid vs Correlative (MESS masked)",
      coast = coast_sf
    ),
    ncol   = 2,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      theme = theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

  save_figure(fig5_cat_masked, "fig5_categorical_maxSSS_MASKED", W2C, HLOW)
}

# ============================================================================ #
# CONTINUOUS DIFFERENCE MAPS (overlay + masked)
# ============================================================================ #

message("\n=== Continuous difference maps ===")

# Use pre-computed difference rasters loaded from disk
# (diff_present_cont and diff_future_cont already loaded above)

rng <- max(
  abs(terra::global(diff_present_cont, "max", na.rm = TRUE)[[1]]),
  abs(terra::global(diff_present_cont, "min", na.rm = TRUE)[[1]]),
  abs(terra::global(diff_future_cont,  "max", na.rm = TRUE)[[1]]),
  abs(terra::global(diff_future_cont,  "min", na.rm = TRUE)[[1]])
)

p_diff_pres <- plot_diff_continuous(
  diff_present_cont, land_mask,
  title = "Present: Hybrid - Correlative",
  lim   = rng,
  coast = coast_sf
)

p_diff_fut <- plot_diff_continuous(
  diff_future_cont, land_mask,
  title = "Future: Hybrid - Correlative (MESS overlay)",
  lim   = rng,
  coast = coast_sf,
  mess     = if (has_mess) r_mess_hyb else NULL,
  mess_thr = MESS_THR
)

fig_diff_cont <- patchwork::wrap_plots(
  p_diff_pres, p_diff_fut,
  ncol   = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    theme = theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
  )

save_figure(fig_diff_cont, "ED_fig_continuous_diff_OVERLAY", W2C, HLOW)

if (has_mess) {
  diff_future_cont_MASKED <- terra::mask(diff_future_cont, r_mess_hyb >= MESS_THR)

  p_diff_fut_mask <- plot_diff_continuous(
    diff_future_cont_MASKED, land_mask,
    title = "Future: Hybrid - Correlative (MESS masked)",
    lim   = rng,
    coast = coast_sf
  )

  fig_diff_cont_mask <- patchwork::wrap_plots(
    p_diff_pres, p_diff_fut_mask,
    ncol   = 2,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      theme = theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

  save_figure(fig_diff_cont_mask, "ED_fig_continuous_diff_MASKED", W2C, HLOW)
}

# ============================================================================ #
# SOFT AGREEMENT MAPS (overlay + masked)
# ============================================================================ #

message("\n=== Soft agreement maps ===")

compute_global_lim <- function(r_list, land_mask) {
  vals <- unlist(lapply(r_list, function(r) {
    v <- terra::values(terra::mask(r, land_mask), mat = FALSE)
    v[is.finite(v)]
  }))
  max(abs(range(vals, na.rm = TRUE)))
}

# Use directly loaded mean suitability rasters
suit_hyb_pres  <- r_hyb_curr_mean
suit_corr_pres <- r_corr_curr_mean
suit_hyb_fut   <- r_hyb_fut_mean
suit_corr_fut  <- r_corr_fut_mean

# Use pre-computed difference rasters
diff_pres <- diff_present_cont
diff_fut  <- diff_future_cont

global_lim <- compute_global_lim(list(diff_pres, diff_fut), land_mask)

p_pres_soft <- plot_diff_soft_agreement(
  r_diff    = diff_pres,
  suit_hyb  = suit_hyb_pres,
  suit_corr = suit_corr_pres,
  land_mask = land_mask,
  title      = "Present: Agreement overlay",
  global_lim = global_lim,
  coast      = coast_sf
)

p_fut_soft <- plot_diff_soft_agreement(
  r_diff    = diff_fut,
  suit_hyb  = suit_hyb_fut,
  suit_corr = suit_corr_fut,
  land_mask = land_mask,
  title      = "Future: Agreement overlay (MESS)",
  global_lim = global_lim,
  coast      = coast_sf,
  mess       = if (has_mess) r_mess_hyb else NULL,
  mess_thr   = MESS_THR
)

fig_soft <- patchwork::wrap_plots(
  p_pres_soft, p_fut_soft,
  ncol   = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    theme = theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
  )

save_figure(fig_soft, "ED_fig_soft_agreement_OVERLAY", W2C, HLOW)

if (has_mess) {
  suit_hyb_fut_MASK  <- terra::mask(suit_hyb_fut,  r_mess_hyb  >= MESS_THR)
  suit_corr_fut_MASK <- terra::mask(suit_corr_fut, r_mess_corr >= MESS_THR)
  diff_fut_MASK      <- suit_hyb_fut_MASK - suit_corr_fut_MASK

  p_fut_soft_mask <- plot_diff_soft_agreement(
    r_diff    = diff_fut_MASK,
    suit_hyb  = suit_hyb_fut_MASK,
    suit_corr = suit_corr_fut_MASK,
    land_mask = land_mask,
    title      = "Future: Agreement overlay (masked)",
    global_lim = global_lim,
    coast      = coast_sf
  )

  fig_soft_mask <- patchwork::wrap_plots(
    p_pres_soft, p_fut_soft_mask,
    ncol   = 2,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      theme = theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

  save_figure(fig_soft_mask, "ED_fig_soft_agreement_MASKED", W2C, HLOW)
}

# ============================================================================ #
# Fig 8: MESS MAPS + HISTOGRAMS (simplified, no percentage text)
# ============================================================================ #

if (has_mess) {
  message("\n=== Figure 8: MESS diagnostic figures ===")

  plot_mess_map <- function(r_mess, title = NULL) {
    df <- terra::as.data.frame(r_mess, xy = TRUE, na.rm = FALSE)
    names(df) <- c("x", "y", "value")
    df$value[!is.finite(df$value)] <- NA

    ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_gradient2(
        low      = "#B2182B",
        mid      = "white",
        high     = "#2166AC",
        midpoint = 0,
        na.value = "grey92",
        name     = "MESS"
      ) +
      coord_equal(expand = FALSE) +
      labs(title = title, x = NULL, y = NULL) +
      axis_x_no_decimals +
      theme_nc()
  }

  p_map_corr <- plot_mess_map(r_mess_corr, "Correlative: MESS (Future vs Present)")
  p_map_hyb  <- plot_mess_map(r_mess_hyb,  "Hybrid: MESS (Future vs Present)")

  # Extract values
  mess_vals_corr <- terra::values(r_mess_corr, mat = FALSE)
  mess_vals_hyb  <- terra::values(r_mess_hyb,  mat = FALSE)

  df_mess <- dplyr::bind_rows(
    data.frame(Model = "Correlative", value = mess_vals_corr),
    data.frame(Model = "Hybrid",      value = mess_vals_hyb)
  ) %>%
    tidyr::drop_na()

  # Calculate extrapolation percentages (for export only, not displayed)
  perc_extrap <- df_mess %>%
    dplyr::mutate(Extrap = value < 0) %>%
    dplyr::group_by(Model) %>%
    dplyr::summarise(Extrap_pct = 100 * mean(Extrap), .groups = "drop")

  message("  Extrapolation summary (MESS < 0):")
  message("    Correlative: ", round(perc_extrap$Extrap_pct[1], 2), "%")
  message("    Hybrid: ", round(perc_extrap$Extrap_pct[2], 2), "%")

  # SIMPLIFIED HISTOGRAM - NO TEXT ANNOTATIONS
  plot_mess_hist_clean <- function(df, xlab = "MESS value") {
    ggplot(df, aes(x = value, fill = Model)) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins     = 60,
        alpha    = 0.6,
        position = "identity"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
      scale_fill_manual(
        values = c("Correlative" = COL_CORR, "Hybrid" = COL_HYB)
      ) +
      labs(x = xlab, y = "Density", fill = NULL) +
      theme_nc() +
      theme(legend.position = "top")
  }

  # Figure 8: Maps + clean histogram
  fig8_mess <- (p_map_corr + p_map_hyb) /
    plot_mess_hist_clean(df_mess) +
    patchwork::plot_layout(heights = c(1, 0.8)) +
    patchwork::plot_annotation(
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

  save_figure(fig8_mess, "fig8_mess_maps_plus_hist", W2C, HHI)

  # Export percentages to CSV
  readr::write_csv(
    perc_extrap,
    file.path(OUTDIR, "figures/mess_extrapolation_percentages.csv")
  )

  # Grouped histogram (values < -10 binned)
  cutoff <- -10
  df_mess_bucket <- df_mess %>%
    dplyr::mutate(value = ifelse(value < cutoff, cutoff, value))

  perc_lt_cut <- df_mess %>%
    dplyr::group_by(Model) %>%
    dplyr::summarise(lt_cut = 100 * mean(value < cutoff), .groups = "drop")

  message("  Severe extrapolation (MESS < ", cutoff, "):")
  message("    Correlative: ", round(perc_lt_cut$lt_cut[1], 2), "%")
  message("    Hybrid: ", round(perc_lt_cut$lt_cut[2], 2), "%")

  p_hist_bucket_clean <- plot_mess_hist_clean(
    df_mess_bucket,
    xlab = paste0("MESS (values < ", cutoff, " grouped)")
  ) +
    geom_vline(xintercept = cutoff, linetype = "dotted", linewidth = 0.5)

  fig8_mess_bucket <- (p_map_corr + p_map_hyb) / p_hist_bucket_clean +
    patchwork::plot_layout(heights = c(1, 0.8)) +
    patchwork::plot_annotation(
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

  save_figure(fig8_mess_bucket, "fig8c_mess_hist_grouped", W2C, HHI)

  # Export severe extrapolation stats
  readr::write_csv(
    perc_lt_cut,
    file.path(OUTDIR, "figures/mess_severe_extrapolation_percentages.csv")
  )
}

# ============================================================================ #
# LATITUDINAL DELTA-SUITABILITY PROFILES
# ============================================================================ #

message("\n=== Latitudinal delta-suitability profiles ===")

# Compute profiles using pre-computed difference rasters
prof_pres <- lat_profile(
  r_diff    = diff_present_cont,
  land_mask = land_mask,
  lat_min   = 36, lat_max = 43, bin_width = 0.25,
  mess = NULL, apply_mess = FALSE, weight_by_area = FALSE
)

prof_fut <- lat_profile(
  r_diff    = diff_future_cont,
  land_mask = land_mask,
  lat_min   = 36, lat_max = 43, bin_width = 0.25,
  mess      = if (has_mess) r_mess_hyb else NULL,
  mess_thr  = -10,
  apply_mess     = has_mess,
  weight_by_area = FALSE
)

# Degree formatter for latitude axis
fmt_deg <- function(x) paste0(scales::number(x, accuracy = 0.1), "\u00B0N")

# Line + SE ribbon version
p_lat <- ggplot() +
  # Future
  geom_ribbon(
    data = prof_fut,
    aes(x = lat_mid, ymin = mean - se, ymax = mean + se),
    alpha = 0.15, fill = COL_HYB
  ) +
  geom_line(
    data = prof_fut, aes(x = lat_mid, y = mean),
    linewidth = 1.0, color = COL_HYB
  ) +
  # Present
  geom_ribbon(
    data = prof_pres,
    aes(x = lat_mid, ymin = mean - se, ymax = mean + se),
    alpha = 0.15, fill = COL_CORR
  ) +
  geom_line(
    data = prof_pres, aes(x = lat_mid, y = mean),
    linewidth = 1.0, color = COL_CORR
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  scale_x_continuous(
    breaks = seq(36, 43, by = 1), labels = fmt_deg,
    limits = c(36, 43), expand = c(0.01, 0.01)
  ) +
  labs(
    x        = "Latitude",
    y        = expression(Delta ~ "suitability (Hybrid - Correlative)"),
    subtitle = if (has_mess) "Future masked where MESS < -10" else NULL
  ) +
  theme_nc() +
  theme(legend.position = "none")

save_figure(p_lat, "ED_fig_latitudinal_profile", W2C, 70)

# LOESS smoothed version with unified legend
prof_long <- dplyr::bind_rows(
  prof_pres %>% dplyr::mutate(Scenario = "Present"),
  prof_fut  %>% dplyr::mutate(Scenario = "Future")
)

loess_span <- 0.5

p_lat_loess <- ggplot(prof_long, aes(x = lat_mid, y = mean,
                                     color = Scenario, fill = Scenario)) +
  geom_smooth(
    method    = "loess",
    se        = TRUE,
    span      = loess_span,
    linewidth = 1.0,
    alpha     = 0.20,
    level     = 0.95,
    show.legend = TRUE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  scale_color_manual(values = c("Present" = COL_CORR, "Future" = COL_HYB), name = "Scenario") +
  scale_fill_manual(values  = c("Present" = COL_CORR, "Future" = COL_HYB), name = "Scenario") +
  scale_x_continuous(
    breaks = seq(36, 43, by = 1), labels = fmt_deg,
    limits = c(36, 43), expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Latitude",
    y = expression(Delta ~ "suitability (Hybrid - Correlative)")
  ) +
  theme_nc() +
  theme(
    legend.position = "bottom",
    legend.box      = "vertical",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8)
  )

save_figure(p_lat_loess, "ED_fig_latitudinal_profile_LOESS", W2C, 70)

# ============================================================================ #
# DONE
# ============================================================================ #

message("\n", strrep("=", 70))
message("  ALL FIGURES EXPORTED TO: ", file.path(OUTDIR, "figures"))
message(strrep("=", 70))
message("\n  Section 1 (main + ED figures) complete\n")


# ============================================================================ #
# Section 2: Supplementary figures (formerly 16_update_SI_figures.R)           #
#   S2.2: PDP response curves                                                  #
#   S2.3: Model-performance boxplots (50 fits)                                 #
# ============================================================================ #

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

set.seed(MASTER_SEED)

# ---- Theme and colors (matching 07_figures.R) --------------------------------

COL_CORR <- "#0072B2"
COL_HYB  <- "#D55E00"

theme_nc <- function() {
  theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text        = element_text(size = 8),
      axis.title       = element_text(size = 9)
    )
}

save_figure <- function(plot_obj, filename_stub, width_mm, height_mm, dpi = 300) {
  dir.create(file.path(OUTDIR, "figures"), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(OUTDIR, "figures", paste0(filename_stub, ".pdf")),
         plot_obj, width = width_mm, height = height_mm, units = "mm",
         device = cairo_pdf)
  ggsave(file.path(OUTDIR, "figures", paste0(filename_stub, ".png")),
         plot_obj, width = width_mm, height = height_mm, units = "mm",
         dpi = dpi)
  message("  Saved: ", filename_stub, " (.pdf + .png)")
}

# ==============================================================================
# Fig S2.2: Partial dependence plots — mechanistic variables
# ==============================================================================

message("\n=== Fig S2.2: Response curves (PDP) ===")

pdp_deviation <- read.csv(file.path(OUTDIR, "pdp_deviation_mean.csv"))
pdp_activity  <- read.csv(file.path(OUTDIR, "pdp_activity.csv"))

# Convert activity from minutes to hours for display
pdp_activity$var_value_h <- pdp_activity$var_value / 60

fig_s22a <- ggplot(pdp_deviation, aes(x = var_value, y = yhat_mean)) +
  geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper),
              fill = COL_HYB, alpha = 0.2) +
  geom_line(linewidth = 1.2, color = COL_HYB) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  labs(x = expression("Thermoregulatory inaccuracy (" * degree * "C)"),
       y = "Predicted suitability") +
  theme_nc() +
  theme(panel.grid.minor = element_blank())

fig_s22b <- ggplot(pdp_activity, aes(x = var_value_h, y = yhat_mean)) +
  geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper),
              fill = COL_CORR, alpha = 0.2) +
  geom_line(linewidth = 1.2, color = COL_CORR) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  labs(x = expression("Thermoregulatory window (h day"^{-1} * ")"),
       y = "Predicted suitability") +
  theme_nc() +
  theme(panel.grid.minor = element_blank())

fig_s22 <- fig_s22a + fig_s22b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

save_figure(fig_s22, "fig_S2.2_response_curves", 180, 70)

# ==============================================================================
# Fig S2.3: Model performance boxplots — 50 cross-validated fits
# ==============================================================================

message("\n=== Fig S2.3: Model performance (50 fits) ===")

perf_50 <- read.csv(file.path(OUTDIR, "model_performance_by_fold_50fits.csv"))

# Reshape: wide (AUC_corr, TSS_corr, AUC_hyb, TSS_hyb) -> long
perf_long <- bind_rows(
  perf_50 %>%
    dplyr::select(Repeat, fold, AUC = AUC_corr, TSS = TSS_corr) %>%
    mutate(Model = "Correlative"),
  perf_50 %>%
    dplyr::select(Repeat, fold, AUC = AUC_hyb, TSS = TSS_hyb) %>%
    mutate(Model = "Hybrid")
) %>%
  mutate(Model = factor(Model, levels = c("Correlative", "Hybrid"))) %>%
  pivot_longer(cols = c(AUC, TSS), names_to = "Metric", values_to = "Value")

fig_s23 <- ggplot(perf_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, shape = 16) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("Correlative" = COL_CORR, "Hybrid" = COL_HYB)) +
  labs(y = NULL, x = NULL) +
  theme_nc() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_blank())

save_figure(fig_s23, "fig_S2.3_model_performance_50fits", 120, 90)

# ---- Summary stats -----------------------------------------------------------

message("\n--- Performance summary (50 test folds) ---")
perf_long %>%
  group_by(Model, Metric) %>%
  summarise(
    mean = round(mean(Value), 3),
    sd   = round(sd(Value), 3),
    .groups = "drop"
  ) %>%
  print()

message("\n  Section 2 (SI figures) complete\n")


# ============================================================================ #
# Section 3: Raster export bundle (formerly 09_exports.R)                      #
#   GeoTIFF (LZW compressed) + ASCII (.asc.gz) + manifest + README + areas     #
# ============================================================================ #

# ============================================================================ #
# LOAD DATA FROM DISK (rasters + RDS)
# ============================================================================ #

cat("\n=== Loading data from disk ===\n")

# Continuous predictions
r_corr_curr_mean <- terra::rast(file.path(OUTDIR, "rasters/current_correlative_mean.tif"))
r_hyb_curr_mean  <- terra::rast(file.path(OUTDIR, "rasters/current_hybrid_mean.tif"))
r_corr_fut_mean  <- terra::rast(file.path(OUTDIR, "rasters/future_correlative_mean.tif"))
r_hyb_fut_mean   <- terra::rast(file.path(OUTDIR, "rasters/future_hybrid_mean.tif"))

# Binary maps
bin_corr_curr_maxSSS <- terra::rast(file.path(OUTDIR, "rasters/binary_current_correlative_maxSSS_PA.tif"))
bin_hyb_curr_maxSSS  <- terra::rast(file.path(OUTDIR, "rasters/binary_current_hybrid_maxSSS_PA.tif"))
bin_corr_fut_maxSSS  <- terra::rast(file.path(OUTDIR, "rasters/binary_future_correlative_maxSSS_PA.tif"))
bin_hyb_fut_maxSSS   <- terra::rast(file.path(OUTDIR, "rasters/binary_future_hybrid_maxSSS_PA.tif"))

# Pre-computed categorical maps
pres_cat <- terra::rast(file.path(OUTDIR, "rasters/categorical_present.tif"))
fut_cat  <- terra::rast(file.path(OUTDIR, "rasters/categorical_future.tif"))

# Pre-computed difference maps
diff_present_cont <- terra::rast(file.path(OUTDIR, "rasters/diff_continuous_present.tif"))
diff_future_cont  <- terra::rast(file.path(OUTDIR, "rasters/diff_continuous_future.tif"))

# Land mask
land_mask <- terra::rast(file.path(OUTDIR, "rasters/land_mask.tif"))

# Thresholds
thr_data        <- readRDS(file.path(OUTDIR, "performance_metrics/thresholds_maxSSS.rds"))
thr_corr_maxSSS <- thr_data$correlativo_mean
thr_hyb_maxSSS  <- thr_data$hibrido_mean

# Alias continuous predictions for downstream use
pred_corr_curr_mean <- r_corr_curr_mean
pred_hyb_curr_mean  <- r_hyb_curr_mean
pred_corr_fut_mean  <- r_corr_fut_mean
pred_hyb_fut_mean   <- r_hyb_fut_mean

cat("\n=== Export bundle: GeoTIFF + ASCII ===\n")

# ============================================================================ #
# DIRECTORIES
# ============================================================================ #

EXP_DIR <- file.path(OUTDIR, "exports")
ASC_DIR <- file.path(OUTDIR, "exports_asc")
dir.create(EXP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ASC_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================ #
# CONSTANTS
# ============================================================================ #

MESS_THR <- -10

# MESS layers
mess_corr_path <- file.path(OUTDIR, "rasters/mess_correlative_future_pointsRef.tif")
mess_hyb_path  <- file.path(OUTDIR, "rasters/mess_hybrid_future_pointsRef.tif")
has_mess <- file.exists(mess_corr_path) && file.exists(mess_hyb_path)

if (has_mess) {
  r_mess_corr <- terra::rast(mess_corr_path)
  r_mess_hyb  <- terra::rast(mess_hyb_path)
}

# ============================================================================ #
# HELPERS
# ============================================================================ #

# GeoTIFF (LZW compressed)
.write_gtiff <- function(r, fname, dtype = "FLT4S", nodata = NA) {
  terra::writeRaster(
    r, filename = file.path(EXP_DIR, fname),
    overwrite = TRUE, datatype = dtype, NAflag = nodata,
    gdal = c("COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
  cat("  GeoTIFF:", fname, "\n")
}

# ASCII (.asc) with driver fallback + optional gzip
.write_asc <- function(r, fname_base, nodata = -9999, gzip_after = TRUE) {
  fpath <- file.path(ASC_DIR, paste0(fname_base, ".asc"))
  ok <- TRUE
  tryCatch({
    terra::writeRaster(r, filename = fpath, overwrite = TRUE, NAflag = nodata)
  }, error = function(e) {
    ok <<- FALSE
    message("  Warning: ASCII attempt 1 failed: ", conditionMessage(e))
  })
  if (!ok) {
    tryCatch({
      terra::writeRaster(r, filename = fpath, overwrite = TRUE,
                         filetype = "AAIGrid", NAflag = nodata)
    }, error = function(e) {
      stop("ASCII export failed even with AAIGrid: ", conditionMessage(e))
    })
  }
  cat("  ASCII:", basename(fpath), "\n")
  if (isTRUE(gzip_after)) {
    utils::gzip(fpath, destname = paste0(fpath, ".gz"), overwrite = TRUE)
    cat("    -> gzipped:", paste0(basename(fpath), ".gz"), "\n")
  }
}

# ============================================================================ #
# PREPARE RASTERS
# ============================================================================ #

# Future MESS-masked version
fut_cat_MASKED <- NULL
if (isTRUE(has_mess)) {
  fut_mask_corr <- terra::mask(pred_corr_fut_mean, r_mess_corr >= MESS_THR)
  fut_mask_hyb  <- terra::mask(pred_hyb_fut_mean,  r_mess_hyb  >= MESS_THR)
  bin_corr_fut_mask <- terra::ifel(is.na(fut_mask_corr), NA, fut_mask_corr >= thr_corr_maxSSS)
  bin_hyb_fut_mask  <- terra::ifel(is.na(fut_mask_hyb),  NA, fut_mask_hyb  >= thr_hyb_maxSSS)
  fut_cat_MASKED <- cat_map_from(bin_corr_fut_mask, bin_hyb_fut_mask)
  fut_cat_MASKED <- terra::ifel(is.na(fut_cat_MASKED), NA, round(fut_cat_MASKED))
}

# Continuous difference: MESS-masked variant
diff_fut_MASKED <- if (isTRUE(has_mess)) terra::mask(diff_future_cont, r_mess_hyb >= MESS_THR) else NULL

# Clean versions (no MESS overlay/mask)
bin_corr_fut_nomask  <- terra::ifel(is.na(pred_corr_fut_mean), NA, pred_corr_fut_mean >= thr_corr_maxSSS)
bin_hyb_fut_nomask   <- terra::ifel(is.na(pred_hyb_fut_mean),  NA, pred_hyb_fut_mean  >= thr_hyb_maxSSS)
fut_cat_CLEAN        <- cat_map_from(bin_corr_fut_nomask, bin_hyb_fut_nomask)
fut_cat_CLEAN        <- terra::ifel(is.na(fut_cat_CLEAN), NA, round(fut_cat_CLEAN))

bin_corr_pres_nomask <- terra::ifel(is.na(pred_corr_curr_mean), NA, pred_corr_curr_mean >= thr_corr_maxSSS)
bin_hyb_pres_nomask  <- terra::ifel(is.na(pred_hyb_curr_mean),  NA, pred_hyb_curr_mean  >= thr_hyb_maxSSS)
pres_cat_CLEAN       <- cat_map_from(bin_corr_pres_nomask, bin_hyb_pres_nomask)
pres_cat_CLEAN       <- terra::ifel(is.na(pres_cat_CLEAN), NA, round(pres_cat_CLEAN))

# ============================================================================ #
# EXPORT GeoTIFF
# ============================================================================ #

cat("\n--- GeoTIFF exports ---\n")

# Categorical (INT1U, nodata=255)
.write_gtiff(pres_cat, "Fig5_categorical_present.tif", dtype = "INT1U", nodata = 255)
.write_gtiff(fut_cat,  "Fig5_categorical_future_OVERLAY.tif", dtype = "INT1U", nodata = 255)
if (!is.null(fut_cat_MASKED)) {
  .write_gtiff(fut_cat_MASKED, "Fig5_categorical_future_MASKED.tif", dtype = "INT1U", nodata = 255)
}

# Continuous (FLT4S)
.write_gtiff(diff_present_cont, "DeltaSuitability_present_HminusC.tif")
.write_gtiff(diff_future_cont,  "DeltaSuitability_future_OVERLAY_HminusC.tif")
if (!is.null(diff_fut_MASKED)) {
  .write_gtiff(diff_fut_MASKED, "DeltaSuitability_future_MASKED_HminusC.tif")
}

# ============================================================================ #
# EXPORT ASCII (.asc.gz)
# ============================================================================ #

cat("\n--- ASCII exports (with gzip) ---\n")

# Categorical
.write_asc(pres_cat, "Fig5_categorical_present")
.write_asc(fut_cat,  "Fig5_categorical_future_OVERLAY")
if (!is.null(fut_cat_MASKED)) {
  .write_asc(fut_cat_MASKED, "Fig5_categorical_future_MASKED")
}

# Continuous
.write_asc(diff_present_cont, "DeltaSuitability_present_HminusC")
.write_asc(diff_future_cont,  "DeltaSuitability_future_OVERLAY_HminusC")
if (!is.null(diff_fut_MASKED)) {
  .write_asc(diff_fut_MASKED, "DeltaSuitability_future_MASKED_HminusC")
}

# Clean versions (no MESS)
cat("\n--- Clean ASCII exports (no MESS) ---\n")
.write_asc(pres_cat_CLEAN, "Fig5_categorical_present_CLEAN", gzip_after = FALSE)
.write_asc(fut_cat_CLEAN,  "Fig5_categorical_future_CLEAN",  gzip_after = FALSE)
.write_asc(diff_present_cont,  "DeltaSuitability_present_HminusC_CLEAN",  gzip_after = FALSE)
.write_asc(diff_future_cont,   "DeltaSuitability_future_HminusC_CLEAN",   gzip_after = FALSE)

# ============================================================================ #
# README + MANIFEST
# ============================================================================ #

cat("\n--- Writing manifest and README ---\n")

legend_rows <- c(
  "Fig5_categorical_present.tif|Categorical agreement (present): 1=Correlative only, 2=Hybrid only, 3=Both, NA=Neither",
  "Fig5_categorical_future_OVERLAY.tif|Categorical agreement (future, MESS overlay)",
  if (!is.null(fut_cat_MASKED)) "Fig5_categorical_future_MASKED.tif|Categorical agreement (future, MESS masked >= threshold)" else NULL,
  "DeltaSuitability_present_HminusC.tif|Continuous delta suitability present (Hybrid - Correlative)",
  "DeltaSuitability_future_OVERLAY_HminusC.tif|Continuous delta suitability future (Hybrid - Correlative, MESS overlay)",
  if (!is.null(diff_fut_MASKED)) "DeltaSuitability_future_MASKED_HminusC.tif|Continuous delta suitability future (Hybrid - Correlative, MESS masked >= threshold)" else NULL
)
legend_df <- tidyr::separate(
  tibble::tibble(x = legend_rows), x,
  into = c("raster_file", "description"), sep = "\\|", remove = TRUE
)
readr::write_csv(legend_df, file.path(EXP_DIR, "Fig5_export_manifest.csv"))

readme_txt <- paste0(
  "# Fig 5 - Export bundle\n",
  "- Thresholds (maxSSS): correlative = ", signif(thr_corr_maxSSS, 6),
  " ; hybrid = ", signif(thr_hyb_maxSSS, 6), "\n",
  "- Categories: 1=Correlative only, 2=Hybrid only, 3=Both, NA=Neither\n",
  "- Continuous rasters: delta = Hybrid - Correlative (unit: suitability [0-1])\n",
  if (isTRUE(has_mess)) paste0("- MESS overlay/masked: threshold = ", MESS_THR, "\n") else "",
  "- CRS: ", terra::crs(pres_cat), "\n",
  "- GeoTIFF: LZW compressed\n",
  "- ASCII: .asc.gz (NAflag = -9999); clean versions uncompressed\n"
)
writeLines(readme_txt, con = file.path(EXP_DIR, "README_exports.txt"))

readme_asc <- paste0(
  "# Fig 5 - ASCII Exports\n",
  "- Thresholds (maxSSS): correlative = ", signif(thr_corr_maxSSS, 6),
  " ; hybrid = ", signif(thr_hyb_maxSSS, 6), "\n",
  "- Categories: 1=Correlative only, 2=Hybrid only, 3=Both, NA=Neither\n",
  "- Continuous rasters: delta = Hybrid - Correlative (unit: suitability [0-1])\n",
  if (isTRUE(has_mess)) paste0("- MESS threshold = ", MESS_THR, "\n") else "",
  "- CRS: ", terra::crs(pres_cat), "\n",
  "- NAflag: -9999\n",
  "\n# Clean versions (no MESS overlay/mask):\n",
  "- Fig5_categorical_present_CLEAN.asc\n",
  "- Fig5_categorical_future_CLEAN.asc\n",
  "- DeltaSuitability_present_HminusC_CLEAN.asc\n",
  "- DeltaSuitability_future_HminusC_CLEAN.asc\n"
)
writeLines(readme_asc, con = file.path(ASC_DIR, "README_ASCII.txt"))

# ============================================================================ #
# AREA STATISTICS BY CLASS
# ============================================================================ #

cat("\n--- Computing area statistics ---\n")

ref_for_area <- if (grepl("longlat", terra::crs(pres_cat), ignore.case = TRUE)) {
  terra::project(pres_cat, "EPSG:3035", method = "near")
} else pres_cat
cell_area_km2 <- prod(terra::res(ref_for_area)) / 1e6

area_tab <- function(cat_r, scenario_label) {
  rc <- if (terra::crs(cat_r) != terra::crs(ref_for_area)) {
    terra::project(cat_r, terra::crs(ref_for_area), method = "near")
  } else cat_r
  vals <- terra::values(rc, mat = FALSE)
  tibble::tibble(class = vals) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::count(class, name = "n_cells") %>%
    dplyr::mutate(
      area_km2 = n_cells * cell_area_km2,
      label    = dplyr::recode(as.integer(class),
                               `1` = "Correlative only", `2` = "Hybrid only", `3` = "Both"),
      scenario = scenario_label
    )
}

areas_present <- area_tab(pres_cat, "Present")
areas_future  <- area_tab(fut_cat,  "Future (overlay)")
readr::write_csv(
  dplyr::bind_rows(areas_present, areas_future),
  file.path(EXP_DIR, "Fig5_areas_by_class_km2.csv")
)

# Delta suitability stats
summ_delta <- function(r) {
  v <- terra::values(terra::mask(r, land_mask), mat = FALSE)
  v <- v[is.finite(v)]
  tibble::tibble(
    mean   = mean(v), median = stats::median(v),
    sd     = stats::sd(v),
    q05    = stats::quantile(v, 0.05),
    q95    = stats::quantile(v, 0.95),
    n      = length(v)
  )
}

delta_stats <- dplyr::bind_rows(
  summ_delta(diff_present_cont) %>% dplyr::mutate(scenario = "Present"),
  summ_delta(diff_future_cont)  %>% dplyr::mutate(scenario = "Future (overlay)"),
  if (!is.null(diff_fut_MASKED))
    summ_delta(diff_fut_MASKED) %>% dplyr::mutate(scenario = "Future (masked)") else NULL
)
readr::write_csv(delta_stats, file.path(EXP_DIR, "DeltaSuitability_stats.csv"))

cat("\n  Section 3 (raster exports) complete\n")
cat("  GeoTIFF exports:", EXP_DIR, "\n")
cat("  ASCII exports:  ", ASC_DIR, "\n")


# ============================================================================ #
# Reproducibility footer                                                       #
# ============================================================================ #

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "session_info_make_figures.txt"))

cat("\n  OK: RCode_SDM_make_figures.R complete\n")
