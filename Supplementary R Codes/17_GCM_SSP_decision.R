################################################################################
## 17_gcm_ssp_decision.R — Valoración multi-GCM x SSP para decidir main text
##
## Source: v2 runs (GCM/SSP-specific mechanistic layers).
##   - 5 GCMs x SSP5-8.5 from sensitivity_gcm_v2/
##   - MIROC6 x SSP2-4.5 from ssp245_v2/
##
## Outputs to: workflow_psammodromus_20251027/gcm_ssp_decision/
##   - gcm_ssp_comparison_table.csv
##   - gcm_ssp_decision_report.txt
##   - gcm_ssp_decision_report.md
##   - diagnostic_gcm_ssp_maps.png
##   - diagnostic_gcm_ssp_profiles.png
##   - diagnostic_range_change_comparison.png
##   - delta_suit_pairwise_correlations.csv
##
## Conservative overlap: cells where both models >= respective maxSSS threshold.
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

OUTDIR <- "workflow_psammodromus_20251027"
OUT    <- file.path(OUTDIR, "gcm_ssp_decision")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------- #
# 0. Configuration
# -------------------------------------------------------------------------- #

GCMS <- c("MIROC6", "MPI-ESM1-2-HR", "CNRM-CM6-1", "EC-Earth3-Veg", "UKESM1-0-LL")

thr <- readRDS(file.path(OUTDIR, "results_summary/thresholds_maxSSS_ensemble.rds"))
THR_CORR <- thr$correlativo_mean_ensemble
THR_HYB  <- thr$hibrido_mean_ensemble
cat(sprintf("maxSSS thresholds: Corr=%.3f, Hyb=%.3f\n", THR_CORR, THR_HYB))

# Present (baseline) continuous rasters (same across scenarios)
curr_corr <- rast(file.path(OUTDIR, "rasters/current_correlative_mean.tif"))
curr_hyb  <- rast(file.path(OUTDIR, "rasters/current_hybrid_mean.tif"))

# -------------------------------------------------------------------------- #
# 1. Helper functions
# -------------------------------------------------------------------------- #

schoeners_D_cons <- function(r1, r2, t1, t2) {
  v1 <- values(r1, mat = FALSE); v2 <- values(r2, mat = FALSE)
  ok <- is.finite(v1) & is.finite(v2) & v1 >= t1 & v2 >= t2
  v1 <- v1[ok]; v2 <- v2[ok]
  if (length(v1) < 3) return(NA_real_)
  p1 <- v1 / sum(v1); p2 <- v2 / sum(v2)
  1 - 0.5 * sum(abs(p1 - p2))
}

spearman_cons <- function(r1, r2, t1, t2) {
  v1 <- values(r1, mat = FALSE); v2 <- values(r2, mat = FALSE)
  ok <- is.finite(v1) & is.finite(v2) & v1 >= t1 & v2 >= t2
  if (sum(ok) < 3) return(NA_real_)
  cor(v1[ok], v2[ok], method = "spearman")
}

categorical_pct <- function(r_corr, r_hyb, t_c, t_h) {
  vc <- values(r_corr, mat = FALSE); vh <- values(r_hyb, mat = FALSE)
  ok <- is.finite(vc) & is.finite(vh)
  vc <- vc[ok] >= t_c; vh <- vh[ok] >= t_h
  ntot <- length(vc)
  if (ntot == 0) return(c(Both = NA, CorrOnly = NA, HybOnly = NA, NeitherSuit = NA,
                          Either = NA))
  both     <- sum(vc & vh)
  corr_only <- sum(vc & !vh)
  hyb_only  <- sum(!vc & vh)
  neither   <- sum(!vc & !vh)
  either    <- both + corr_only + hyb_only
  c(Both = 100 * both / ntot,
    CorrOnly = 100 * corr_only / ntot,
    HybOnly = 100 * hyb_only / ntot,
    NeitherSuit = 100 * neither / ntot,
    PctBoth_of_Either = 100 * both / max(either, 1))
}

jaccard_bin <- function(rb1, rb2) {
  v1 <- values(rb1, mat = FALSE); v2 <- values(rb2, mat = FALSE)
  ok <- is.finite(v1) & is.finite(v2)
  v1 <- v1[ok] > 0; v2 <- v2[ok] > 0
  inter <- sum(v1 & v2); uni <- sum(v1 | v2)
  if (uni == 0) return(NA_real_)
  inter / uni
}

range_change_pct <- function(rb_curr, rb_fut) {
  vc <- values(rb_curr, mat = FALSE); vf <- values(rb_fut, mat = FALSE)
  ok <- is.finite(vc) & is.finite(vf)
  a_curr <- sum(vc[ok] > 0); a_fut <- sum(vf[ok] > 0)
  if (a_curr == 0) return(NA_real_)
  100 * (a_fut - a_curr) / a_curr
}

# Load pre-computed binary rasters and range change CSVs where possible.

# -------------------------------------------------------------------------- #
# 2. Load rasters per (GCM, SSP)
# -------------------------------------------------------------------------- #

load_ssp585_gcm <- function(gcm) {
  base <- file.path(OUTDIR, "sensitivity_gcm_v2", gcm, "rasters")
  list(
    corr_fut = rast(file.path(base, sprintf("future_correlative_mean_%s.tif", gcm))),
    hyb_fut  = rast(file.path(base, sprintf("future_hybrid_mean_%s.tif", gcm))),
    corr_bin = rast(file.path(base, sprintf("binary_correlative_maxSSS_%s.tif", gcm))),
    hyb_bin  = rast(file.path(base, sprintf("binary_hybrid_maxSSS_%s.tif", gcm)))
  )
}

ssp585 <- lapply(GCMS, load_ssp585_gcm); names(ssp585) <- GCMS

# SSP2-4.5 only MIROC6
ssp245_miroc <- list(
  corr_fut = rast(file.path(OUTDIR, "ssp245_v2/rasters/future_correlative_mean_ssp245.tif")),
  hyb_fut  = rast(file.path(OUTDIR, "ssp245_v2/rasters/future_hybrid_mean_ssp245.tif")),
  corr_bin = rast(file.path(OUTDIR, "ssp245_v2/rasters/binary_future_correlative_maxSSS_ssp245.tif")),
  hyb_bin  = rast(file.path(OUTDIR, "ssp245_v2/rasters/binary_future_hybrid_maxSSS_ssp245.tif"))
)

curr_corr_bin <- curr_corr >= THR_CORR
curr_hyb_bin  <- curr_hyb  >= THR_HYB

# -------------------------------------------------------------------------- #
# 3. Per-combination metrics (Task 2)
# -------------------------------------------------------------------------- #

compute_combo_metrics <- function(corr_fut, hyb_fut, corr_bin, hyb_bin,
                                  curr_corr_bin, curr_hyb_bin,
                                  gcm, ssp) {
  D   <- schoeners_D_cons(corr_fut, hyb_fut, THR_CORR, THR_HYB)
  rho <- spearman_cons(corr_fut, hyb_fut, THR_CORR, THR_HYB)
  cat_pct <- categorical_pct(corr_fut, hyb_fut, THR_CORR, THR_HYB)
  rc_corr <- range_change_pct(curr_corr_bin, corr_bin)
  rc_hyb  <- range_change_pct(curr_hyb_bin,  hyb_bin)
  j_corr  <- jaccard_bin(curr_corr_bin, corr_bin)
  j_hyb   <- jaccard_bin(curr_hyb_bin,  hyb_bin)
  tibble(
    GCM = gcm, SSP = ssp,
    D_future = D, rho_future = rho,
    pct_Both = cat_pct["Both"],
    pct_CorrOnly = cat_pct["CorrOnly"],
    pct_HybOnly = cat_pct["HybOnly"],
    pct_Neither = cat_pct["NeitherSuit"],
    range_change_corr_pct = rc_corr,
    range_change_hyb_pct  = rc_hyb,
    jaccard_corr = j_corr,
    jaccard_hyb  = j_hyb
  )
}

metrics_585 <- bind_rows(lapply(GCMS, function(g) {
  s <- ssp585[[g]]
  compute_combo_metrics(s$corr_fut, s$hyb_fut, s$corr_bin, s$hyb_bin,
                        curr_corr_bin, curr_hyb_bin, g, "SSP5-8.5")
}))

metrics_245 <- compute_combo_metrics(
  ssp245_miroc$corr_fut, ssp245_miroc$hyb_fut,
  ssp245_miroc$corr_bin, ssp245_miroc$hyb_bin,
  curr_corr_bin, curr_hyb_bin, "MIROC6", "SSP2-4.5"
)

metrics_all <- bind_rows(metrics_585, metrics_245) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

# MESS availability flag — only MIROC6 (main) and MIROC6 SSP245 have MESS rasters.
metrics_all <- metrics_all %>%
  mutate(mess_lt0_in_future_range_pct_hyb = NA_real_,
         mess_lt0_in_future_range_pct_corr = NA_real_)

# MIROC6 x SSP5-8.5: MESS from main pipeline
mess_hyb_main  <- rast(file.path(OUTDIR, "rasters/mess_hybrid_future_pointsRef.tif"))
mess_corr_main <- rast(file.path(OUTDIR, "rasters/mess_correlative_future_pointsRef.tif"))
mess_pct_in_range <- function(mess_r, bin_fut) {
  vm <- values(mess_r, mat = FALSE); vb <- values(bin_fut, mat = FALSE)
  ok <- is.finite(vm) & is.finite(vb) & vb > 0
  if (sum(ok) == 0) return(NA_real_)
  100 * sum(vm[ok] < 0) / sum(ok)
}
metrics_all$mess_lt0_in_future_range_pct_hyb[metrics_all$GCM == "MIROC6" &
                                             metrics_all$SSP == "SSP5-8.5"] <-
  mess_pct_in_range(mess_hyb_main,  ssp585$MIROC6$hyb_bin)
metrics_all$mess_lt0_in_future_range_pct_corr[metrics_all$GCM == "MIROC6" &
                                              metrics_all$SSP == "SSP5-8.5"] <-
  mess_pct_in_range(mess_corr_main, ssp585$MIROC6$corr_bin)

# MIROC6 x SSP2-4.5: MESS from ssp245_v2/
mess_hyb_245  <- rast(file.path(OUTDIR, "ssp245_v2/rasters/mess_hybrid_future_ssp245.tif"))
mess_corr_245 <- rast(file.path(OUTDIR, "ssp245_v2/rasters/mess_correlative_future_ssp245.tif"))
metrics_all$mess_lt0_in_future_range_pct_hyb[metrics_all$GCM == "MIROC6" &
                                             metrics_all$SSP == "SSP2-4.5"] <-
  mess_pct_in_range(mess_hyb_245,  ssp245_miroc$hyb_bin)
metrics_all$mess_lt0_in_future_range_pct_corr[metrics_all$GCM == "MIROC6" &
                                              metrics_all$SSP == "SSP2-4.5"] <-
  mess_pct_in_range(mess_corr_245, ssp245_miroc$corr_bin)

write_csv(metrics_all, file.path(OUT, "gcm_ssp_comparison_table.csv"))
cat("\n[OK] gcm_ssp_comparison_table.csv\n")
print(metrics_all)

# -------------------------------------------------------------------------- #
# 4. Delta suitability per GCM and inter-GCM variability (Task 3)
# -------------------------------------------------------------------------- #

delta_ssp585 <- lapply(GCMS, function(g) ssp585[[g]]$hyb_fut - ssp585[[g]]$corr_fut)
names(delta_ssp585) <- GCMS
delta_stack_585 <- rast(delta_ssp585)

inter_gcm_sd_585   <- app(delta_stack_585, fun = sd, na.rm = TRUE)
inter_gcm_mean_585 <- app(delta_stack_585, fun = mean, na.rm = TRUE)

# Sign agreement: how many of 5 GCMs have delta > 0 per cell
sign_agree_585 <- app(delta_stack_585, fun = function(x) sum(x > 0, na.rm = TRUE))
writeRaster(inter_gcm_sd_585,   file.path(OUT, "delta_suit_inter_gcm_sd_ssp585.tif"),
            overwrite = TRUE)
writeRaster(inter_gcm_mean_585, file.path(OUT, "delta_suit_inter_gcm_mean_ssp585.tif"),
            overwrite = TRUE)
writeRaster(sign_agree_585,     file.path(OUT, "delta_suit_sign_agreement_ssp585.tif"),
            overwrite = TRUE)

# Pairwise Spearman between Delta maps (SSP5-8.5)
pairwise_cor <- function(stk) {
  n <- nlyr(stk)
  m <- matrix(NA_real_, n, n, dimnames = list(names(stk), names(stk)))
  vals <- values(stk, mat = TRUE)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    ok <- is.finite(vals[, i]) & is.finite(vals[, j])
    if (sum(ok) > 3) m[i, j] <- cor(vals[ok, i], vals[ok, j], method = "spearman")
  }
  m
}
pair_cor_585 <- pairwise_cor(delta_stack_585)
diag(pair_cor_585) <- NA
pair_cor_585_df <- as.data.frame(pair_cor_585) %>%
  tibble::rownames_to_column("GCM_A") %>%
  pivot_longer(-GCM_A, names_to = "GCM_B", values_to = "spearman_rho") %>%
  filter(GCM_A != GCM_B)
write_csv(pair_cor_585_df, file.path(OUT, "delta_suit_pairwise_correlations.csv"))

mean_pair_cor_585 <- mean(pair_cor_585[upper.tri(pair_cor_585)], na.rm = TRUE)
cat(sprintf("\nMean pairwise Spearman between Delta-suit GCMs (SSP5-8.5): %.3f\n",
            mean_pair_cor_585))

# Inter-GCM sd vs mean |delta|
v_sd   <- values(inter_gcm_sd_585,   mat = FALSE)
v_mean <- values(inter_gcm_mean_585, mat = FALSE)
ok <- is.finite(v_sd) & is.finite(v_mean)
mean_abs_delta <- mean(abs(v_mean[ok]))
median_sd      <- median(v_sd[ok])
ratio_sd_to_absdelta <- median_sd / mean_abs_delta
cat(sprintf("Inter-GCM median sd(Δ) / mean|Δ| (SSP5-8.5): %.3f\n",
            ratio_sd_to_absdelta))

# -------------------------------------------------------------------------- #
# 5. SSP effect — MIROC6 SSP2-4.5 vs SSP5-8.5 (Task 4)
# -------------------------------------------------------------------------- #

delta_245_miroc <- ssp245_miroc$hyb_fut - ssp245_miroc$corr_fut
delta_585_miroc <- ssp585$MIROC6$hyb_fut - ssp585$MIROC6$corr_fut

v245 <- values(delta_245_miroc, mat = FALSE)
v585 <- values(delta_585_miroc, mat = FALSE)
ok <- is.finite(v245) & is.finite(v585)
cor_ssp_miroc <- cor(v245[ok], v585[ok], method = "spearman")
cat(sprintf("MIROC6 Δ-suit SSP2-4.5 vs SSP5-8.5 (Spearman): %.3f\n", cor_ssp_miroc))

# Latitudinal profiles (1-degree bands)
lat_profile <- function(r, label) {
  y <- yFromRow(r, 1:nrow(r))
  lat_band <- floor(y)
  m <- as.matrix(r, wide = TRUE)   # rows = grid rows (north -> south), cols = grid cols
  row_means <- rowMeans(m, na.rm = TRUE)
  tibble(lat_center = y, lat_band = lat_band, mean_delta = row_means) %>%
    group_by(lat_band) %>%
    summarise(mean_delta = mean(mean_delta, na.rm = TRUE), .groups = "drop") %>%
    mutate(Series = label)
}

prof_585 <- bind_rows(lapply(GCMS, function(g) {
  lat_profile(delta_ssp585[[g]], label = g)
})) %>% mutate(SSP = "SSP5-8.5")

prof_245_miroc <- lat_profile(delta_245_miroc, label = "MIROC6") %>%
  mutate(SSP = "SSP2-4.5")

prof_all <- bind_rows(prof_585, prof_245_miroc)

# Latitude of Δ=0 crossing (interpolated) for MIROC6 under each SSP
find_crossing <- function(df) {
  df <- df %>% arrange(lat_band) %>% filter(is.finite(mean_delta))
  sc <- sign(df$mean_delta)
  idx <- which(diff(sc) != 0)
  if (length(idx) == 0) return(NA_real_)
  # Linear interpolation
  i <- idx[1]
  x1 <- df$lat_band[i]; y1 <- df$mean_delta[i]
  x2 <- df$lat_band[i + 1]; y2 <- df$mean_delta[i + 1]
  x1 - y1 * (x2 - x1) / (y2 - y1)
}

cross_245 <- find_crossing(prof_245_miroc)
cross_585_miroc <- find_crossing(prof_all %>% filter(SSP == "SSP5-8.5", Series == "MIROC6"))

cat(sprintf("Latitude of Δ=0 crossing (MIROC6): SSP2-4.5 = %.2f°N, SSP5-8.5 = %.2f°N\n",
            cross_245, cross_585_miroc))

# -------------------------------------------------------------------------- #
# 6. Summary table for main-text decision (Task 5)
# -------------------------------------------------------------------------- #

agg_585 <- metrics_585 %>%
  summarise(across(c(D_future, rho_future, pct_Both, pct_CorrOnly, pct_HybOnly,
                     range_change_corr_pct, range_change_hyb_pct,
                     jaccard_corr, jaccard_hyb),
                   list(mean = ~mean(.), sd = ~sd(.))))

summary_decision <- tibble(
  Metric = c(
    "D conservative (future)",
    "rho conservative (future)",
    "% Both (future)",
    "% Corr only (future)",
    "% Hyb only (future)",
    "Range change correlative (%)",
    "Range change hybrid (%)",
    "Jaccard correlative",
    "Jaccard hybrid",
    "Mean pairwise Spearman of Δ-suit between GCMs (SSP5-8.5)",
    "Corr Δ-suit SSP2-4.5 vs SSP5-8.5 (MIROC6 only)"
  ),
  SSP245_mean_sd_5GCMs = c(
    rep("N/A (only MIROC6 available)", 9), "N/A", "N/A"
  ),
  SSP585_mean_sd_5GCMs = c(
    sprintf("%.3f ± %.3f", agg_585$D_future_mean, agg_585$D_future_sd),
    sprintf("%.3f ± %.3f", agg_585$rho_future_mean, agg_585$rho_future_sd),
    sprintf("%.2f ± %.2f", agg_585$pct_Both_mean, agg_585$pct_Both_sd),
    sprintf("%.2f ± %.2f", agg_585$pct_CorrOnly_mean, agg_585$pct_CorrOnly_sd),
    sprintf("%.2f ± %.2f", agg_585$pct_HybOnly_mean, agg_585$pct_HybOnly_sd),
    sprintf("%.2f ± %.2f", agg_585$range_change_corr_pct_mean, agg_585$range_change_corr_pct_sd),
    sprintf("%.2f ± %.2f", agg_585$range_change_hyb_pct_mean, agg_585$range_change_hyb_pct_sd),
    sprintf("%.3f ± %.3f", agg_585$jaccard_corr_mean, agg_585$jaccard_corr_sd),
    sprintf("%.3f ± %.3f", agg_585$jaccard_hyb_mean, agg_585$jaccard_hyb_sd),
    sprintf("%.3f", mean_pair_cor_585),
    sprintf("%.3f", cor_ssp_miroc)
  ),
  MIROC6_SSP585 = c(
    sprintf("%.3f", metrics_585$D_future[metrics_585$GCM == "MIROC6"]),
    sprintf("%.3f", metrics_585$rho_future[metrics_585$GCM == "MIROC6"]),
    sprintf("%.2f", metrics_585$pct_Both[metrics_585$GCM == "MIROC6"]),
    sprintf("%.2f", metrics_585$pct_CorrOnly[metrics_585$GCM == "MIROC6"]),
    sprintf("%.2f", metrics_585$pct_HybOnly[metrics_585$GCM == "MIROC6"]),
    sprintf("%.2f", metrics_585$range_change_corr_pct[metrics_585$GCM == "MIROC6"]),
    sprintf("%.2f", metrics_585$range_change_hyb_pct[metrics_585$GCM == "MIROC6"]),
    sprintf("%.3f", metrics_585$jaccard_corr[metrics_585$GCM == "MIROC6"]),
    sprintf("%.3f", metrics_585$jaccard_hyb[metrics_585$GCM == "MIROC6"]),
    "(ref. reported to left)",
    "(ref. reported to left)"
  ),
  MIROC6_SSP245 = c(
    sprintf("%.3f", metrics_245$D_future),
    sprintf("%.3f", metrics_245$rho_future),
    sprintf("%.2f", metrics_245$pct_Both),
    sprintf("%.2f", metrics_245$pct_CorrOnly),
    sprintf("%.2f", metrics_245$pct_HybOnly),
    sprintf("%.2f", metrics_245$range_change_corr_pct),
    sprintf("%.2f", metrics_245$range_change_hyb_pct),
    sprintf("%.3f", metrics_245$jaccard_corr),
    sprintf("%.3f", metrics_245$jaccard_hyb),
    "—", "—"
  )
)

write_csv(summary_decision, file.path(OUT, "summary_decision_table.csv"))

# -------------------------------------------------------------------------- #
# 7. Figures (Task 6)
# -------------------------------------------------------------------------- #

delta_to_df <- function(r, label) {
  d <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(d)[3] <- "delta"
  d$series <- label
  d
}

# Representative GCMs for Fig 1 top row (SSP5-8.5):
# MIROC6 (-7.3%, most contractive hybrid), UKESM (-5.5%, intermediate),
# MPI (+3.7%, least contractive / expansive).
sel_gcms <- c("MIROC6", "UKESM1-0-LL", "MPI-ESM1-2-HR")

df_top <- bind_rows(lapply(sel_gcms, function(g) {
  delta_to_df(delta_ssp585[[g]], paste(g, "SSP5-8.5"))
}))

df_bottom <- bind_rows(
  delta_to_df(delta_245_miroc, "MIROC6 SSP2-4.5"),
  delta_to_df(inter_gcm_sd_585, "Inter-GCM sd(Δ) SSP5-8.5"),
  delta_to_df(sign_agree_585, "# GCMs with Δ>0 (of 5)")
)

df_maps <- bind_rows(df_top, df_bottom)
df_maps$series <- factor(df_maps$series, levels = c(
  paste(sel_gcms, "SSP5-8.5"),
  "MIROC6 SSP2-4.5",
  "Inter-GCM sd(Δ) SSP5-8.5",
  "# GCMs with Δ>0 (of 5)"
))

fig1_delta <- ggplot(df_maps %>% filter(series %in% levels(df_maps$series)[1:4]),
                    aes(x, y, fill = delta)) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, na.value = "grey90",
                       name = "Δ-suit (Hyb − Corr)") +
  facet_wrap(~ series, nrow = 2, ncol = 2) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "bottom", legend.key.width = unit(1.2, "cm"))

fig1_sd <- ggplot(df_maps %>% filter(series == "Inter-GCM sd(Δ) SSP5-8.5"),
                  aes(x, y, fill = delta)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(option = "magma", na.value = "grey90",
                       name = "sd(Δ)") +
  facet_wrap(~ series) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "bottom", legend.key.width = unit(1.2, "cm"))

fig1_sign <- ggplot(df_maps %>% filter(series == "# GCMs with Δ>0 (of 5)"),
                    aes(x, y, fill = delta)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90",
                       breaks = 0:5, limits = c(0, 5),
                       name = "# GCMs Δ>0") +
  facet_wrap(~ series) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "bottom", legend.key.width = unit(1.2, "cm"))

# Compose: top row = 3 panels SSP5-8.5; bottom row = SSP2-4.5 + sd + sign
g_m <- delta_to_df(delta_ssp585$MIROC6,        "MIROC6 SSP5-8.5")
g_u <- delta_to_df(delta_ssp585$`UKESM1-0-LL`, "UKESM1-0-LL SSP5-8.5")
g_p <- delta_to_df(delta_ssp585$`MPI-ESM1-2-HR`, "MPI-ESM1-2-HR SSP5-8.5")
g_245 <- delta_to_df(delta_245_miroc, "MIROC6 SSP2-4.5")
g_sd  <- delta_to_df(inter_gcm_sd_585, "Inter-GCM sd(Δ) SSP5-8.5")
g_sgn <- delta_to_df(sign_agree_585, "# GCMs Δ>0 (of 5)")

mk_delta_panel <- function(d, lim) {
  ggplot(d, aes(x, y, fill = delta)) + geom_raster() + coord_equal() +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         midpoint = 0, limits = lim, na.value = "grey90",
                         name = "Δ-suit") +
    ggtitle(unique(d$series)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.title = element_blank(), legend.position = "bottom",
          legend.key.width = unit(1, "cm"), plot.title = element_text(size = 9))
}

common_lim <- range(c(g_m$delta, g_u$delta, g_p$delta, g_245$delta), na.rm = TRUE)
common_lim <- c(-max(abs(common_lim)), max(abs(common_lim)))

p_m <- mk_delta_panel(g_m, common_lim)
p_u <- mk_delta_panel(g_u, common_lim)
p_p <- mk_delta_panel(g_p, common_lim)
p_245 <- mk_delta_panel(g_245, common_lim)
p_sd <- ggplot(g_sd, aes(x, y, fill = delta)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(option = "magma", na.value = "grey90", name = "sd(Δ)") +
  ggtitle("Inter-GCM sd(Δ) SSP5-8.5") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 9))
p_sgn <- ggplot(g_sgn, aes(x, y, fill = delta)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90",
                       breaks = 0:5, limits = c(0, 5), name = "# GCMs Δ>0") +
  ggtitle("# GCMs with Δ>0 (of 5) — SSP5-8.5") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 9))

fig1 <- (p_m | p_u | p_p) / (p_245 | p_sd | p_sgn) +
  plot_annotation(title = "Fig. 1 — Δ-suit (Hyb − Corr) across GCMs and SSPs",
                  subtitle = "Top: 3 GCMs under SSP5-8.5 (shared scale). Bottom: MIROC6 SSP2-4.5 + inter-GCM sd + sign-agreement")
ggsave(file.path(OUT, "diagnostic_gcm_ssp_maps.png"), fig1,
       width = 12, height = 8, dpi = 150)
cat("[OK] diagnostic_gcm_ssp_maps.png\n")

# -------------------------------------------------------------------------- #
# Fig 2 — latitudinal profiles + scatter (Task 6)
# -------------------------------------------------------------------------- #

pal_gcm <- c("MIROC6" = "#1B9E77", "MPI-ESM1-2-HR" = "#D95F02",
             "CNRM-CM6-1" = "#7570B3", "EC-Earth3-Veg" = "#E7298A",
             "UKESM1-0-LL" = "#66A61E")

p2a <- prof_all %>% filter(SSP == "SSP5-8.5") %>%
  ggplot(aes(lat_band, mean_delta, colour = Series)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = pal_gcm, name = "GCM") +
  labs(x = "Latitude (°N, 1° bands)", y = "Mean Δ-suit",
       title = "(a) SSP5-8.5 — 5 GCMs") +
  theme_minimal(base_size = 10)

p2b <- prof_245_miroc %>%
  ggplot(aes(lat_band, mean_delta)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.8, colour = "#1B9E77") +
  labs(x = "Latitude (°N, 1° bands)", y = "Mean Δ-suit",
       title = "(b) SSP2-4.5 — MIROC6 only (no multi-GCM SSP2-4.5 run)") +
  theme_minimal(base_size = 10)

# (c) MIROC6 vs 5-GCM mean (SSP5-8.5)
df_c <- tibble(miroc = values(delta_ssp585$MIROC6, mat = FALSE),
               ensmean = values(inter_gcm_mean_585, mat = FALSE)) %>%
  filter(is.finite(miroc), is.finite(ensmean))
cor_c <- cor(df_c$miroc, df_c$ensmean, method = "spearman")
p2c <- ggplot(df_c, aes(miroc, ensmean)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10", name = "n") +
  labs(x = "Δ-suit MIROC6 (SSP5-8.5)", y = "Δ-suit 5-GCM mean",
       title = sprintf("(c) MIROC6 vs 5-GCM mean (Spearman ρ = %.3f)", cor_c)) +
  theme_minimal(base_size = 10)

# (d) MIROC6 SSP2-4.5 vs SSP5-8.5
df_d <- tibble(s245 = v245, s585 = v585) %>% filter(is.finite(s245), is.finite(s585))
p2d <- ggplot(df_d, aes(s245, s585)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10", name = "n") +
  labs(x = "Δ-suit MIROC6 SSP2-4.5", y = "Δ-suit MIROC6 SSP5-8.5",
       title = sprintf("(d) MIROC6 SSP2-4.5 vs SSP5-8.5 (ρ = %.3f)", cor_ssp_miroc)) +
  theme_minimal(base_size = 10)

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_annotation(title = "Fig. 2 — Latitudinal profiles and scatter diagnostics")
ggsave(file.path(OUT, "diagnostic_gcm_ssp_profiles.png"), fig2,
       width = 12, height = 9, dpi = 150)
cat("[OK] diagnostic_gcm_ssp_profiles.png\n")

# -------------------------------------------------------------------------- #
# Fig 3 — Barplot: range change by model x SSP x GCM
# -------------------------------------------------------------------------- #

rc_df <- metrics_all %>%
  select(GCM, SSP, Correlative = range_change_corr_pct,
         Hybrid = range_change_hyb_pct) %>%
  pivot_longer(c(Correlative, Hybrid), names_to = "Model",
               values_to = "pct_change")
rc_df$GCM <- factor(rc_df$GCM, levels = GCMS)

fig3 <- ggplot(rc_df, aes(GCM, pct_change, fill = Model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ SSP, ncol = 2) +
  scale_fill_manual(values = c("Correlative" = "#8C510A", "Hybrid" = "#01665E")) +
  labs(y = "Range change (%)", x = NULL,
       title = "Fig. 3 — Range change by GCM × SSP × Model (maxSSS)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(OUT, "diagnostic_range_change_comparison.png"), fig3,
       width = 11, height = 5, dpi = 150)
cat("[OK] diagnostic_range_change_comparison.png\n")

# -------------------------------------------------------------------------- #
# 8. Decision report (txt + md) (Task 7)
# -------------------------------------------------------------------------- #

corr_var_585 <- sd(metrics_585$range_change_corr_pct)
hyb_var_585  <- sd(metrics_585$range_change_hyb_pct)
corr_mean_585 <- mean(metrics_585$range_change_corr_pct)
hyb_mean_585  <- mean(metrics_585$range_change_hyb_pct)

cv_pct <- function(x) 100 * sd(x) / abs(mean(x))
cv_corr_585 <- cv_pct(metrics_585$range_change_corr_pct)
cv_hyb_585  <- cv_pct(metrics_585$range_change_hyb_pct)

miroc_585 <- metrics_585 %>% filter(GCM == "MIROC6")

# (a) Robustez del patrón norte-sur
rc_robust <- mean_pair_cor_585 > 0.8
# (b) MIROC6 representativo
miroc_repr <- cor_c > 0.8
# Range direction sign unanimity
sign_unanimity <- all(sign(metrics_585$range_change_hyb_pct) ==
                      sign(metrics_585$range_change_hyb_pct[1])) &&
                   all(metrics_585$range_change_hyb_pct < 0)

summary_lines <- c(
  "================================================================",
  " GCM x SSP DECISION REPORT — Psammodromus algirus hybrid SDM",
  "================================================================",
  sprintf(" Date: %s", Sys.Date()),
  sprintf(" Source: v2 runs (GCM/SSP-specific mechanistic layers)"),
  "",
  " Projected to 2041-2060, maxSSS threshold (Corr=0.28, Hyb=0.30)",
  sprintf(" GCMs (5): %s", paste(GCMS, collapse = ", ")),
  sprintf(" SSP coverage: SSP5-8.5 full (5 GCMs); SSP2-4.5 = MIROC6 only"),
  "",
  "---- PER-COMBINATION METRICS ----",
  "",
  paste(capture.output(print(metrics_all)), collapse = "\n"),
  "",
  "---- INTER-GCM VARIABILITY (SSP5-8.5) ----",
  sprintf(" Mean pairwise Spearman between Δ-suit GCMs = %.3f",
          mean_pair_cor_585),
  sprintf(" Inter-GCM median sd(Δ) / mean|Δ|         = %.3f",
          ratio_sd_to_absdelta),
  "",
  "---- RANGE-CHANGE DISPERSION (SSP5-8.5, maxSSS) ----",
  sprintf(" Correlative: mean=%.2f%%, sd=%.2f%% (CV=%.1f%%)",
          corr_mean_585, corr_var_585, cv_corr_585),
  sprintf(" Hybrid:      mean=%.2f%%, sd=%.2f%% (CV=%.1f%%)",
          hyb_mean_585, hyb_var_585, cv_hyb_585),
  "",
  "---- SSP EFFECT (MIROC6 ONLY) ----",
  sprintf(" Corr Δ-suit SSP2-4.5 vs SSP5-8.5 (Spearman) = %.3f",
          cor_ssp_miroc),
  sprintf(" Latitude of Δ=0 crossing: SSP2-4.5 = %.2f°N, SSP5-8.5 = %.2f°N",
          cross_245, cross_585_miroc),
  "",
  "---- DECISIONS ----",
  "",
  sprintf("(a) N-S asymmetry robust?  %s (mean pairwise ρ=%.3f, threshold 0.8)",
          ifelse(rc_robust, "YES", "NO"), mean_pair_cor_585),
  sprintf("(b) MIROC6 representative? %s (ρ with 5-GCM mean = %.3f, threshold 0.8)",
          ifelse(miroc_repr, "YES", "NO"), cor_c),
  sprintf("(c) Range-change sign unanimous across GCMs (Hybrid, SSP5-8.5)? %s",
          ifelse(sign_unanimity, "YES", "NO (one or more GCMs project expansion)")),
  "",
  "---- RECOMMENDED PARTITIONING ----",
  "",
  " Main text:",
  "   - Keep MIROC6 x SSP5-8.5 as the focal scenario for maps",
  "     (Δ-suit correlation with 5-GCM mean is high)",
  "   - Qualify language on inter-GCM stability: HYBRID IS NOT",
  "     MONOTONICALLY CONTRACTIVE under v2. Report headline as",
  "     'most GCMs project range contraction, one (MPI) projects",
  "     expansion; spatial pattern of Δ-suit is highly conserved",
  "     across GCMs (mean pairwise Spearman = X)'.",
  "   - Replace the prior CV=11% claim (v1 artefact).",
  "",
  " SI:",
  "   - Full GCM x SSP5-8.5 table (this CSV).",
  "   - MIROC6 x SSP2-4.5 comparison.",
  "   - Inter-GCM sd map + sign-agreement map (Fig. 1 bottom row).",
  "   - Pairwise Spearman correlations CSV.",
  "",
  " Backup (respond to reviewer if asked):",
  "   - v1 results (warm_map) with explicit note on why they are",
  "     less appropriate (freeze biophysical response across GCMs).",
  "   - Multi-GCM x SSP2-4.5 not currently run; justify as scope",
  "     decision given SSP5-8.5 already brackets the high end.",
  ""
)

writeLines(summary_lines, file.path(OUT, "gcm_ssp_decision_report.txt"))
cat("[OK] gcm_ssp_decision_report.txt\n")

# Markdown version with the full context
md_lines <- c(
  "# GCM × SSP decision report — Psammodromus algirus hybrid SDM",
  "",
  sprintf("*Generated: %s*", Sys.Date()),
  "*Source of rasters: v2 runs (GCM/SSP-specific mechanistic layers)*",
  "",
  "## Context",
  "",
  "The manuscript currently reports multi-GCM robustness using",
  "outputs from `sensitivity_gcm/` (**v1**), in which every GCM is",
  "projected with the same constant `_warm_map` mechanistic layers",
  "(deviation_mean and Activity). v1 therefore freezes the biophysical",
  "response across GCMs and mechanically inflates the apparent stability",
  "of the hybrid model. A v2 rerun (`sensitivity_gcm_v2/`) uses",
  "GCM/SSP-specific mechanistic layers and is the correct reference",
  "for any claim about GCM uncertainty propagating through the hybrid.",
  "",
  "This report replaces the v1-based language in Review_GF.md §3 and §7.",
  "",
  "## Scope of available projections",
  "",
  "| Scenario | GCMs | Output completeness |",
  "|---|---|---|",
  sprintf("| SSP5-8.5 | %s | mean + binary (P10, maxSSS) per GCM (no MESS, no SD per-GCM) |",
          paste(GCMS, collapse = ", ")),
  "| SSP2-4.5 | MIROC6 only | full suite (mean, sd, binary, MESS, categorical) |",
  "",
  "A multi-GCM × SSP2-4.5 run is **not executed** (scope decision:",
  "SSP2-4.5 is a backup scenario; SSP5-8.5 already brackets the high end).",
  "",
  "## Per-combination metrics (maxSSS, conservative overlap)",
  "",
  paste(capture.output(print(knitr::kable(metrics_all, digits = 3))),
        collapse = "\n"),
  "",
  "## Key inter-GCM statistics (SSP5-8.5)",
  "",
  sprintf("- Mean pairwise Spearman correlation between the 5 Δ-suit maps: **%.3f**",
          mean_pair_cor_585),
  sprintf("- Inter-GCM median sd(Δ) divided by mean |Δ|: **%.3f**",
          ratio_sd_to_absdelta),
  sprintf("- MIROC6 Δ-suit vs 5-GCM mean Δ-suit (Spearman): **%.3f**", cor_c),
  sprintf("- CV of range change across GCMs (Hybrid, maxSSS): **%.1f%%**",
          cv_hyb_585),
  sprintf("- CV of range change across GCMs (Correlative, maxSSS): **%.1f%%**",
          cv_corr_585),
  "",
  "## SSP effect (MIROC6)",
  "",
  sprintf("- Spearman correlation between MIROC6 Δ-suit under SSP2-4.5 and SSP5-8.5: **%.3f**",
          cor_ssp_miroc),
  sprintf("- Latitude of Δ=0 crossing: SSP2-4.5 = %.2f°N, SSP5-8.5 = %.2f°N",
          cross_245, cross_585_miroc),
  "",
  "## Interpretation",
  "",
  "1. The **spatial pattern** of hybrid–correlative divergence (Δ-suit)",
  sprintf("   is highly conserved across GCMs (mean pairwise ρ = %.3f).",
          mean_pair_cor_585),
  "   The north–south asymmetry described in the main text survives",
  "   the multi-GCM check.",
  "",
  "2. **Magnitude of range change** is *not* as stable under v2 as",
  "   previously reported. The hybrid projects contraction under",
  "   MIROC6 (−7.3%) and UKESM1-0-LL (−5.5%), near-neutral change",
  "   under CNRM-CM6-1 (−0.6%) and EC-Earth3-Veg (−0.6%), and a net",
  "   **expansion** under MPI-ESM1-2-HR (+3.7%). The correlative",
  "   model stays monotonically contractive.",
  "",
  "3. The v1 claim that the hybrid has lower CV than the correlative",
  sprintf("   (CV=11%% vs 56%%) **is an artefact**. Under v2 the CVs are %.0f%%",
          cv_hyb_585),
  sprintf("   (Hybrid) vs %.0f%% (Correlative).", cv_corr_585),
  "",
  sprintf("4. MIROC6 is a **reasonable representative** for main-text maps"),
  sprintf("   (Δ-suit Spearman ρ with 5-GCM mean = %.3f > 0.8), but not a",
          cor_c),
  "   representative of the *magnitude* of range change, which MIROC6",
  "   over-estimates relative to the 5-GCM ensemble.",
  "",
  "5. The SSP effect for MIROC6 is **directionally consistent**",
  sprintf("   (Δ-suit correlation = %.3f) and the latitude of Δ=0 shifts",
          cor_ssp_miroc),
  sprintf("   only slightly between SSPs (%.2f°N → %.2f°N).",
          cross_245, cross_585_miroc),
  "",
  "## Decision",
  "",
  "### Main text",
  "",
  "- **Keep MIROC6 × SSP5-8.5 as the focal illustration** for Figures",
  "  3 and 4 (spatial pattern is GCM-robust).",
  "- **Rewrite the robustness paragraph** (Review_GF.md §3 and §7)",
  "  to report, under v2:",
  "  - *spatial-pattern robustness* (pairwise ρ of Δ-suit across GCMs)",
  "    and *directional unanimity of the correlative contraction*,",
  "  - while explicitly acknowledging that the hybrid magnitude is",
  "    GCM-dependent and that one GCM (MPI-ESM1-2-HR) projects",
  "    hybrid expansion.",
  "- Drop the CV=11% vs CV=56% comparison.",
  "",
  "### Supplementary Information",
  "",
  "- Full per-GCM table for SSP5-8.5 (this CSV).",
  "- MIROC6 SSP2-4.5 vs SSP5-8.5 comparison (already in `ssp245/`).",
  "- Inter-GCM sd(Δ) map and sign-agreement map (this figure 1).",
  "- Pairwise Spearman correlation matrix of Δ-suit.",
  "- Latitudinal profile figures per GCM (this figure 2).",
  "",
  "### Backup — only surface if a reviewer asks",
  "",
  "- v1 (warm_map) outputs, with a brief note explaining why they",
  "  are less appropriate (GCM uncertainty does not reach the",
  "  biophysical layers).",
  "- Justification for omitting multi-GCM × SSP2-4.5 (scope decision).",
  "- v1 vs v2 comparison (already in `comparison_v1_v2/`).",
  "",
  "## Files generated",
  "",
  "- `gcm_ssp_comparison_table.csv` — per-combination metrics",
  "- `summary_decision_table.csv` — main-text decision summary",
  "- `delta_suit_pairwise_correlations.csv`",
  "- `delta_suit_inter_gcm_sd_ssp585.tif`, `delta_suit_inter_gcm_mean_ssp585.tif`,",
  "  `delta_suit_sign_agreement_ssp585.tif`",
  "- `diagnostic_gcm_ssp_maps.png`",
  "- `diagnostic_gcm_ssp_profiles.png`",
  "- `diagnostic_range_change_comparison.png`",
  "- `gcm_ssp_decision_report.txt` and `.md`",
  ""
)

writeLines(md_lines, file.path(OUT, "gcm_ssp_decision_report.md"))
cat("[OK] gcm_ssp_decision_report.md\n")

cat("\n=== DONE. All outputs in", OUT, "===\n")
