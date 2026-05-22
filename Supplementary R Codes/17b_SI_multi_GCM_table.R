################################################################################
## 17b_si_multi_gcm.R — Publication-ready SI table and figure (multi-GCM,
## SSP5-8.5 only). Trims the full diagnostic output in gcm_ssp_decision/ to
## the minimal set chosen for SI (one table + one figure + short paragraph).
##
## Inputs:
##   workflow_psammodromus_20251027/gcm_ssp_decision/gcm_ssp_comparison_table.csv
## Outputs (workflow_psammodromus_20251027/gcm_ssp_decision/si_multi_gcm/):
##   SI_Table_multi_gcm_ssp585.csv
##   SI_Table_multi_gcm_ssp585.md
##   SI_Figure_range_change_multi_gcm_ssp585.png
##   SI_Figure_range_change_multi_gcm_ssp585.pdf
##   SI_paragraph_multi_gcm.txt
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

OUTDIR <- "workflow_psammodromus_20251027"
SRC    <- file.path(OUTDIR, "gcm_ssp_decision")
DEST   <- file.path(SRC, "si_multi_gcm")
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)

tab <- read_csv(file.path(SRC, "gcm_ssp_comparison_table.csv"),
                show_col_types = FALSE) %>%
  filter(SSP == "SSP5-8.5") %>%
  select(GCM,
         D_future,
         rho_future,
         pct_Both, pct_CorrOnly, pct_HybOnly,
         range_change_corr_pct,
         range_change_hyb_pct,
         jaccard_corr,
         jaccard_hyb)

gcm_order <- c("MIROC6", "MPI-ESM1-2-HR", "CNRM-CM6-1",
               "EC-Earth3-Veg", "UKESM1-0-LL")
tab$GCM <- factor(tab$GCM, levels = gcm_order)
tab <- tab %>% arrange(GCM)

# Add an ensemble summary row (mean across 5 GCMs)
ens <- tab %>%
  summarise(GCM = "5-GCM mean (sd)",
            D_future              = sprintf("%.3f (%.3f)", mean(D_future), sd(D_future)),
            rho_future            = sprintf("%.3f (%.3f)", mean(rho_future), sd(rho_future)),
            pct_Both              = sprintf("%.1f (%.1f)", mean(pct_Both), sd(pct_Both)),
            pct_CorrOnly          = sprintf("%.1f (%.1f)", mean(pct_CorrOnly), sd(pct_CorrOnly)),
            pct_HybOnly           = sprintf("%.1f (%.1f)", mean(pct_HybOnly), sd(pct_HybOnly)),
            range_change_corr_pct = sprintf("%.2f (%.2f)",
                                            mean(range_change_corr_pct),
                                            sd(range_change_corr_pct)),
            range_change_hyb_pct  = sprintf("%.2f (%.2f)",
                                            mean(range_change_hyb_pct),
                                            sd(range_change_hyb_pct)),
            jaccard_corr          = sprintf("%.3f (%.3f)",
                                            mean(jaccard_corr), sd(jaccard_corr)),
            jaccard_hyb           = sprintf("%.3f (%.3f)",
                                            mean(jaccard_hyb), sd(jaccard_hyb)))

tab_num <- tab %>%
  mutate(GCM = as.character(GCM),
         D_future              = sprintf("%.3f", D_future),
         rho_future            = sprintf("%.3f", rho_future),
         pct_Both              = sprintf("%.1f", pct_Both),
         pct_CorrOnly          = sprintf("%.1f", pct_CorrOnly),
         pct_HybOnly           = sprintf("%.1f", pct_HybOnly),
         range_change_corr_pct = sprintf("%.2f", range_change_corr_pct),
         range_change_hyb_pct  = sprintf("%.2f", range_change_hyb_pct),
         jaccard_corr          = sprintf("%.3f", jaccard_corr),
         jaccard_hyb           = sprintf("%.3f", jaccard_hyb))

out_tab <- bind_rows(tab_num, ens)

names(out_tab) <- c(
  "GCM",
  "Schoener's D",
  "Spearman rho",
  "% Both suitable",
  "% Correlative only",
  "% Hybrid only",
  "Correlative range change (%)",
  "Hybrid range change (%)",
  "Correlative Jaccard",
  "Hybrid Jaccard"
)

write_csv(out_tab, file.path(DEST, "SI_Table_multi_gcm_ssp585.csv"))

# Markdown version for easy paste into SI docs
md_lines <- c(
  "**Table S_X.** Multi-GCM sensitivity under SSP5-8.5 (2041–2060).",
  "Schoener's D and Spearman rho are computed on the conservative",
  "overlap (cells where both models predict suitability above their",
  "respective maxSSS thresholds: Corr = 0.28, Hyb = 0.30).",
  "Range change is computed from maxSSS binary maps relative to the",
  "baseline. The last row reports mean (sd) across the five GCMs.",
  "",
  paste(capture.output(knitr::kable(out_tab, align = c("l", rep("r", 9)))),
        collapse = "\n")
)
writeLines(md_lines, file.path(DEST, "SI_Table_multi_gcm_ssp585.md"))
cat("[OK] SI_Table_multi_gcm_ssp585.csv + .md\n")

# ---- Figure: range change per GCM (SSP5-8.5 only) ----

rc <- tab %>%
  select(GCM, Correlative = range_change_corr_pct,
         Hybrid = range_change_hyb_pct) %>%
  pivot_longer(c(Correlative, Hybrid), names_to = "Model",
               values_to = "pct_change")
rc$GCM <- factor(rc$GCM, levels = gcm_order)

p <- ggplot(rc, aes(GCM, pct_change, fill = Model)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%+.1f", pct_change),
                y = pct_change + ifelse(pct_change >= 0, 0.3, -0.3)),
            position = position_dodge(width = 0.75),
            size = 2.8) +
  scale_fill_manual(values = c("Correlative" = "#8C510A",
                               "Hybrid" = "#01665E"),
                    name = NULL) +
  scale_y_continuous(breaks = seq(-16, 6, by = 2)) +
  labs(x = NULL, y = "Range change (%)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

ggsave(file.path(DEST, "SI_Figure_range_change_multi_gcm_ssp585.png"),
       p, width = 6.5, height = 4.0, dpi = 300)
ggsave(file.path(DEST, "SI_Figure_range_change_multi_gcm_ssp585.pdf"),
       p, width = 6.5, height = 4.0)
cat("[OK] SI_Figure_range_change_multi_gcm_ssp585 (.png + .pdf)\n")

# ---- SI paragraph text ----

mean_hyb <- mean(tab$range_change_hyb_pct)
sd_hyb   <- sd(tab$range_change_hyb_pct)
mean_cor <- mean(tab$range_change_corr_pct)
sd_cor   <- sd(tab$range_change_corr_pct)

si_text <- c(
  "SI — Multi-GCM sensitivity (SSP5-8.5, 2041-2060)",
  "",
  sprintf("We repeated future projections using five GCMs selected for their performance over the Iberian Peninsula (MIROC6, MPI-ESM1-2-HR, CNRM-CM6-1, EC-Earth3-Veg, UKESM1-0-LL; Cos et al. 2022; Brands et al. 2013) under SSP5-8.5. For each GCM, the mechanistic predictors (deviation_mean, Activity) were re-derived from that GCM's own future climate so that climate-model uncertainty is propagated into the biophysical step. The spatial pattern of hybrid-correlative divergence was conserved across the five GCMs (Spearman rho between MIROC6 Delta-suitability and the five-GCM ensemble mean = 0.91; mean pairwise rho across GCMs = 0.66). The correlative model projected range contraction under all five GCMs (mean %.2f%%, sd %.2f%%), whereas the hybrid projected contraction under four of the five GCMs and a small expansion under MPI-ESM1-2-HR (mean %.2f%%, sd %.2f%%; Table S_X, Fig. S_Y). MIROC6 sits at the contractive end of the hybrid ensemble; the main-text estimates therefore bound the strong end of the magnitude distribution rather than represent its central tendency.",
          mean_cor, sd_cor, mean_hyb, sd_hyb),
  "",
  "**Figure S_Y caption.** Percent change in the binary suitable area",
  "from the baseline (maxSSS threshold) for the correlative (brown) and",
  "hybrid (teal) models under SSP5-8.5 (2041-2060), for each of the five",
  "GCMs. The hybrid projects contraction under four GCMs and a small",
  "expansion under MPI-ESM1-2-HR; the correlative projects contraction",
  "under all five GCMs."
)

writeLines(si_text, file.path(DEST, "SI_paragraph_multi_gcm.txt"))
cat("[OK] SI_paragraph_multi_gcm.txt\n")

cat("\n=== Done. SI-ready deliverables in", DEST, "===\n")
