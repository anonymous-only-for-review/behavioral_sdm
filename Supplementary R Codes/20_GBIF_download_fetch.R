################################################################################
## 20_gbif_download_fetch.R — Poll the GBIF download, retrieve the zip,
## apply the local cutoff (month <= 6 for year 2025), and compare to the
## 29-Jan-2025 snapshot in data/gbif_lizard_cleaned.RData.
##
## Reads downloadKey from data/gbif_download_2026-04-22/download_key.txt.
################################################################################

suppressPackageStartupMessages({
  library(rgbif)
  library(dplyr)
  library(readr)
  library(sf)
})
select <- dplyr::select

OUT <- "data/gbif_download_2026-04-22"
stopifnot(file.exists(file.path(OUT, "download_key.txt")))
info <- readLines(file.path(OUT, "download_key.txt"))
key  <- sub("downloadKey: ", "", grep("^downloadKey:", info, value = TRUE))
cat("Polling downloadKey:", key, "\n")

# -------------------------------------------------------------------------- #
# 1. Wait until ready
# -------------------------------------------------------------------------- #

meta <- occ_download_wait(key, status_ping = 30, curlopts = list(timeout_ms = 600000))
cat("Status:", meta$status, "\n")
stopifnot(meta$status == "SUCCEEDED")

doi      <- meta$doi
zip_url  <- meta$downloadLink
n_records <- meta$totalRecords
cat(sprintf("DOI: %s\nTotal records: %d\n", doi, n_records))

# Persist the final metadata
writeLines(c(
  sprintf("downloadKey: %s",  key),
  sprintf("doi: %s",          doi),
  sprintf("link: %s",         zip_url),
  sprintf("status: %s",       meta$status),
  sprintf("totalRecords: %d", n_records),
  sprintf("created: %s",      meta$created),
  sprintf("modified: %s",     meta$modified)
), file.path(OUT, "download_metadata.txt"))

# -------------------------------------------------------------------------- #
# 2. Download the zip and unpack
# -------------------------------------------------------------------------- #

zip_path <- occ_download_get(key, path = OUT, overwrite = TRUE)
cat("Zip saved to:", as.character(zip_path), "\n")

# Import as data frame
df <- occ_download_import(zip_path)
cat(sprintf("Imported %d rows, %d columns\n", nrow(df), ncol(df)))

# -------------------------------------------------------------------------- #
# 3. Apply local cutoff: drop records with year==2025 & month > 6
# -------------------------------------------------------------------------- #

n_pre <- nrow(df)
df_cut <- df %>%
  filter(!(year == 2025 & month > 6))
n_post <- nrow(df_cut)
n_dropped_jul2025 <- n_pre - n_post
cat(sprintf("\nCutoff (month <= 6 for year 2025): dropped %d records\n",
            n_dropped_jul2025))

# -------------------------------------------------------------------------- #
# 4. Compare with the 29-Jan-2025 snapshot
# -------------------------------------------------------------------------- #

cat("\n=== Comparison vs data/gbif_lizard_cleaned.RData (2025-01-29) ===\n")
load("data/gbif_lizard_cleaned.RData")
old <- gbif_lizard_cleaned

# Key comparisons
n_old <- nrow(old)
n_new <- nrow(df_cut)

cat(sprintf("N (old, 2025-01-29 cleaned): %d\n", n_old))
cat(sprintf("N (new, 2026-04-22 + cutoff 2025-06): %d\n", n_new))
cat(sprintf("Difference: %+d (%.1f%%)\n", n_new - n_old,
            100 * (n_new - n_old) / n_old))

# Year distribution
yr_old <- as.data.frame(table(old$year))
names(yr_old) <- c("year", "n_old")
yr_new <- as.data.frame(table(df_cut$year))
names(yr_new) <- c("year", "n_new")
yr_cmp <- full_join(yr_old, yr_new, by = "year") %>%
  mutate(year = as.integer(as.character(year)),
         n_old = tidyr::replace_na(n_old, 0L),
         n_new = tidyr::replace_na(n_new, 0L),
         delta = n_new - n_old) %>%
  arrange(year)

cat("\nYear-by-year delta (top 5 gains and losses):\n")
print(yr_cmp %>% arrange(desc(delta)) %>% head(5))
cat("\n")
print(yr_cmp %>% arrange(delta) %>% head(5))

# Spatial overlap: compare by gbifID / key
key_col_old <- if ("key" %in% names(old)) "key" else "gbifID"
ids_old <- as.character(old[[key_col_old]])

key_col_new <- if ("gbifID" %in% names(df_cut)) "gbifID" else
                 if ("key" %in% names(df_cut)) "key" else NA
if (!is.na(key_col_new)) {
  ids_new <- as.character(df_cut[[key_col_new]])
  common  <- intersect(ids_old, ids_new)
  only_old <- setdiff(ids_old, ids_new)
  only_new <- setdiff(ids_new, ids_old)
  cat(sprintf("\nRecord-ID overlap:\n"))
  cat(sprintf("  in both      : %d\n", length(common)))
  cat(sprintf("  only old     : %d (%.1f%%)\n", length(only_old),
              100 * length(only_old) / n_old))
  cat(sprintf("  only new     : %d (%.1f%%)\n", length(only_new),
              100 * length(only_new) / n_new))
}

# Lat/lon range
cat(sprintf("\nLat range old: [%.2f, %.2f]\n",
            min(old$decimalLatitude, na.rm = TRUE),
            max(old$decimalLatitude, na.rm = TRUE)))
cat(sprintf("Lat range new: [%.2f, %.2f]\n",
            min(df_cut$decimalLatitude, na.rm = TRUE),
            max(df_cut$decimalLatitude, na.rm = TRUE)))
cat(sprintf("Lon range old: [%.2f, %.2f]\n",
            min(old$decimalLongitude, na.rm = TRUE),
            max(old$decimalLongitude, na.rm = TRUE)))
cat(sprintf("Lon range new: [%.2f, %.2f]\n",
            min(df_cut$decimalLongitude, na.rm = TRUE),
            max(df_cut$decimalLongitude, na.rm = TRUE)))

# Save cleaned subset as RData + CSV for later use
save(df_cut, file = file.path(OUT, "gbif_2026-04-22_cut2025-06.RData"))
write_csv(df_cut %>% select(any_of(c("gbifID", "scientificName",
                                     "decimalLatitude", "decimalLongitude",
                                     "coordinateUncertaintyInMeters",
                                     "eventDate", "year", "month", "day",
                                     "basisOfRecord", "datasetKey"))),
          file.path(OUT, "gbif_2026-04-22_cut2025-06_coords.csv"))

# Comparison table
write_csv(yr_cmp, file.path(OUT, "year_comparison.csv"))

cat("\n[OK] All outputs in", OUT, "\n")
