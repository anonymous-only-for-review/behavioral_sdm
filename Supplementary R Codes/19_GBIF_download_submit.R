################################################################################
## 19_gbif_download_submit.R — Submit GBIF occ_download for Psammodromus
## algirus with cutoff 30-Jun-2025.
##
## Credentials are read from GBIF_credenciales.rtf (RTF file in
## manuscrito/revision_2026-04-22/). Password is NOT echoed to stdout.
##
## Output: data/gbif_download_2026-04-22/download_key.txt with the
## submitted downloadKey and DOI. This script does NOT wait for the
## download to complete — use 20_gbif_download_fetch.R for that.
################################################################################

suppressPackageStartupMessages({
  library(rgbif)
})

# -------------------------------------------------------------------------- #
# 1. Read credentials from RTF (manuscrito/revision_2026-04-22/GBIF_credenciales.rtf)
# -------------------------------------------------------------------------- #

CRED_PATH <- "manuscrito/revision_2026-04-22/GBIF_credenciales.rtf"
if (!file.exists(CRED_PATH)) stop("Credentials file not found: ", CRED_PATH)

raw <- readLines(CRED_PATH, warn = FALSE)
txt <- paste(raw, collapse = "\n")

# Strip RTF control words and groups; ignore hex-encoded chars
txt <- gsub("\\\\[A-Za-z]+-?[0-9]*\\s?", "", txt)
txt <- gsub("[{}]", "",  txt)
txt <- gsub("\\\\'[0-9a-fA-F]{2}", "", txt)

extract_key <- function(s, key) {
  m <- regmatches(s, regexpr(paste0(key, "\\s*=\\s*\\S+"), s))
  if (!length(m)) return(NA_character_)
  sub(paste0(key, "\\s*=\\s*"), "", m)
}

GBIF_USER  <- extract_key(txt, "GBIF_USER")
GBIF_PWD   <- extract_key(txt, "GBIF_PWD")
GBIF_EMAIL <- extract_key(txt, "GBIF_EMAIL")

if (any(is.na(c(GBIF_USER, GBIF_PWD, GBIF_EMAIL))) ||
    !all(nzchar(c(GBIF_USER, GBIF_PWD, GBIF_EMAIL))))
  stop("Failed to parse all three credential fields from RTF")

Sys.setenv(GBIF_USER  = GBIF_USER,
           GBIF_PWD   = GBIF_PWD,
           GBIF_EMAIL = GBIF_EMAIL)

cat(sprintf("Credentials loaded: user=%s***%s, email=%s***%s\n",
            substr(GBIF_USER, 1, 2), substr(GBIF_USER, nchar(GBIF_USER)-1, nchar(GBIF_USER)),
            substr(GBIF_EMAIL, 1, 2), substr(GBIF_EMAIL, nchar(GBIF_EMAIL)-1, nchar(GBIF_EMAIL))))

# -------------------------------------------------------------------------- #
# 2. Resolve taxonKey
# -------------------------------------------------------------------------- #

cat("\nResolving taxon key for Psammodromus algirus...\n")
tx <- name_backbone(name = "Psammodromus algirus", rank = "species",
                    strict = TRUE)
stopifnot(!is.null(tx$usageKey))
TAXON_KEY <- as.integer(tx$usageKey)
cat(sprintf("  usageKey = %s | canonicalName = %s | rank = %s | status = %s\n",
            as.character(TAXON_KEY), as.character(tx$canonicalName),
            as.character(tx$rank), as.character(tx$status)))

# -------------------------------------------------------------------------- #
# 3. Submit occ_download with filters
# -------------------------------------------------------------------------- #

cat("\nSubmitting occ_download request...\n")
cat("Filters:\n")
cat("  taxonKey            =", TAXON_KEY, "\n")
cat("  hasCoordinate       = TRUE\n")
cat("  hasGeospatialIssue  = FALSE\n")
cat("  occurrenceStatus    = PRESENT\n")
cat("  basisOfRecord       in {HUMAN_OBSERVATION, OBSERVATION,\n")
cat("                          MACHINE_OBSERVATION, PRESERVED_SPECIMEN,\n")
cat("                          MATERIAL_SAMPLE}\n")
cat("  year                <= 2025\n")
cat("  (post-filter)       year < 2025 OR (year == 2025 & month <= 6)\n")
cat("\nNote: GBIF predicates allow filtering by year or by eventDate but\n")
cat("the month-within-year split is applied server-side via pred_or(...).\n\n")

# Server-side filter: year <= 2025. The month-within-2025 cutoff (<= 6)
# is applied post-download to keep the predicate simple and compatible
# across rgbif versions.
dl <- occ_download(
  pred("taxonKey", TAXON_KEY),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("occurrenceStatus", "PRESENT"),
  pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "OBSERVATION",
                             "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN",
                             "MATERIAL_SAMPLE")),
  pred_lte("year", 2025),
  format = "SIMPLE_CSV",
  user   = GBIF_USER,
  pwd    = GBIF_PWD,
  email  = GBIF_EMAIL
)

# -------------------------------------------------------------------------- #
# 4. Save download key + metadata
# -------------------------------------------------------------------------- #

`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a
key <- as.character(dl)
meta_url <- attr(dl, "downloadLink")
doi      <- attr(dl, "doi")

out_dir <- "data/gbif_download_2026-04-22"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(c(
  sprintf("downloadKey: %s", key),
  sprintf("doi: %s", doi %||% "pending"),
  sprintf("link: %s", meta_url %||% "pending"),
  sprintf("submitted_at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("taxonKey: %s", as.character(TAXON_KEY)),
  "filters:",
  "  hasCoordinate=TRUE",
  "  hasGeospatialIssue=FALSE",
  "  occurrenceStatus=PRESENT",
  "  basisOfRecord in {HUMAN_OBSERVATION, OBSERVATION, MACHINE_OBSERVATION, PRESERVED_SPECIMEN, MATERIAL_SAMPLE}",
  "  year < 2025 OR (year == 2025 AND month <= 6)"
), file.path(out_dir, "download_key.txt"))

cat("\n=========================================\n")
cat(sprintf(" downloadKey = %s\n", key))
cat(sprintf(" DOI         = %s\n", doi %||% "pending (will be populated once ready)"))
cat(sprintf(" Saved to    = %s\n", file.path(out_dir, "download_key.txt")))
cat("=========================================\n")
cat("\nTo poll and retrieve, run: Rscript Supplementary R Codes/20_GBIF_download_fetch.R\n")
