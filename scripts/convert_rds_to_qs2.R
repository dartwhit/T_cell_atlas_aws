#!/usr/bin/env Rscript

# Convert Seurat/gene-list .rds files to .qs2 for faster app loads.
#
# Why: readRDS on gzip-compressed Seurat objects is single-threaded and sits on
# the critical path of "Load Data" (5-14s in production). qs2 (zstd) loads the
# same object ~7x faster while keeping files small. The app's read_object()
# helper (app_code/setup.R) automatically prefers a .qs2 sibling when present
# and falls back to .rds, so this conversion is safe to run incrementally and
# the originals can be kept as a fallback.
#
# Usage:
#   Rscript scripts/convert_rds_to_qs2.R [DATA_DIR]
#
#   DATA_DIR  Root to scan recursively for *.rds (default: app_code/data).
#
# Env:
#   QS2_NTHREADS   Threads for qs_save (default: all available cores).
#   OVERWRITE      "1" to re-convert even if an up-to-date .qs2 exists.
#
# Notes:
#   - DE/VAM tables are plain text (.txt/.csv) and are NOT touched.
#   - The .qs2 files must be written by the same major qs2 version the app runs
#     with. Convert with the same qs2 used in the Docker image to be safe.

suppressPackageStartupMessages({
  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("Package 'qs2' is not installed. Run install.packages('qs2') first.")
  }
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "app_code/data"
overwrite <- identical(Sys.getenv("OVERWRITE"), "1")
nthreads <- {
  n <- suppressWarnings(as.integer(Sys.getenv("QS2_NTHREADS")))
  if (is.na(n) || n < 1) max(1L, parallel::detectCores()) else n
}

if (!dir.exists(data_dir)) stop("Data dir not found: ", data_dir)

rds_files <- list.files(data_dir, pattern = "\\.rds$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
if (length(rds_files) == 0) {
  cat("No .rds files found under", data_dir, "\n")
  quit(status = 0)
}

cat(sprintf("Found %d .rds file(s) under %s  (qs2 threads: %d)\n",
            length(rds_files), data_dir, nthreads))

total_rds <- 0
total_qs2 <- 0
converted <- 0
skipped <- 0
failed <- 0

for (rds in rds_files) {
  qs2 <- sub("\\.rds$", ".qs2", rds, ignore.case = TRUE)

  # Skip if an up-to-date .qs2 already exists (newer than the source).
  if (!overwrite && file.exists(qs2) &&
      file.info(qs2)$mtime >= file.info(rds)$mtime) {
    cat(sprintf("  skip   %s (.qs2 up to date)\n", basename(rds)))
    skipped <- skipped + 1
    next
  }

  res <- tryCatch({
    t_read <- system.time(obj <- readRDS(rds))["elapsed"]
    qs2::qs_save(obj, qs2, nthreads = nthreads)
    # Round-trip sanity check: confirm the file reads back.
    t_qs <- system.time(invisible(qs2::qs_read(qs2, nthreads = nthreads)))["elapsed"]
    list(t_read = t_read, t_qs = t_qs,
         sz_rds = file.info(rds)$size, sz_qs = file.info(qs2)$size)
  }, error = function(e) {
    cat(sprintf("  FAIL   %s : %s\n", basename(rds), conditionMessage(e)))
    if (file.exists(qs2)) unlink(qs2)  # don't leave a half-written file
    NULL
  })

  if (is.null(res)) { failed <- failed + 1; next }

  total_rds <- total_rds + res$sz_rds
  total_qs2 <- total_qs2 + res$sz_qs
  converted <- converted + 1
  cat(sprintf("  ok     %-40s  load %.2fs -> %.2fs  (%.0f -> %.0f MB)\n",
              basename(rds), res$t_read, res$t_qs,
              res$sz_rds / 1e6, res$sz_qs / 1e6))
}

cat(sprintf("\nDone. converted=%d skipped=%d failed=%d\n",
            converted, skipped, failed))
if (converted > 0) {
  cat(sprintf("Total size: %.0f MB (.rds) -> %.0f MB (.qs2)\n",
              total_rds / 1e6, total_qs2 / 1e6))
}
if (failed > 0) quit(status = 1)
