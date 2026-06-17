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
#   Rscript scripts/convert_rds_to_qs2.R <DATA_DIR>
#
#   DATA_DIR  Root to scan recursively for *.rds. Required.
#             Local:  app_code/data
#             EC2:    the mounted data path, e.g. /srv/shiny-server/atlas/data
#
# Env:
#   QS2_NTHREADS   Threads for qs_save (default: all available cores).
#   OVERWRITE      "1" to re-convert even if an up-to-date .qs2 exists.
#   SLIM           "1" to drop app-unused assays/layers (`integrated`,
#                  `scale.data`) from non-spatial Seurat objects before saving.
#                  ~10x smaller files; requires the 'Seurat' package.
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
if (length(args) < 1 || !nzchar(args[[1]])) {
  cat("Usage: Rscript scripts/convert_rds_to_qs2.R <DATA_DIR>\n",
      "  <DATA_DIR>  Root to scan recursively for *.rds.\n",
      "              On EC2 pass the mounted data path, e.g.\n",
      "              /srv/shiny-server/atlas/data\n",
      "Env: QS2_NTHREADS (save threads), OVERWRITE=1 (re-convert).\n",
      sep = "")
  quit(status = 2)
}
data_dir <- args[[1]]
overwrite <- identical(Sys.getenv("OVERWRITE"), "1")
nthreads <- {
  cores <- parallel::detectCores()
  if (is.na(cores) || cores < 1) cores <- 1L
  n <- suppressWarnings(as.integer(Sys.getenv("QS2_NTHREADS")))
  if (is.na(n) || n < 1) cores else min(n, cores)
}

# Optional slimming: drop assays/layers the app never reads (the `integrated`
# assay and all `scale.data` layers are ~80% of the object's weight but unused
# by the Shiny app). Cuts file size ~10x, which is the dominant cost on EC2's
# mounted volume. Spatial objects are left untouched (the spatial viewer needs
# their full structure). Enable with SLIM=1.
slim_enabled <- identical(Sys.getenv("SLIM"), "1")

# Assays/reductions/layers the app actually uses (intersected per-object, so
# it's safe across studies/levels that have different subsets).
SEURAT_KEEP_ASSAYS <- c("RNA", "SCT", "VAMdist", "VAMcdf")
SEURAT_KEEP_REDUCS <- c("pca", "umap")
SEURAT_KEEP_LAYERS <- c("counts", "data")

if (slim_enabled && !requireNamespace("Seurat", quietly = TRUE)) {
  stop("SLIM=1 requires the 'Seurat' package.")
}

# Returns a slimmed copy for non-spatial Seurat objects; everything else
# (gene-list lists, spatial Seurat objects, unrecognized layouts) is returned
# unchanged. Falls back to the original object if DietSeurat errors.
slim_seurat <- function(obj) {
  if (!slim_enabled || !inherits(obj, "Seurat")) return(obj)
  if (length(obj@images) > 0 || "Spatial" %in% SeuratObject::Assays(obj)) {
    return(obj)  # spatial: keep whole
  }
  assays <- intersect(SEURAT_KEEP_ASSAYS, SeuratObject::Assays(obj))
  if (length(assays) == 0) return(obj)  # unfamiliar layout: don't risk it
  reducs <- intersect(SEURAT_KEEP_REDUCS, SeuratObject::Reductions(obj))
  tryCatch(
    Seurat::DietSeurat(obj, layers = SEURAT_KEEP_LAYERS,
                       assays = assays, dimreducs = reducs),
    error = function(e) {
      cat(sprintf("    (slim skipped: %s)\n", conditionMessage(e)))
      obj
    }
  )
}

if (!dir.exists(data_dir)) stop("Data dir not found: ", data_dir)

rds_files <- list.files(data_dir, pattern = "\\.rds$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
if (length(rds_files) == 0) {
  cat("No .rds files found under", data_dir, "\n")
  quit(status = 0)
}

cat(sprintf("Found %d .rds file(s) under %s  (qs2 threads: %d, slim: %s)\n",
            length(rds_files), data_dir, nthreads, if (slim_enabled) "on" else "off"))

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
    obj <- slim_seurat(obj)
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
  cat(sprintf("Converted files only: %.0f MB (.rds) -> %.0f MB (.qs2)\n",
              total_rds / 1e6, total_qs2 / 1e6))
}
if (failed > 0) quit(status = 1)
