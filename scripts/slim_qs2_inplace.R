#!/usr/bin/env Rscript

# Slim existing .qs2 Seurat objects IN PLACE (for when the source .rds files
# are already gone). Drops assays/layers the Shiny app never reads -- the
# `integrated` assay and all `scale.data` layers -- which are ~80% of a typical
# object's weight. This shrinks files ~10x, the dominant load cost on EC2's
# mounted volume.
#
# Safe to re-run: objects already slim (no `integrated` assay, no `scale.data`)
# are skipped. Spatial objects, gene-list lists, and unfamiliar layouts are left
# untouched.
#
# Usage:
#   Rscript scripts/slim_qs2_inplace.R <DATA_DIR>
#
#   DATA_DIR  Root to scan recursively for *.qs2. Required.
#             EC2: the mounted data path, e.g. /srv/shiny-server/atlas/data
#
# Env:
#   QS2_NTHREADS   Threads for qs_read/qs_save (default: detected cores).
#   BACKUP         "1" to keep the original as <file>.qs2.bak before overwrite.
#
# NOTE: writes are atomic (temp file + rename) so an interrupted run can't leave
# a half-written .qs2 in place of a good one.

suppressPackageStartupMessages({
  if (!requireNamespace("qs2", quietly = TRUE))
    stop("Package 'qs2' is not installed.")
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Package 'Seurat' is not installed.")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[[1]])) {
  cat("Usage: Rscript scripts/slim_qs2_inplace.R <DATA_DIR>\n",
      "  <DATA_DIR>  Root to scan recursively for *.qs2 (e.g.\n",
      "              /srv/shiny-server/atlas/data on EC2).\n",
      "Env: QS2_NTHREADS, BACKUP=1 (keep .qs2.bak).\n", sep = "")
  quit(status = 2)
}
data_dir <- args[[1]]
backup <- identical(Sys.getenv("BACKUP"), "1")
nthreads <- {
  cores <- parallel::detectCores()
  if (is.na(cores) || cores < 1) cores <- 1L
  n <- suppressWarnings(as.integer(Sys.getenv("QS2_NTHREADS")))
  if (is.na(n) || n < 1) cores else min(n, cores)
}
if (!dir.exists(data_dir)) stop("Data dir not found: ", data_dir)

# Assays/reductions/layers the app actually uses (intersected per-object).
KEEP_ASSAYS <- c("RNA", "SCT", "VAMdist", "VAMcdf")
KEEP_REDUCS <- c("pca", "umap")
KEEP_LAYERS <- c("counts", "data")

# Does this object still carry app-unused weight worth dropping?
needs_slim <- function(obj) {
  if ("integrated" %in% SeuratObject::Assays(obj)) return(TRUE)
  for (a in SeuratObject::Assays(obj)) {
    if ("scale.data" %in% SeuratObject::Layers(obj[[a]])) return(TRUE)
  }
  FALSE
}

qs2_files <- list.files(data_dir, pattern = "\\.qs2$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
if (length(qs2_files) == 0) {
  cat("No .qs2 files found under", data_dir, "\n"); quit(status = 0)
}
cat(sprintf("Found %d .qs2 file(s) under %s  (threads: %d)\n",
            length(qs2_files), data_dir, nthreads))

# Process one file. Returns a one-line status string; on a successful slim it
# also writes the new file and updates the size accumulators in the caller.
# Defined as a function so early-exit `return()`s for the skip cases work.
process_file <- function(f) {
  obj <- qs2::qs_read(f, nthreads = nthreads)

  if (!inherits(obj, "Seurat")) return("skip (not a Seurat object)")
  if (length(obj@images) > 0 || "Spatial" %in% SeuratObject::Assays(obj))
    return("skip (spatial)")
  if (!needs_slim(obj)) return("skip (already slim)")
  assays <- intersect(KEEP_ASSAYS, SeuratObject::Assays(obj))
  if (length(assays) == 0) return("skip (unfamiliar layout)")

  reducs <- intersect(KEEP_REDUCS, SeuratObject::Reductions(obj))
  slim <- Seurat::DietSeurat(obj, layers = KEEP_LAYERS,
                             assays = assays, dimreducs = reducs)

  before <- file.info(f)$size
  if (backup) file.copy(f, paste0(f, ".bak"), overwrite = TRUE)
  tmp <- paste0(f, ".tmp")
  qs2::qs_save(slim, tmp, nthreads = nthreads)
  invisible(qs2::qs_read(tmp, nthreads = nthreads))  # round-trip check
  file.rename(tmp, f)                                # atomic replace
  after <- file.info(f)$size

  sz_before <<- sz_before + before
  sz_after  <<- sz_after  + after
  sprintf("slimmed  %.0f -> %.0f MB", before / 1e6, after / 1e6)
}

slimmed <- 0; skipped <- 0; failed <- 0; sz_before <- 0; sz_after <- 0

for (f in qs2_files) {
  res <- tryCatch(process_file(f), error = function(e) {
    tmp <- paste0(f, ".tmp"); if (file.exists(tmp)) unlink(tmp)
    paste("FAIL:", conditionMessage(e))
  })

  if (grepl("^slimmed", res))      slimmed <- slimmed + 1
  else if (grepl("^FAIL", res))    failed  <- failed  + 1
  else                             skipped <- skipped + 1
  cat(sprintf("  %-50s %s\n", basename(f), res))
}

cat(sprintf("\nDone. slimmed=%d skipped=%d failed=%d\n",
            slimmed, skipped, failed))
if (slimmed > 0)
  cat(sprintf("Slimmed files: %.0f MB -> %.0f MB\n",
              sz_before / 1e6, sz_after / 1e6))
if (failed > 0) quit(status = 1)
