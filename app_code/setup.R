inDir <- "data/"
# DE/VAM tables now live in data/<study>/de/ ; everything else in data/<study>/

# Resolve a data file to its per-dataset folder. `de = TRUE` for DE/VAM tables.
data_path <- function(study, name, de = FALSE) {
  rel <- if (de) file.path(study, "de", name) else file.path(study, name)
  paste0(inDir, rel)
}

# Load a serialized R object, preferring a fast qs2 sibling when present.
# Given a path like ".../foo.rds", if ".../foo.qs2" exists it is read with
# qs2::qs_read() (zstd: ~7x faster to load than gzipped RDS, small on disk).
# Otherwise falls back to readRDS(), so the app works before/without conversion.
# `qs2_available` is checked once at startup to avoid a requireNamespace() call
# on every load.
qs2_available <- requireNamespace("qs2", quietly = TRUE)
# Threads for qs_read. Reads are short bursts, so a small count gives a big
# speedup (0.9s -> 0.2s on a 230MB object) without meaningfully oversubscribing
# cores shared by the other Shiny workers. Default 4, but never more than the
# available cores; override with QS2_READ_THREADS (also capped to cores).
qs2_read_threads <- {
  cores <- parallel::detectCores()
  if (is.na(cores) || cores < 1) cores <- 1L
  n <- suppressWarnings(as.integer(Sys.getenv("QS2_READ_THREADS")))
  if (is.na(n) || n < 1) min(4L, cores) else min(n, cores)
}
read_object <- function(path) {
  if (length(path) != 1 || is.na(path) || !nzchar(path)) {
    stop("read_object(): 'path' must be a single non-empty file path; got ",
         deparse(path))
  }
  if (qs2_available) {
    qs2_path <- sub("\\.rds$", ".qs2", path, ignore.case = TRUE)
    if (qs2_path != path && file.exists(qs2_path)) {
      cat("[read_object] qs2:", qs2_path, "\n", file = stderr())
      return(qs2::qs_read(qs2_path, nthreads = qs2_read_threads))
    }
  }
  cat("[read_object] rds:", path, "\n", file = stderr())
  readRDS(path)
}

cat("Working directory is:", getwd(), "\n", file = stderr())

# Default comparison metadata for datasets without explicit configuration
DEFAULT_COMPARISON_TYPE <- "disease"
DEFAULT_COMPARISON_LABEL <- "SSc vs Healthy"




dataset_meta <- read.delim(
  file.path("config", "datasets.tsv"),
  sep = "\t",
  stringsAsFactors = FALSE,
  comment.char = "#"
)
dataset_meta$has_scrna <- as.logical(dataset_meta$has_scrna)
dataset_meta$has_spatial <- as.logical(dataset_meta$has_spatial)

# Load dataset details from JSON, if available, and build dynamic data structures.
# This allows for complex configurations with multiple data levels and file paths per dataset.
# If a dataset is not found in the JSON file, it defaults to a minimal configuration.

library(jsonlite)
dataset_details_path <- file.path("config", "dataset_details.json")
dataset_details <- if (file.exists(dataset_details_path)) {
  fromJSON(dataset_details_path)
} else {
  list()
}

# Initialize lists
dataset_files <- list()
data_level_choices <- list()
dataset_comparison_type <- list()
dataset_comparison_label <- list()
dataset_metadata_file <- list()
dataset_comparison_column <- list()

# Populate lists based on dataset_meta and dataset_details
for (i in 1:nrow(dataset_meta)) {
  dataset_id <- dataset_meta$id[i]

  if (dataset_id %in% names(dataset_details)) {
    # If details are available in the JSON file, use them to build the data structures.
    dataset_files[[dataset_id]] <- dataset_details[[dataset_id]]$files
    data_level_choices[[dataset_id]] <- unlist(dataset_details[[dataset_id]]$data_levels)
    dataset_comparison_type[[dataset_id]] <- dataset_details[[dataset_id]]$comparison_type %||% DEFAULT_COMPARISON_TYPE
    dataset_comparison_label[[dataset_id]] <- dataset_details[[dataset_id]]$comparison_label %||% DEFAULT_COMPARISON_LABEL
    dataset_metadata_file[[dataset_id]] <- dataset_details[[dataset_id]]$files$meta %||% NULL
    dataset_comparison_column[[dataset_id]] <- dataset_details[[dataset_id]]$comparison_column %||% "Disease"
  } else {
    # If no details are found for the dataset in the JSON file, fall back to a minimal default structure.
    # This ensures that all datasets listed in datasets.tsv are available in the app, even without a detailed configuration.
    dataset_files[[dataset_id]] <- list(full = list(seurat = dataset_meta$file_path[i]))
    data_level_choices[[dataset_id]] <- c("Full" = "full")
    dataset_comparison_type[[dataset_id]] <- DEFAULT_COMPARISON_TYPE
    dataset_comparison_label[[dataset_id]] <- DEFAULT_COMPARISON_LABEL
    dataset_metadata_file[[dataset_id]] <- NULL
    dataset_comparison_column[[dataset_id]] <- "Disease"
  }
}

# Create a named vector of dataset choices for use in UI dropdowns
dataset_choices <- setNames(dataset_meta$id, dataset_meta$name)

# scRNA-only choices for the Explore sidebar (excludes spatial-only datasets)
scrna_dataset_choices <- setNames(
  dataset_meta$id[dataset_meta$has_scrna],
  dataset_meta$name[dataset_meta$has_scrna]
)