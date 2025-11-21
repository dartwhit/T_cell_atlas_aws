inDir <- "data/"
DE_dir <- "DE_dfs/"
cat("Working directory is:", getwd(), "\n", file = stderr())




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

# Populate lists based on dataset_meta and dataset_details
for (i in 1:nrow(dataset_meta)) {
  dataset_id <- dataset_meta$id[i]
  
  if (dataset_id %in% names(dataset_details)) {
    # If details are available in the JSON file, use them to build the data structures.
    dataset_files[[dataset_id]] <- dataset_details[[dataset_id]]$files
    data_level_choices[[dataset_id]] <- unlist(dataset_details[[dataset_id]]$data_levels)
  } else {
    # If no details are found for the dataset in the JSON file, fall back to a minimal default structure.
    # This ensures that all datasets listed in datasets.tsv are available in the app, even without a detailed configuration.
    dataset_files[[dataset_id]] <- list(full = list(seurat = dataset_meta$file_path[i]))
    data_level_choices[[dataset_id]] <- c("Full" = "full")
  }
}

# Create a named vector of dataset choices for use in UI dropdowns
dataset_choices <- setNames(dataset_meta$id, dataset_meta$name)