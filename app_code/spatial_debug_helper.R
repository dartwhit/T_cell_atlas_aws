# Spatial Module Debug Helper
# Run this in your R console to check spatial data configuration

library(jsonlite)

# Load configuration
dataset_details <- fromJSON("config/dataset_details.json")
dataset_meta <- read.delim("config/datasets.tsv", sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# Check which datasets have spatial data
cat("\n=== Datasets with Spatial Data ===\n")
for (dataset_id in names(dataset_details)) {
  spatial_file <- dataset_details[[dataset_id]]$files$spatial_seurat
  if (!is.null(spatial_file) && spatial_file != "") {
    cat(sprintf("  %s: %s\n", dataset_id, spatial_file))
    
    # Check if file exists
    full_path <- file.path("data", spatial_file)
    if (file.exists(full_path)) {
      cat(sprintf("    ✓ File exists at: %s\n", full_path))
      
      # Try to load it and check structure
      tryCatch({
        obj <- readRDS(full_path)
        cat(sprintf("    ✓ Object loaded successfully\n"))
        cat(sprintf("    - Images: %s\n", paste(names(obj@images), collapse = ", ")))
        cat(sprintf("    - Assays: %s\n", paste(names(obj@assays), collapse = ", ")))
        cat(sprintf("    - Reductions: %s\n", paste(names(obj@reductions), collapse = ", ")))
        
        # Check for metadata columns
        if ("cluster" %in% colnames(obj@meta.data)) {
          cat("    ✓ 'cluster' column found in metadata\n")
        } else {
          cat("    ✗ 'cluster' column NOT found in metadata\n")
          cat("      Available columns:", paste(colnames(obj@meta.data)[1:min(10, ncol(obj@meta.data))], collapse = ", "), "\n")
        }
        
        if ("Key_Regions" %in% colnames(obj@meta.data)) {
          cat("    ✓ 'Key_Regions' column found in metadata\n")
        } else {
          cat("    ✗ 'Key_Regions' column NOT found in metadata\n")
        }
        
      }, error = function(e) {
        cat(sprintf("    ✗ Error loading file: %s\n", e$message))
      })
    } else {
      cat(sprintf("    ✗ File NOT FOUND at: %s\n", full_path))
    }
  }
}

cat("\n=== Testing dataset_files structure (as built by setup.R) ===\n")
source("setup.R")
studies_with_spatial <- names(dataset_files)[sapply(dataset_files, function(x) !is.null(x$spatial_seurat))]
cat("Studies with spatial data:", paste(studies_with_spatial, collapse = ", "), "\n")

for (study in studies_with_spatial) {
  path <- paste0("data/", dataset_files[[study]][["spatial_seurat"]])
  cat(sprintf("\n  %s:\n", study))
  cat(sprintf("    Path: %s\n", path))
  cat(sprintf("    Exists: %s\n", file.exists(path)))
}
