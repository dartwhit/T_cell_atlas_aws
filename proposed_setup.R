inDir <- "data/"
DE_dir <- "DE_dfs/"
cat("Working directory is:", getwd(), "\n", file = stderr())

# Load dataset metadata from single TSV file
dataset_meta <- read.delim(
  file.path("config", "datasets.tsv"),
  sep = "\t",
  stringsAsFactors = FALSE
)

# Convert boolean columns
dataset_meta$has_scrna <- as.logical(dataset_meta$has_scrna)
dataset_meta$has_spatial <- as.logical(dataset_meta$has_spatial)

# Parse study levels into list format
dataset_meta$study_levels_list <- lapply(dataset_meta$study_levels, function(x) {
  if (is.na(x) || x == "") return("full")
  strsplit(x, ",")[[1]]
})

# Helper function to generate file paths based on naming patterns
generate_file_paths <- function(dataset_id, file_pattern, study_level) {
  base_files <- list()
  
  if (study_level == "full") {
    # Full dataset files
    base_files <- list(
      seurat = paste0(file_pattern, "_filtered_reduced_", format(Sys.Date(), "%Y-%m-%d"), "_VAM.rds"),
      gene_list = paste0(file_pattern, "_gene_list.rds"),
      DEGs_auto = paste0(file_pattern, "_auto_DEGs_sig.txt"),
      DEGs_broad = paste0(file_pattern, "_broad_DEGs_sig.txt"),
      VAM_df = paste0(file_pattern, "_full_filtered_VAM_top.txt"),
      DE_by_disease_broad = paste0(file_pattern, "_full_broad_DE_by_disease.txt"),
      DE_by_disease_auto = paste0(file_pattern, "_full_auto_DE_by_disease.txt"),
      VAM_by_disease_broad = paste0(file_pattern, "_broad_VAMcdf_DE_by_disease.txt"),
      VAM_by_disease_auto = paste0(file_pattern, "_auto_VAMcdf_DE_by_disease.txt")
    )
    
    # Try alternative naming patterns if primary doesn't exist
    alt_seurat <- paste0(file_pattern, "_full_filtered_VAM.rds")
    if (file.exists(file.path(inDir, alt_seurat)) && !file.exists(file.path(inDir, base_files$seurat))) {
      base_files$seurat <- alt_seurat
    }
    
  } else {
    # Subset-specific files
    base_files <- list(
      seurat = paste0(file_pattern, "_", study_level, "_filtered_VAM.rds"),
      gene_list = paste0(file_pattern, "_", study_level, "_gene_list.rds"),
      DEGs = paste0(file_pattern, "_", study_level, "_DEGs_sig.txt"),
      VAM_df = paste0(file_pattern, "_", study_level, "_filtered_VAM_top.txt"),
      DE_by_disease = paste0(file_pattern, "_", study_level, "_DE_by_disease.txt"),
      VAM_by_disease = paste0(file_pattern, "_", study_level, "_VAMcdf_DE_by_disease.txt")
    )
  }
  
  return(base_files)
}

# Generate dataset_files dynamically
dataset_files <- list()
for (i in 1:nrow(dataset_meta)) {
  row <- dataset_meta[i, ]
  dataset_id <- row$id
  file_pattern <- row$file_patterns
  study_levels <- row$study_levels_list[[1]]
  
  dataset_files[[dataset_id]] <- list()
  
  # Add metadata file if it exists
  if (!is.na(row$meta_file) && row$meta_file != "") {
    dataset_files[[dataset_id]][["meta"]] <- row$meta_file
  }
  
  # Generate files for each study level
  for (level in study_levels) {
    dataset_files[[dataset_id]][[level]] <- generate_file_paths(dataset_id, file_pattern, level)
  }
  
  # Add spatial data if available
  if (row$has_spatial && !is.na(row$spatial_file) && row$spatial_file != "") {
    dataset_files[[dataset_id]][["spatial_seurat"]] <- row$spatial_file
  } else {
    dataset_files[[dataset_id]][["spatial_seurat"]] <- NULL
  }
}

# Generate data_level_choices dynamically
data_level_choices <- list()
level_labels <- list(
  full = "Full",
  fib = "Fibroblasts", 
  immune = "Immune cells",
  mye = "Myeloid cells"
)

for (i in 1:nrow(dataset_meta)) {
  row <- dataset_meta[i, ]
  dataset_id <- row$id
  study_levels <- row$study_levels_list[[1]]
  
  # Create named vector for choices
  choices <- character()
  for (level in study_levels) {
    label <- ifelse(level %in% names(level_labels), level_labels[[level]], tools::toTitleCase(level))
    choices[label] <- level
  }
  
  data_level_choices[[dataset_id]] <- choices
}

# Create datasets object for sidebar module
datasets <- dataset_meta

cat("✅ Dynamic setup completed. Loaded", nrow(dataset_meta), "datasets\n", file = stderr())
cat("📁 Dataset files structure generated for:", paste(names(dataset_files), collapse = ", "), "\n", file = stderr())