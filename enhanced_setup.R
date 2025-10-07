inDir <- "data/"
DE_dir <- "DE_dfs/"
cat("Working directory is:", getwd(), "\n", file = stderr())

# Load dataset metadata from enhanced TSV file
dataset_meta <- read.delim(
  file.path("config", "datasets.tsv"),
  sep = "\t",
  stringsAsFactors = FALSE
)

# Convert boolean columns
dataset_meta$has_scrna <- as.logical(dataset_meta$has_scrna)
dataset_meta$has_spatial <- as.logical(dataset_meta$has_spatial)

# Parse study levels with labels
parse_study_levels <- function(study_levels_str) {
  if (is.na(study_levels_str) || study_levels_str == "") {
    return(list(levels = "full", labels = c("Full" = "full")))
  }
  
  pairs <- strsplit(study_levels_str, ";")[[1]]
  levels <- character()
  labels <- character()
  
  for (pair in pairs) {
    parts <- strsplit(pair, ":")[[1]]
    if (length(parts) == 2) {
      level_code <- parts[1]
      level_label <- parts[2]
      levels <- c(levels, level_code)
      labels[level_label] <- level_code
    }
  }
  
  return(list(levels = levels, labels = labels))
}

# Add parsed study levels to metadata
dataset_meta$parsed_levels <- lapply(dataset_meta$study_levels, parse_study_levels)

# Helper function to find actual file based on patterns
find_file <- function(pattern, data_dir = inDir) {
  # First try exact pattern
  if (file.exists(file.path(data_dir, pattern))) {
    return(pattern)
  }
  
  # Try with wildcard substitution (replace * with actual date patterns)
  if (grepl("\\*", pattern)) {
    # Look for files matching the pattern
    pattern_regex <- gsub("\\*", ".*", pattern)
    pattern_regex <- gsub("\\.", "\\\\.", pattern_regex)
    
    files <- list.files(data_dir, pattern = pattern_regex, full.names = FALSE)
    if (length(files) > 0) {
      # Return the most recent file if multiple matches
      return(files[length(files)])
    }
  }
  
  return(NULL)
}

# Helper function to generate file paths with intelligent fallback
generate_file_paths <- function(dataset_row, study_level) {
  dataset_id <- dataset_row$id
  file_base <- dataset_row$file_base_pattern
  date_pattern <- dataset_row$date_pattern
  alt_patterns <- dataset_row$alt_patterns
  
  base_files <- list()
  
  if (study_level == "full") {
    # Try different naming patterns for full dataset
    seurat_candidates <- c()
    
    if (!is.na(date_pattern) && date_pattern != "") {
      seurat_candidates <- c(
        paste0(file_base, "_filtered_reduced_", date_pattern, "_VAM.rds"),
        paste0(file_base, "_full_filtered_", date_pattern, "_VAM.rds")
      )
    }
    
    seurat_candidates <- c(seurat_candidates,
      paste0(file_base, "_full_filtered_VAM.rds"),
      paste0(file_base, "_filtered_VAM.rds")
    )
    
    if (!is.na(alt_patterns) && alt_patterns != "") {
      seurat_candidates <- c(seurat_candidates, alt_patterns)
    }
    
    # Find the first existing file
    seurat_file <- NULL
    for (candidate in seurat_candidates) {
      found_file <- find_file(candidate)
      if (!is.null(found_file)) {
        seurat_file <- found_file
        break
      }
    }
    
    base_files <- list(
      seurat = seurat_file,
      gene_list = paste0(file_base, "_gene_list.rds"),
      DEGs_auto = paste0(file_base, "_auto_DEGs_sig.txt"),
      DEGs_broad = paste0(file_base, "_broad_DEGs_sig.txt"),
      VAM_df = paste0(file_base, "_full_filtered_VAM_top.txt"),
      DE_by_disease_broad = paste0(file_base, "_full_broad_DE_by_disease.txt"),
      DE_by_disease_auto = paste0(file_base, "_full_auto_DE_by_disease.txt"),
      VAM_by_disease_broad = paste0(file_base, "_broad_VAMcdf_DE_by_disease.txt"),
      VAM_by_disease_auto = paste0(file_base, "_auto_VAMcdf_DE_by_disease.txt")
    )
    
    # Alternative VAM_df naming
    alt_vam <- paste0(file_base, "_filtered_reduced_", date_pattern, "_VAM_top.txt")
    if (!is.na(date_pattern) && find_file(alt_vam)) {
      base_files$VAM_df <- alt_vam
    }
    
  } else {
    # Subset-specific files
    base_files <- list(
      seurat = paste0(file_base, "_", study_level, "_filtered_VAM.rds"),
      gene_list = paste0(file_base, "_", study_level, "_gene_list.rds"),
      DEGs = paste0(file_base, "_", study_level, "_DEGs_sig.txt"),
      VAM_df = paste0(file_base, "_", study_level, "_filtered_VAM_top.txt"),
      DE_by_disease = paste0(file_base, "_", study_level, "_DE_by_disease.txt"),
      VAM_by_disease = paste0(file_base, "_", study_level, "_VAMcdf_DE_by_disease.txt")
    )
  }
  
  # Validate files exist and warn about missing ones
  for (file_type in names(base_files)) {
    if (!is.null(base_files[[file_type]]) && !find_file(base_files[[file_type]])) {
      cat("⚠️ Missing file for", dataset_id, study_level, file_type, ":", base_files[[file_type]], "\n", file = stderr())
    }
  }
  
  return(base_files)
}

# Generate dataset_files dynamically
dataset_files <- list()
data_level_choices <- list()

for (i in 1:nrow(dataset_meta)) {
  row <- dataset_meta[i, ]
  dataset_id <- row$id
  parsed_levels <- row$parsed_levels[[1]]
  
  dataset_files[[dataset_id]] <- list()
  
  # Add metadata file if it exists
  if (!is.na(row$meta_file) && row$meta_file != "" && row$meta_file != "NA") {
    dataset_files[[dataset_id]][["meta"]] <- row$meta_file
  }
  
  # Generate files for each study level
  for (level in parsed_levels$levels) {
    dataset_files[[dataset_id]][[level]] <- generate_file_paths(row, level)
  }
  
  # Add spatial data if available
  if (row$has_spatial && !is.na(row$spatial_file) && row$spatial_file != "NA") {
    dataset_files[[dataset_id]][["spatial_seurat"]] <- row$spatial_file
  } else {
    dataset_files[[dataset_id]][["spatial_seurat"]] <- NULL
  }
  
  # Set data level choices
  data_level_choices[[dataset_id]] <- parsed_levels$labels
}

# Create datasets object for sidebar module (backward compatibility)
datasets <- dataset_meta

# Print summary
cat("✅ Enhanced setup completed. Loaded", nrow(dataset_meta), "datasets\n", file = stderr())
cat("📁 Dataset files structure generated for:", paste(names(dataset_files), collapse = ", "), "\n", file = stderr())
cat("🔍 Data level choices:", capture.output(str(data_level_choices)), "\n", file = stderr())