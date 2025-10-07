#!/usr/bin/env Rscript

# Test script for dynamic dataset configuration
# This script validates that the new dynamic setup works correctly

cat("🧪 Testing Dynamic Dataset Configuration\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Set working directory to app_code for testing
setwd("app_code")

# Source the new setup
tryCatch({
  source("setup.R")
  cat("✅ Setup script loaded successfully\n")
}, error = function(e) {
  cat("❌ Failed to load setup script:", e$message, "\n")
  quit(status = 1)
})

# Test 1: Check dataset_meta is loaded
cat("\n📊 Test 1: Dataset Metadata Loading\n")
if (exists("dataset_meta") && is.data.frame(dataset_meta)) {
  cat("✅ dataset_meta loaded -", nrow(dataset_meta), "datasets\n")
  cat("📋 Datasets:", paste(dataset_meta$id, collapse = ", "), "\n")
} else {
  cat("❌ dataset_meta not properly loaded\n")
}

# Test 2: Check data_level_choices structure
cat("\n🔧 Test 2: Data Level Choices\n")
if (exists("data_level_choices") && is.list(data_level_choices)) {
  cat("✅ data_level_choices created\n")
  for (dataset_id in names(data_level_choices)) {
    levels <- data_level_choices[[dataset_id]]
    cat("  📌", dataset_id, ":", paste(names(levels), collapse = ", "), "\n")
  }
} else {
  cat("❌ data_level_choices not properly created\n")
}

# Test 3: Check dataset_files structure
cat("\n📁 Test 3: Dataset Files Structure\n")
if (exists("dataset_files") && is.list(dataset_files)) {
  cat("✅ dataset_files created\n")
  for (dataset_id in names(dataset_files)[1:2]) {  # Test first 2 datasets
    cat("  📂", dataset_id, "\n")
    for (level in names(dataset_files[[dataset_id]])) {
      if (level %in% c("meta", "spatial_seurat")) {
        cat("    🗂️", level, ":", dataset_files[[dataset_id]][[level]], "\n")
      } else {
        cat("    🗂️", level, "- seurat:", dataset_files[[dataset_id]][[level]][["seurat"]], "\n")
      }
    }
  }
} else {
  cat("❌ dataset_files not properly created\n")
}

# Test 4: Check study levels parsing
cat("\n🎯 Test 4: Study Levels Parsing\n")
test_cases <- c(
  "full:Full;fib:Fibroblasts;immune:Immune cells",
  "full:Full;mye:Myeloid cells;fib:Fibroblasts",
  ""
)

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

for (i in seq_along(test_cases)) {
  result <- parse_study_levels(test_cases[i])
  cat("  📝 Test", i, ":", test_cases[i], "\n")
  cat("    Levels:", paste(result$levels, collapse = ", "), "\n")
  cat("    Labels:", paste(names(result$labels), collapse = ", "), "\n")
}

# Test 5: File existence validation
cat("\n📋 Test 5: File Existence Validation\n")
if (exists("missing_files") && is.list(missing_files)) {
  if (length(missing_files) == 0) {
    cat("✅ All expected files found\n")
  } else {
    cat("⚠️", length(missing_files), "missing files detected\n")
    for (i in 1:min(3, length(missing_files))) {
      cat("  ❌", names(missing_files)[i], "\n")
    }
  }
}

# Test 6: Backward compatibility check
cat("\n🔄 Test 6: Backward Compatibility\n")
required_objects <- c("dataset_files", "data_level_choices", "datasets", "inDir", "DE_dir")
missing_objects <- c()

for (obj in required_objects) {
  if (exists(obj)) {
    cat("✅", obj, "exists\n")
  } else {
    missing_objects <- c(missing_objects, obj)
    cat("❌", obj, "missing\n")
  }
}

if (length(missing_objects) == 0) {
  cat("🎉 All backward compatibility checks passed!\n")
} else {
  cat("⚠️ Some objects missing for backward compatibility\n")
}

# Test 7: Dynamic vs Hardcoded comparison (if we have sample data)
cat("\n⚖️ Test 7: Dynamic Configuration Validation\n")
expected_datasets <- c("tmkmh", "tabib", "gur", "ma", "khanna")
actual_datasets <- names(dataset_files)

cat("Expected datasets:", paste(expected_datasets, collapse = ", "), "\n")
cat("Actual datasets:  ", paste(actual_datasets, collapse = ", "), "\n")

if (all(expected_datasets %in% actual_datasets)) {
  cat("✅ All expected datasets found\n")
} else {
  missing <- expected_datasets[!expected_datasets %in% actual_datasets]
  cat("❌ Missing datasets:", paste(missing, collapse = ", "), "\n")
}

# Test 8: Study levels validation
cat("\n🔍 Test 8: Study Levels Validation\n")
expected_study_levels <- list(
  tmkmh = c("full", "fib", "immune"),
  tabib = c("full", "fib", "immune"),
  gur = c("full", "fib", "immune"),
  ma = c("full", "mye", "fib"),
  khanna = c("full", "fib", "mye")
)

for (dataset_id in names(expected_study_levels)) {
  if (dataset_id %in% names(data_level_choices)) {
    actual_levels <- as.character(data_level_choices[[dataset_id]])
    expected_levels <- expected_study_levels[[dataset_id]]
    
    if (all(expected_levels %in% actual_levels)) {
      cat("✅", dataset_id, "- study levels correct\n")
    } else {
      missing <- expected_levels[!expected_levels %in% actual_levels]
      cat("❌", dataset_id, "- missing study levels:", paste(missing, collapse = ", "), "\n")
    }
  } else {
    cat("❌", dataset_id, "- dataset not found in data_level_choices\n")
  }
}

cat("\n🎊 Testing Complete!\n")
cat(paste(rep("=", 50), collapse = ""), "\n")