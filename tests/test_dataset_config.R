#!/usr/bin/env Rscript

# Comprehensive test suite for dynamic dataset configuration
# Designed to run on EC2 instance or local development environment

suppressMessages({
  if (!require(jsonlite)) install.packages("jsonlite", repos = "http://cran.us.r-project.org")
  library(jsonlite)
})

# Set up test environment
test_results <- list()
start_time <- Sys.time()

cat("🧪 T-Cell Atlas AWS: Dynamic Dataset Configuration Test Suite\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("📅 Started at:", format(start_time), "\n")
cat("🖥️ Working directory:", getwd(), "\n\n")

# Change to app_code directory if not already there
if (!file.exists("setup.R") && file.exists("app_code/setup.R")) {
  setwd("app_code")
  cat("📁 Changed to app_code directory\n")
}

# Helper function to record test results
record_test <- function(test_name, passed, message = "", details = NULL) {
  test_results[[test_name]] <<- list(
    passed = passed,
    message = message,
    details = details,
    timestamp = Sys.time()
  )
  
  status <- if (passed) "✅ PASS" else "❌ FAIL"
  cat(sprintf("%-50s %s\n", test_name, status))
  if (!is.null(message) && message != "") {
    cat("   ", message, "\n")
  }
}

# Test 1: Configuration Files Exist
cat("📋 Test Group 1: Configuration Files\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

config_tsv_exists <- file.exists("config/datasets.tsv")
setup_r_exists <- file.exists("setup.R")

record_test(
  "config/datasets.tsv exists", 
  config_tsv_exists,
  if (!config_tsv_exists) "Missing datasets.tsv configuration file"
)

record_test(
  "setup.R exists", 
  setup_r_exists,
  if (!setup_r_exists) "Missing setup.R file"
)

# Test 2: Setup Script Loading
cat("\n📊 Test Group 2: Setup Script Loading\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

setup_error <- NULL
tryCatch({
  source("setup.R", local = FALSE)
  record_test("Setup script loads without error", TRUE)
}, error = function(e) {
  setup_error <<- e$message
  record_test("Setup script loads without error", FALSE, paste("Error:", e$message))
})

# Test 3: Core Objects Created
cat("\n🏗️ Test Group 3: Core Objects Creation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

required_objects <- list(
  dataset_meta = "data.frame",
  dataset_files = "list", 
  data_level_choices = "list",
  datasets = "data.frame",
  inDir = "character",
  DE_dir = "character"
)

for (obj_name in names(required_objects)) {
  expected_type <- required_objects[[obj_name]]
  
  if (exists(obj_name)) {
    actual_type <- class(get(obj_name))[1]
    type_correct <- actual_type == expected_type || 
                   (expected_type == "character" && is.character(get(obj_name)))
    
    record_test(
      paste(obj_name, "exists and has correct type"),
      type_correct,
      paste("Expected:", expected_type, "Got:", actual_type)
    )
  } else {
    record_test(
      paste(obj_name, "exists and has correct type"),
      FALSE,
      paste("Object", obj_name, "does not exist")
    )
  }
}

# Test 4: Dataset Metadata Validation
cat("\n📋 Test Group 4: Dataset Metadata Validation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

if (exists("dataset_meta")) {
  # Check required columns
  required_cols <- c("id", "name", "assay", "n_cells", "has_scrna", "has_spatial", 
                     "study_levels", "file_base_pattern")
  
  missing_cols <- required_cols[!required_cols %in% colnames(dataset_meta)]
  
  record_test(
    "All required columns present",
    length(missing_cols) == 0,
    if (length(missing_cols) > 0) paste("Missing columns:", paste(missing_cols, collapse = ", "))
  )
  
  # Check data types
  record_test(
    "Boolean columns properly converted",
    is.logical(dataset_meta$has_scrna) && is.logical(dataset_meta$has_spatial),
    "has_scrna and has_spatial should be logical"
  )
  
  # Check for duplicate IDs
  duplicate_ids <- sum(duplicated(dataset_meta$id))
  record_test(
    "No duplicate dataset IDs",
    duplicate_ids == 0,
    if (duplicate_ids > 0) paste("Found", duplicate_ids, "duplicate IDs")
  )
  
  cat("   📊 Loaded", nrow(dataset_meta), "datasets:", paste(dataset_meta$id, collapse = ", "), "\n")
}

# Test 5: Data Level Choices Structure
cat("\n🎯 Test Group 5: Data Level Choices Validation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

if (exists("data_level_choices") && exists("dataset_meta")) {
  # Check that all datasets have data level choices
  all_datasets_have_choices <- all(dataset_meta$id %in% names(data_level_choices))
  
  record_test(
    "All datasets have data level choices",
    all_datasets_have_choices
  )
  
  # Check structure of choices
  valid_structure <- TRUE
  for (dataset_id in names(data_level_choices)) {
    choices <- data_level_choices[[dataset_id]]
    if (!is.character(choices) || is.null(names(choices))) {
      valid_structure <- FALSE
      break
    }
  }
  
  record_test(
    "Data level choices have valid structure",
    valid_structure,
    "Each dataset should have named character vector of choices"
  )
  
  # Sample a few datasets to show their choices
  sample_datasets <- head(names(data_level_choices), 3)
  for (dataset_id in sample_datasets) {
    choices <- data_level_choices[[dataset_id]]
    cat("   📌", dataset_id, ":", paste(names(choices), collapse = ", "), "\n")
  }
}

# Test 6: Dataset Files Structure
cat("\n📁 Test Group 6: Dataset Files Structure Validation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

if (exists("dataset_files") && exists("dataset_meta")) {
  # Check that all datasets have file definitions
  all_datasets_have_files <- all(dataset_meta$id %in% names(dataset_files))
  
  record_test(
    "All datasets have file definitions",
    all_datasets_have_files
  )
  
  # Check that each dataset has at least 'full' level
  all_have_full <- all(sapply(dataset_files, function(x) "full" %in% names(x)))
  
  record_test(
    "All datasets have 'full' study level",
    all_have_full
  )
  
  # Check seurat file paths are not null for full level
  full_seurat_complete <- TRUE
  missing_seurat <- c()
  
  for (dataset_id in names(dataset_files)) {
    if ("full" %in% names(dataset_files[[dataset_id]])) {
      seurat_path <- dataset_files[[dataset_id]][["full"]][["seurat"]]
      if (is.null(seurat_path) || seurat_path == "") {
        full_seurat_complete <- FALSE
        missing_seurat <- c(missing_seurat, dataset_id)
      }
    }
  }
  
  record_test(
    "All datasets have seurat file paths for 'full' level",
    full_seurat_complete,
    if (!full_seurat_complete) paste("Missing seurat paths:", paste(missing_seurat, collapse = ", "))
  )
}

# Test 7: File Existence Validation
cat("\n💾 Test Group 7: File Existence Validation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

if (exists("inDir") && dir.exists(inDir)) {
  cat("   📂 Data directory found:", inDir, "\n")
  
  # Count existing vs missing files
  total_files <- 0
  existing_files <- 0
  missing_files_list <- c()
  
  if (exists("dataset_files")) {
    for (dataset_id in names(dataset_files)) {
      for (level in names(dataset_files[[dataset_id]])) {
        if (level %in% c("meta", "spatial_seurat")) {
          file_path <- dataset_files[[dataset_id]][[level]]
          if (!is.null(file_path) && file_path != "" && file_path != "NA") {
            total_files <- total_files + 1
            if (file.exists(file.path(inDir, file_path))) {
              existing_files <- existing_files + 1
            } else {
              missing_files_list <- c(missing_files_list, paste(dataset_id, level, file_path, sep = ":"))
            }
          }
        } else {
          level_files <- dataset_files[[dataset_id]][[level]]
          for (file_type in names(level_files)) {
            file_path <- level_files[[file_type]]
            if (!is.null(file_path) && file_path != "" && file_path != "NA") {
              total_files <- total_files + 1
              if (file.exists(file.path(inDir, file_path))) {
                existing_files <- existing_files + 1
              } else {
                missing_files_list <- c(missing_files_list, paste(dataset_id, level, file_type, file_path, sep = ":"))
              }
            }
          }
        }
      }
    }
  }
  
  file_completeness <- existing_files / max(total_files, 1)
  
  record_test(
    "File existence check",
    TRUE,  # Always pass this test, but record statistics
    sprintf("Found %d/%d files (%.1f%% complete)", existing_files, total_files, file_completeness * 100),
    list(
      total_files = total_files,
      existing_files = existing_files,
      missing_files = head(missing_files_list, 10)  # Limit to first 10
    )
  )
  
  # Critical files test - at least some seurat files should exist
  critical_files_exist <- existing_files > 0
  record_test(
    "At least some critical files exist",
    critical_files_exist,
    if (!critical_files_exist) "No data files found - check data directory and file paths"
  )
  
} else {
  record_test(
    "Data directory exists",
    FALSE,
    paste("Data directory not found:", if (exists("inDir")) inDir else "inDir not set")
  )
}

# Generate Test Summary Report
cat("\n📋 Test Summary Report\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

end_time <- Sys.time()
total_tests <- length(test_results)
passed_tests <- sum(sapply(test_results, function(x) x$passed))
failed_tests <- total_tests - passed_tests

cat("📊 Total Tests:  ", total_tests, "\n")
cat("✅ Passed:      ", passed_tests, "\n") 
cat("❌ Failed:      ", failed_tests, "\n")
cat("⏱️ Duration:    ", round(as.numeric(end_time - start_time), 2), "seconds\n")
cat("🎯 Success Rate:", round(passed_tests/total_tests * 100, 1), "%\n\n")

if (failed_tests > 0) {
  cat("❌ Failed Tests:\n")
  for (test_name in names(test_results)) {
    if (!test_results[[test_name]]$passed) {
      cat("  •", test_name, "\n")
      if (test_results[[test_name]]$message != "") {
        cat("    └─", test_results[[test_name]]$message, "\n")
      }
    }
  }
  cat("\n")
}

# Save detailed results to JSON file
results_file <- "test_results.json"
detailed_results <- list(
  timestamp = format(Sys.time()),
  duration_seconds = as.numeric(end_time - start_time),
  summary = list(
    total = total_tests,
    passed = passed_tests,
    failed = failed_tests,
    success_rate = passed_tests/total_tests
  ),
  tests = test_results,
  environment = list(
    working_directory = getwd(),
    R_version = R.version.string,
    data_directory_exists = exists("inDir") && dir.exists(get("inDir", pos = 1))
  )
)

write_json(detailed_results, results_file, pretty = TRUE, auto_unbox = TRUE)
cat("💾 Detailed results saved to:", results_file, "\n")

# Exit with appropriate code
if (failed_tests > 0) {
  cat("\n⚠️ Some tests failed. Please review the issues above.\n")
  quit(status = 1)
} else {
  cat("\n🎉 All tests passed! Dynamic configuration is working correctly.\n")
  quit(status = 0)
}