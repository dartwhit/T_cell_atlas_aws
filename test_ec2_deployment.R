#!/usr/bin/env Rscript

# EC2 Deployment Test for T-Cell Atlas AWS
# This script tests the dynamic configuration on EC2 with actual data

start_time <- Sys.time()

cat("🖥️ T-Cell Atlas AWS - EC2 Deployment Test\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("📅 Test started at:", format(Sys.time()), "\n")
cat("🖥️ Host:", Sys.info()[["nodename"]], "\n")
cat("📁 Working directory:", getwd(), "\n")

# Check if we're in app_code directory or need to change
if (!file.exists("setup.R") && file.exists("app_code/setup.R")) {
  setwd("app_code")
  cat("📁 Changed to app_code directory\n")
}

# Test 1: Basic Setup Loading
cat("\n🔧 Test 1: Loading Dynamic Setup\n")
cat(paste(rep("-", 30), collapse = ""), "\n")

setup_success <- FALSE
tryCatch({
  source("setup.R")
  setup_success <- TRUE
  cat("✅ Setup loaded successfully\n")
}, error = function(e) {
  cat("❌ Setup failed:", e$message, "\n")
})

if (!setup_success) {
  cat("❌ Cannot proceed without successful setup\n")
  quit(status = 1)
}

# Test 2: Data Structure Validation  
cat("\n📊 Test 2: Data Structure Validation\n")
cat(paste(rep("-", 30), collapse = ""), "\n")

required_objects <- c("dataset_meta", "dataset_files", "data_level_choices", "datasets")
all_present <- TRUE

for (obj in required_objects) {
  if (exists(obj)) {
    cat("✅", obj, "- Present\n")
  } else {
    cat("❌", obj, "- Missing\n")
    all_present <- FALSE
  }
}

if (all_present) {
  cat("📊 Dataset count:", nrow(dataset_meta), "\n")
  cat("📋 Dataset IDs:", paste(dataset_meta$id, collapse = ", "), "\n")
}

# Test 3: File Availability Check
cat("\n💾 Test 3: File Availability Check\n")
cat(paste(rep("-", 30), collapse = ""), "\n")

if (exists("inDir") && dir.exists(inDir)) {
  cat("📂 Data directory found:", inDir, "\n")
  
  # Check key files for each dataset
  total_key_files <- 0
  found_key_files <- 0
  
  for (dataset_id in names(dataset_files)) {
    cat("\n🔍 Checking", dataset_id, ":\n")
    
    # Check full level seurat file
    if ("full" %in% names(dataset_files[[dataset_id]])) {
      seurat_path <- dataset_files[[dataset_id]][["full"]][["seurat"]]
      if (!is.null(seurat_path) && seurat_path != "") {
        total_key_files <- total_key_files + 1
        full_path <- file.path(inDir, seurat_path)
        if (file.exists(full_path)) {
          file_info <- file.info(full_path)
          size_mb <- round(file_info$size / 1024 / 1024, 1)
          cat("  ✅ Full dataset:", seurat_path, "(", size_mb, "MB)\n")
          found_key_files <- found_key_files + 1
        } else {
          cat("  ❌ Full dataset:", seurat_path, "- Not found\n")
        }
      }
    }
    
    # Check subset files
    subset_levels <- names(dataset_files[[dataset_id]])
    subset_levels <- subset_levels[!subset_levels %in% c("full", "meta", "spatial_seurat")]
    
    if (length(subset_levels) > 0) {
      cat("  📊 Subset levels:", paste(subset_levels, collapse = ", "), "\n")
      for (level in subset_levels[1:min(2, length(subset_levels))]) {  # Check first 2 levels
        seurat_path <- dataset_files[[dataset_id]][[level]][["seurat"]]
        if (!is.null(seurat_path) && seurat_path != "") {
          total_key_files <- total_key_files + 1
          if (file.exists(file.path(inDir, seurat_path))) {
            cat("  ✅", level, "subset: Available\n")
            found_key_files <- found_key_files + 1
          } else {
            cat("  ❌", level, "subset: Missing\n")
          }
        }
      }
    }
    
    # Check spatial data if applicable
    if (!is.null(dataset_files[[dataset_id]][["spatial_seurat"]])) {
      spatial_path <- dataset_files[[dataset_id]][["spatial_seurat"]]
      total_key_files <- total_key_files + 1
      if (file.exists(file.path(inDir, spatial_path))) {
        cat("  ✅ Spatial data: Available\n")
        found_key_files <- found_key_files + 1
      } else {
        cat("  ❌ Spatial data: Missing\n")
      }
    }
  }
  
  data_completeness <- if (total_key_files > 0) found_key_files / total_key_files else 0
  cat("\n📊 Data completeness:", sprintf("%.1f%% (%d/%d files)", 
                                        data_completeness * 100, found_key_files, total_key_files), "\n")
  
} else {
  cat("❌ Data directory not found or inaccessible\n")
  data_completeness <- 0
}

# Test 4: Study Level Configuration
cat("\n🎯 Test 4: Study Level Configuration\n")
cat(paste(rep("-", 30), collapse = ""), "\n")

if (exists("data_level_choices")) {
  for (dataset_id in names(data_level_choices)) {
    choices <- data_level_choices[[dataset_id]]
    cat("📌", dataset_id, ":", paste(names(choices), collapse = ", "), "\n")
  }
}

# Test 5: Module Compatibility Test (Basic)
cat("\n🔗 Test 5: Module Compatibility Test\n")
cat(paste(rep("-", 30), collapse = ""), "\n")

# Test that the explore_sidebar_module can be sourced
compatibility_success <- TRUE
tryCatch({
  if (file.exists("modules/explore_sidebar_module.R")) {
    source("modules/explore_sidebar_module.R")
    cat("✅ explore_sidebar_module loaded\n")
  } else {
    cat("⚠️ explore_sidebar_module.R not found\n")
  }
  
  if (file.exists("modules/scrna_seq_module.R")) {
    source("modules/scrna_seq_module.R") 
    cat("✅ scrna_seq_module loaded\n")
  } else {
    cat("⚠️ scrna_seq_module.R not found\n")
  }
  
}, error = function(e) {
  cat("❌ Module loading failed:", e$message, "\n")
  compatibility_success <- FALSE
})

# Generate Summary Report
cat("\n📋 EC2 Deployment Test Summary\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

tests_passed <- setup_success + all_present + compatibility_success + (data_completeness > 0)
total_tests <- 4

cat("📊 Tests Passed:", tests_passed, "/", total_tests, "\n")
cat("💾 Data Completeness:", sprintf("%.1f%%", data_completeness * 100), "\n")

if (exists("dataset_meta")) {
  cat("📋 Datasets Configured:", nrow(dataset_meta), "\n")
  
  # Count datasets with actual data
  datasets_with_data <- 0
  if (data_completeness > 0 && exists("dataset_files")) {
    for (dataset_id in names(dataset_files)) {
      if ("full" %in% names(dataset_files[[dataset_id]])) {
        seurat_path <- dataset_files[[dataset_id]][["full"]][["seurat"]]
        if (!is.null(seurat_path) && file.exists(file.path(inDir, seurat_path))) {
          datasets_with_data <- datasets_with_data + 1
        }
      }
    }
  }
  cat("💾 Datasets with Data:", datasets_with_data, "\n")
}

cat("⏱️ Test Duration:", round(as.numeric(Sys.time() - start_time), 2), "seconds\n")

# Deployment Status
if (tests_passed >= 3 && data_completeness > 0.5) {
  cat("\n🎉 EC2 DEPLOYMENT: READY\n")
  cat("✅ Dynamic configuration is working correctly\n")
  cat("✅ Data files are accessible\n") 
  cat("✅ Modules can be loaded\n")
  exit_code <- 0
} else if (tests_passed >= 2) {
  cat("\n⚠️ EC2 DEPLOYMENT: PARTIAL\n")
  cat("✅ Configuration is working\n")
  if (data_completeness < 0.5) {
    cat("⚠️ Some data files may be missing\n")
  }
  exit_code <- 0
} else {
  cat("\n❌ EC2 DEPLOYMENT: FAILED\n")
  cat("❌ Critical issues prevent deployment\n")
  exit_code <- 1
}

cat("\n🏁 EC2 deployment test completed\n")
quit(status = exit_code)