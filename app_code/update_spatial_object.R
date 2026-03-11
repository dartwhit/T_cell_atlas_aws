# Script to update spatial Seurat object for compatibility
# Run this locally before deploying to fix spatial image alignment issues
#
# Background: Seurat v5 changed how spatial image data is stored.
# Objects created in Seurat v4.x need to be updated before spatial plotting works in v5.
#
# This issue affects:
# - Seurat v4.0 - v4.3 objects used with Seurat v5.0+
# - Spatial images may not align correctly without updating
#
# Solution: Run UpdateSeuratObject() and re-save the RDS file

library(Seurat)

# Configuration
input_file <- "data/2025-07-07_MaSSc_Visium_PRECAST_SingleCellPredicted_RegionsNamed_CARD.rds"
output_file <- "data/2025-07-07_MaSSc_Visium_PRECAST_SingleCellPredicted_RegionsNamed_CARD_updated.rds"

cat("=== Updating Spatial Seurat Object ===\n\n")
cat("Current Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n\n")

# Check if file exists
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Load the object
cat("Loading Seurat object...\n")
obj <- readRDS(input_file)

cat("Original object info:\n")
cat("  - Images:", paste(names(obj@images), collapse = ", "), "\n")
cat("  - Assays:", paste(names(obj@assays), collapse = ", "), "\n")
cat("  - Seurat version (stored):", as.character(slot(obj, "version")), "\n\n")

# Update the object
cat("Updating Seurat object...\n")
obj_updated <- UpdateSeuratObject(obj)

cat("\nUpdated object info:\n")
cat("  - Images:", paste(names(obj_updated@images), collapse = ", "), "\n")
cat("  - Assays:", paste(names(obj_updated@assays), collapse = ", "), "\n")
cat("  - Seurat version (updated):", as.character(slot(obj_updated, "version")), "\n\n")

# Verify images are accessible
cat("Verifying spatial images...\n")
for (img_name in names(obj_updated@images)) {
  img <- obj_updated@images[[img_name]]
  cat("  ✓", img_name, "- spots:", nrow(img@coordinates), "\n")
}

# Save the updated object
cat("\nSaving updated object to:", output_file, "\n")
saveRDS(obj_updated, file = output_file)

cat("\n✓ Update complete!\n")
cat("\nNext steps:\n")
cat("1. Test the updated object locally:\n")
cat("   obj <- readRDS('", output_file, "')\n", sep = "")
cat("   SpatialDimPlot(obj, images = 'A')\n")
cat("2. If tests pass, replace the original file:\n")
cat("   file.rename('", output_file, "', '", input_file, "')\n", sep = "")
cat("3. Commit and deploy the updated RDS file\n")
