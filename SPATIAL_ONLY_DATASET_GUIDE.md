# Adding Spatial-Only Datasets Guide

## Overview
This guide explains how to add datasets that contain **only spatial data** (no scRNA-seq component) to the T-Cell Atlas AWS application.

## Steps to Add a Spatial-Only Dataset

### 1. Prepare Your Data Files
Ensure your spatial data is in Seurat format and includes:
- **Spatial coordinates** in `@images` slot
- **Gene expression data** in appropriate assay (e.g., `SCT`, `RNA`)
- **Metadata** with sample information
- **Optional**: Cluster annotations, cell type predictions

### 2. Update datasets.tsv
Add a new row to `app_code/config/datasets.tsv` with these specifications:

```tsv
spatial_only_study	Spatial Only Study 2024	Spatial	25000	imgs/spatial_study.png	Pure spatial transcriptomics dataset	FALSE	TRUE	spatial:Spatial Analysis	SpatialStudy	SpatialStudy_processed_data.rds	SpatialStudy_metadata.txt	2024-11-15	NA
```

**Key Parameters for Spatial-Only Datasets:**
- `has_scrna`: **FALSE** (no single-cell RNA-seq data)  
- `has_spatial`: **TRUE** (has spatial data)
- `assay`: **Spatial** (or **RNA** if using standard RNA assay)
- `study_levels`: **spatial:Spatial Analysis** (or custom levels if you have subsets)
- `spatial_file`: **YourData.rds** (path to your spatial Seurat object)

### 3. Upload Data Files
Place your files in the appropriate directories:

```
app_code/
├── data/
│   ├── SpatialStudy_processed_data.rds    # Main spatial Seurat object
│   └── SpatialStudy_metadata.txt          # Optional metadata file
└── www/
    └── imgs/
        └── spatial_study.png               # Dataset preview image
```

### 4. File Naming Convention
For spatial-only datasets, the system expects:
- **Main file**: `{file_base_pattern}_{spatial_identifier}.rds`
- **Metadata**: `{file_base_pattern}_metadata.txt` (optional)

### 5. Dataset Configuration Example

```tsv
id: spatial_only_study
name: Spatial Only Study 2024  
assay: Spatial
n_cells: 25000
image: imgs/spatial_study.png
desc: Pure spatial transcriptomics dataset without scRNA-seq component
has_scrna: FALSE
has_spatial: TRUE  
study_levels: spatial:Spatial Analysis
file_base_pattern: SpatialStudy
spatial_file: SpatialStudy_processed_data.rds
meta_file: SpatialStudy_metadata.txt
date_pattern: 2024-11-15
alt_patterns: NA
```

## How It Works

### 1. **Dataset Gallery**
- Spatial-only datasets will show **"Explore Spatial"** button only
- No scRNA-seq exploration option will be available
- Datasets are tagged with **"Spatial"** badge

### 2. **Navigation** 
- Clicking "Explore Spatial" navigates directly to the Spatial tab
- Study selector automatically updates to show your dataset

### 3. **Spatial Analysis**
- Full spatial visualization capabilities available
- UMAP overview (if dimensional reduction exists)
- Feature expression mapping
- Sample-by-sample analysis
- All standard spatial plot types

## Data Requirements

### Required Seurat Object Structure
```r
# Your spatial Seurat object should have:
spatial_seurat@images          # Spatial coordinates and images
spatial_seurat@assays          # Gene expression (RNA/SCT)
spatial_seurat@meta.data       # Cell/spot metadata
spatial_seurat@reductions      # Optional: UMAP/PCA for overview plots
```

### Recommended Metadata Columns
```r
# Useful metadata columns:
- orig.ident       # Sample identifier
- SampleID         # Sample ID  
- cluster          # Cluster assignments
- predicted_labels # Cell type predictions (if available)
```

## Testing Your Dataset

### 1. **Validate Configuration**
```bash
# Run the test suite
Rscript tests/test_dataset_config.R
```

### 2. **Check Data Loading**
```r
# Test in R console
source("app_code/setup.R")
spatial_datasets <- dataset_meta[dataset_meta$has_spatial, ]
print(spatial_datasets$id)  # Should include your dataset
```

### 3. **Manual Testing**
- Start the Shiny app
- Navigate to Datasets tab
- Look for your spatial-only dataset
- Click "Explore Spatial" 
- Verify spatial plots load correctly

## Troubleshooting

### Common Issues

1. **Dataset Not Appearing**
   - Check TSV syntax (tabs, not spaces)
   - Verify `has_spatial = TRUE`
   - Ensure file paths are correct

2. **Spatial Plots Not Loading**
   - Verify Seurat object has `@images` slot populated
   - Check that spatial coordinates exist
   - Ensure gene expression data is in expected assay

3. **No Features Available**
   - Verify gene names in rownames of assays
   - Check assay names (`RNA`, `SCT`, etc.)
   - Ensure expression data is properly normalized

### File Path Debugging
```r
# Check if spatial file is detected
source("app_code/setup.R")
dataset_files[["your_dataset_id"]][["spatial_seurat"]]
```

## Advanced Configuration

### Multiple Spatial Samples
If your dataset has multiple spatial samples, ensure:
```r
# Multiple images in Seurat object
names(spatial_seurat@images)  # Should list all sample names
```

### Custom Study Levels
For spatial datasets with subsets:
```tsv
study_levels: spatial:All Samples;sample1:Sample 1;sample2:Sample 2
```

### Integration with Existing Datasets
To add spatial data to existing scRNA-seq datasets:
```tsv
has_scrna: TRUE
has_spatial: TRUE
study_levels: full:Full;fib:Fibroblasts;spatial:Spatial Analysis
```

## Example Files Structure

```
your_spatial_dataset/
├── processed_data.rds           # Main Seurat object
├── metadata.txt                 # Optional metadata  
├── sample_info.csv             # Sample information
└── images/
    ├── sample1_image.jpg       # Spatial images (if separate)
    └── sample2_image.jpg
```

The dynamic configuration system will automatically:
✅ Detect your spatial-only dataset  
✅ Generate appropriate UI components  
✅ Handle file loading and validation  
✅ Enable spatial-specific analysis tools  
✅ Provide download capabilities  

## Summary

Adding spatial-only datasets is now as simple as:
1. **Add one row** to `datasets.tsv`
2. **Upload your data files**
3. **Set `has_spatial=TRUE, has_scrna=FALSE`**
4. **Restart the application**

The dynamic configuration system handles everything else automatically! 🎉