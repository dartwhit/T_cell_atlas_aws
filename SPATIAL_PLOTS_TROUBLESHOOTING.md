# Spatial Plots Troubleshooting Guide

## Issue Summary
Spatial plots work in `test_spat.R` but fail in the deployed UI with sample-specific spatial plots not showing and erroring out.

## Changes Made

### 1. Fixed UI Structure (`ui.R`)
**Problem**: The sidebar controls were structured correctly to match the test app pattern.

**Fix Applied**: Verified the layout_sidebar structure matches test_spat.R pattern where controls are namespaced with `NS("sp1", ...)`.

### 2. Added Error Handling (`modules/spatial_unit.R`)
**Problem**: No error messages when plots fail, making debugging impossible.

**Fixes Applied**:
- Added `tryCatch` blocks to all plot rendering functions (UMAP, FeaturePlot, SpatialDimPlot, SpatialFeaturePlot)
- Error messages now display in the plot area instead of silent failure
- Errors are logged to stderr for server-side debugging

### 3. Added Diagnostic Logging
**Locations**:
- `modules/spatial_unit.R`: Logs when loading Seurat object, reports images/assays available
- `server.R`: Logs spatial data path construction and file existence checks

**How to use**: Check your Shiny app logs (stderr) for messages like:
```
Loading spatial object from path: data/2025-07-07_MaSSc_Visium...
✓ Spatial object loaded successfully
  - Images: MaSSc_1, MaSSc_2
  - Assays: RNA, SCT, Spatial
```

### 4. Created Debug Helper Script
**File**: `spatial_debug_helper.R`

**Usage**:
```r
# Navigate to app_code directory
setwd("app_code")
source("spatial_debug_helper.R")
```

This will:
- List all datasets with spatial data
- Check if spatial RDS files exist
- Load each spatial object and verify structure
- Check for expected metadata columns ("cluster", "Key_Regions")
- Show available metadata columns if expected ones are missing

## Common Issues & Solutions

### Issue 1: "Error: Column 'cluster' not found"
**Cause**: Spatial Seurat object metadata doesn't have a column named "cluster"

**Solution**:
1. Run `spatial_debug_helper.R` to see available columns
2. Either:
   - Option A: Update your Seurat object to rename the appropriate column to "cluster"
   - Option B: Update `ui.R` line 254-257 to use the actual column name

Example for Option B:
```r
# If your actual column is "seurat_clusters" instead of "cluster"
selectInput(NS("sp1", "group_by"), 
            "Group by (metadata)", 
            choices = c("Seurat cluster" = "seurat_clusters",  # Changed here
            "Key regions" =  "Key_Regions"),
          selected = "seurat_clusters")  # And here
```

### Issue 2: "Error: Assay 'SCT' not found"
**Cause**: The spatial object doesn't have an SCT assay (line 77, 165 in spatial_unit.R)

**Solution**: 
If your spatial object uses "Spatial" or another assay instead:
```r
# In spatial_unit.R, replace "SCT" with your actual assay name
DefaultAssay(plot_obj) <- "Spatial"  # or "RNA" or whatever you have
```

### Issue 3: "Images not showing"
**Cause**: Image names in the Seurat object don't match sample names

**Debug**:
The debug helper will show: `Images: MaSSc_1, MaSSc_2`
These names must match what you select in the "Select samples to view" checkbox.

### Issue 4: File path issues
**Symptoms**: "File not found" errors in logs

**Check**:
1. In deployed app, verify `inDir` is set correctly in `setup.R` (currently "data/")
2. Docker volume mounts - ensure data directory is mounted
3. File permissions - ensure R can read the .rds files

**Docker-specific**: If running in Docker, check:
```bash
docker exec <container-name> ls -la /srv/shiny-server/data/
```

### Issue 5: Sample-specific plots fail but main plots work
**Cause**: Usually due to:
- Image-specific issues in the Seurat object
- Mismatch between image names and sample IDs
- Point size calculation failing for some images

**Debug**: Check the stderr logs for messages like:
```
Error in SpatialDimPlot for sample MaSSc_1: <error message>
```

## Testing Workflow

### 1. Test Locally First
```r
# Run the test app
source("modules/test_spat.R")
shinyApp(test_ui, test_server)
```

### 2. Run Diagnostics
```r
source("spatial_debug_helper.R")
```

### 3. Check Deployed App Logs
Monitor stderr when loading spatial data:
```r
# In server.R, already added:
cat("Spatial data requested for study:", study, "\n", file = stderr())
```

### 4. Enable Shiny Trace (Temporarily)
In `server.R`, the line exists but is set to FALSE:
```r
options(shiny.trace = TRUE)  # Change to TRUE for detailed debugging
```

## Key Differences: Test vs Production

| Aspect | test_spat.R | Production (app.R) |
|--------|-------------|-------------------|
| Data source | Hardcoded path | Dynamic from JSON config |
| Study selection | None (single dataset) | Dropdown selector |
| Input controls | In parent UI with NS() | Same structure ✓ |
| Error handling | None | Now added ✓ |
| Logging | None | Now added ✓ |

## Metadata Column Requirements

The spatial module expects these exact column names in `obj@meta.data`:

1. **For grouping dropdown**: 
   - `cluster` (mapped to "Seurat cluster")
   - `Key_Regions` (mapped to "Key regions")

2. **Must exist in the Seurat object**:
   - At least one reduction (UMAP or umap)
   - SCT assay (or modify code to use different assay)
   - At least one image in `obj@images`

## Next Steps

1. **Run the debug helper** to identify exact issue
2. **Check stderr logs** in your deployed app
3. **Verify file paths** and data mounting (especially in Docker)
4. **Check metadata columns** match expectations
5. **Test with simplified data** if issues persist

## Contact/Notes
- Original working test: `modules/test_spat.R`
- Main spatial module: `modules/spatial_unit.R`
- UI integration: `ui.R` lines 249-270
- Server integration: `server.R` lines 56-74
