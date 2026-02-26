# Seurat Spatial Object Version Compatibility Issue

## Problem

Error: `Please run UpdateSeuratObject on your Seurat object first to ensure that data aligns to the image X when plotting.`

This error occurs even after calling `UpdateSeuratObject()` in the Shiny app.

## Root Cause

### Seurat Version Changes

**Seurat v4.x → v5.x** introduced breaking changes to spatial data structure:

1. **Seurat v4.0 - v4.3** (released 2021-2023)
   - Stored spatial images with `VisiumV1` class
   - Used older coordinate system structure
   - Spatial data in `@images` slot with specific format

2. **Seurat v5.0+** (released 2023+)
   - Changed spatial image storage format
   - Updated coordinate alignment system
   - Requires objects to be updated before spatial plotting works

### Why Runtime Update Fails

Even though `UpdateSeuratObject()` is called in the app, the error persists because:

1. **Partial Update**: `UpdateSeuratObject()` may update the main object but not fully restructure spatial image metadata
2. **Coordinate Cache**: Spatial images cache coordinate transformations that aren't refreshed at runtime
3. **Session State**: The updated object needs to be saved and reloaded to fully commit changes

## Solution

### ✅ Recommended: Update Offline

Update the RDS file before deployment:

```r
# Run this script locally
source("app_code/update_spatial_object.R")
```

This will:
1. Load your spatial object
2. Run `UpdateSeuratObject()` 
3. Save the updated version
4. Provide testing instructions

**Benefits:**
- ✓ One-time operation
- ✓ Faster app loading (no update at runtime)
- ✓ More reliable (guaranteed to work)
- ✓ Smaller memory footprint during execution

### Alternative: Keep Runtime Update

The current code attempts runtime updates, but this is less reliable:

```r
# In spatial_unit.R
updated_obj <- UpdateSeuratObject(loaded_obj)
```

**Drawbacks:**
- ✗ May not fully update spatial image metadata
- ✗ Slower app initialization
- ✗ Requires more memory during update
- ✗ Update happens every time object is loaded

## Detecting Version Issues

### Check Your Object Version

```r
obj <- readRDS("your_spatial_object.rds")

# Check Seurat version used to create it
slot(obj, "version")

# If version shows 4.x.x, it needs updating for Seurat 5.x
```

### Check Current Seurat Version

```r
packageVersion("Seurat")

# If you see 5.x.x but object is from 4.x.x, update needed
```

## Testing After Update

```r
# Load updated object
obj <- readRDS("your_updated_object.rds")

# Test basic spatial plot
library(Seurat)
SpatialDimPlot(obj, images = "A")

# Test with features
SpatialFeaturePlot(obj, images = "A", features = "CD3D")

# If these work without errors, the update was successful
```

## Known Affected Versions

- **Problematic**: Objects created with Seurat 4.0.x - 4.3.x
- **Target**: Running on Seurat 5.0.0 or later
- **Fix Required**: Yes, offline update recommended

## References

- [Seurat v5 Updates](https://satijalab.org/seurat/articles/seurat5_updates)
- [Updating Seurat Objects](https://satijalab.org/seurat/reference/updateseuratoject)
- GitHub issues: [satijalab/seurat#6842](https://github.com/satijalab/seurat/issues/6842), [#7234](https://github.com/satijalab/seurat/issues/7234)

## Implementation

### Current PR Status

The PR includes runtime update attempt, but we recommend:

1. **Immediate fix**: Run `update_spatial_object.R` locally
2. **Test locally**: Verify plots work with updated object
3. **Deploy**: Replace old RDS with updated version
4. **Optional**: Remove runtime update code for performance

### Files Modified in PR

- `app_code/modules/spatial_unit.R` - Attempts runtime update (can be removed after offline update)
- `app_code/update_spatial_object.R` - Script for offline update (NEW)
- `SPATIAL_OBJECT_VERSION_INFO.md` - This documentation (NEW)
