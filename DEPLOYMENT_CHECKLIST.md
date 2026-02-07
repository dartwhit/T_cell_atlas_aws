# Deployment Checklist for DE Comparison UI Changes

## Prerequisites

Before deploying these changes, ensure the following data files exist in the deployment environment:

### Khanna Timepoint Differential Expression Files

These files need to be generated and placed in `data/DE_dfs/` directory:

1. **Full dataset - Auto annotation:**
   - File: `Khanna_auto_DE_by_timepoint.txt`
   - Content: DE results comparing "After Treatment" vs "Before Treatment" for auto-annotated clusters
   - Format: Tab-separated values with columns: cluster, gene, avg_log2FC, p_val, p_val_adj, etc.

2. **Full dataset - Broad annotation:**
   - File: `Khanna_broad_DE_by_timepoint.txt`
   - Content: DE results comparing "After Treatment" vs "Before Treatment" for broad clusters
   - Format: Same as above

3. **Fibroblasts subset:**
   - File: `Khanna_fib_DE_by_timepoint.txt`
   - Content: DE results for fibroblast cells only
   - Format: Same as above

4. **Myeloid cells subset:**
   - File: `Khanna_mye_DE_by_timepoint.txt`
   - Content: DE results for myeloid cells only
   - Format: Same as above

## Testing Procedure

### 1. Visual Testing - Khanna Dataset

**a. Initial Load**
- [ ] Select "Khanna et al. 2022" from dataset dropdown
- [ ] Navigate to "DEGs table" tab
- [ ] Click "Load Data" button
- [ ] Verify checkbox shows "Compare by timepoint" (not "Compare by disease")
- [ ] Verify orange info box appears showing "Cluster markers (all cells)"

**b. Timepoint Comparison Mode**
- [ ] Check the "Compare by timepoint" checkbox
- [ ] Verify blue info box appears showing "After vs Before Treatment"
- [ ] Verify DEG table updates to show timepoint comparison results
- [ ] Select different cell clusters and verify data updates

**c. Toggle Back**
- [ ] Uncheck the "Compare by timepoint" checkbox
- [ ] Verify info box changes back to orange with "Cluster markers (all cells)"
- [ ] Verify DEG table updates to show cluster markers

### 2. Visual Testing - Disease Datasets (Tabib, Gur, Ma, TMKMH)

**a. Initial Load**
- [ ] Select one of: "Tabib et al. 2021", "Gur et al. 2022", "Ma et al. 2024", or "TMKMH"
- [ ] Navigate to "DEGs table" tab
- [ ] Click "Load Data" button
- [ ] Verify checkbox shows "Compare by disease"
- [ ] Verify orange info box appears showing "Cluster markers (all cells)"

**b. Disease Comparison Mode**
- [ ] Check the "Compare by disease" checkbox
- [ ] Verify blue info box appears showing "SSc vs Healthy"
- [ ] Verify DEG table updates to show disease comparison results
- [ ] Select different cell clusters and verify data updates

**c. Toggle Back**
- [ ] Uncheck the "Compare by disease" checkbox
- [ ] Verify info box changes back to orange with "Cluster markers (all cells)"
- [ ] Verify DEG table updates to show cluster markers

### 3. Data Level Testing

For each dataset, test all available data levels:

**Khanna:**
- [ ] Full
- [ ] Myeloid cells
- [ ] Fibroblasts

**Other datasets:**
- [ ] Full
- [ ] Fibroblasts
- [ ] Immune cells (where available)
- [ ] Myeloid cells (where available)

Verify that:
- [ ] Info box appears/disappears appropriately based on data availability
- [ ] Checkbox state persists when switching data levels
- [ ] Correct comparison label is shown

### 4. Edge Cases

- [ ] Switch between datasets while comparison mode is active
- [ ] Switch data levels while comparison mode is active
- [ ] Change annotation type (Original seurat clusters / OFF) while comparison mode is active
- [ ] Rapidly toggle checkbox on/off
- [ ] Test with slow network conditions

### 5. Cross-Browser Testing

Test on:
- [ ] Chrome/Edge (Chromium-based)
- [ ] Firefox
- [ ] Safari
- [ ] Mobile browsers (optional)

## Known Limitations

1. **File Naming Convention:** While Khanna files are named `*_DE_by_timepoint.txt`, they are accessed via the same `DE_by_disease_*` keys in the configuration for backward compatibility with existing code structure.

2. **Checkbox Input Name:** The checkbox input remains `by_disease` internally, even though the label changes to "Compare by timepoint" for Khanna. This is intentional to avoid breaking existing code dependencies.

3. **Missing Data:** If timepoint DE files are not present, the checkbox and info box will not appear for Khanna dataset when in comparison mode.

## Rollback Plan

If issues are discovered after deployment:

1. **Quick Fix:** Temporarily remove `comparison_type` and `comparison_label` from `dataset_details.json` for Khanna. The app will fall back to showing "Compare by disease".

2. **Full Rollback:** Revert to commit `7d9c997` (before these changes).

## Post-Deployment Verification

- [ ] Check application logs for any errors related to missing files
- [ ] Verify all datasets load correctly
- [ ] Spot-check a few DEG visualizations to ensure data is correct
- [ ] Monitor user feedback for any confusion about the new UI elements

## Data Generation Script (If Needed)

If the Khanna timepoint DE files don't exist yet, they need to be generated using the same pipeline as the other datasets' `DE_by_disease` files, but with timepoint as the comparison variable instead of disease status.

Example R pseudocode:
```r
# Load Khanna Seurat object
khanna <- readRDS("Khanna_filtered_reduced_2025-01-28_VAM.rds")

# Perform differential expression by timepoint
# (Assuming timepoint metadata is stored in a column like "timepoint" or "treatment_status")
de_results <- FindMarkers(
  khanna,
  ident.1 = "After",      # After treatment
  ident.2 = "Before",     # Before treatment
  group.by = "timepoint", # or whatever the metadata column is named
  test.use = "wilcox",    # or appropriate test
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Filter and save results
write.table(de_results, "Khanna_auto_DE_by_timepoint.txt", sep = "\t", quote = FALSE)
```

Repeat for:
- Broad annotations
- Fibroblast subset
- Myeloid subset

## Support Contacts

If issues arise:
- Technical issues: [Repository maintainer]
- Data generation: [Bioinformatics team]
- UI/UX feedback: [Product team]
