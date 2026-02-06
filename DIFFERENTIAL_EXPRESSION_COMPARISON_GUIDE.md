# Differential Expression Comparison Guide

## Overview

The T Cell Atlas application now supports different types of differential expression comparisons for different datasets. This guide explains how the comparison types are configured and displayed in the UI.

## Comparison Types

### 1. Disease-Based Comparison (Default)
Used by most datasets to compare gene expression between disease states.

**Datasets:** Tabib, Gur, Ma, TMKMH  
**Comparison:** SSc (Systemic Sclerosis) vs Healthy  
**Checkbox Label:** "Compare by disease"

### 2. Timepoint-Based Comparison
Used by the Khanna dataset to compare gene expression across treatment timepoints.

**Dataset:** Khanna  
**Comparison:** After Treatment vs Before Treatment  
**Checkbox Label:** "Compare by timepoint"

## Configuration

### Dataset Metadata (dataset_details.json)

Each dataset in `app_code/config/dataset_details.json` now includes:

```json
{
  "dataset_id": {
    "comparison_type": "disease" | "timepoint",
    "comparison_label": "Description of comparison",
    "data_levels": { ... },
    "files": { ... }
  }
}
```

**Example - Disease comparison (Tabib):**
```json
{
  "tabib": {
    "comparison_type": "disease",
    "comparison_label": "SSc vs Healthy",
    ...
  }
}
```

**Example - Timepoint comparison (Khanna):**
```json
{
  "khanna": {
    "comparison_type": "timepoint",
    "comparison_label": "After vs Before Treatment",
    ...
  }
}
```

### Required Files

For comparison mode to work, datasets must have corresponding DE files:

**Full data level:**
- `DE_by_disease_auto` - DE results for auto-annotated clusters
- `DE_by_disease_broad` - DE results for broad clusters

**Other data levels (fib, mye, immune):**
- `DE_by_disease` - DE results for that specific data level

**Note:** For Khanna dataset, these are named `*_DE_by_timepoint.txt` to reflect the actual comparison type, but are accessed via the same `DE_by_disease_*` keys in the configuration for backward compatibility.

## UI Components

### Dynamic Checkbox Label
The checkbox label automatically changes based on the dataset's comparison type:
- **Disease datasets:** "Compare by disease"
- **Khanna dataset:** "Compare by timepoint"

### Comparison Information Display
A colored information box appears below the checkbox showing what comparison is currently displayed:

**When comparison mode is active (checkbox checked):**
- 🔵 Blue box displaying the specific comparison (e.g., "SSc vs Healthy" or "After vs Before Treatment")

**When showing cluster markers (checkbox unchecked):**
- 🟠 Orange box displaying "Cluster markers (all cells)"

**When comparison data is not available:**
- No information box is shown (checkbox may be hidden or disabled)

## Implementation Details

### Files Modified

1. **`app_code/config/dataset_details.json`**
   - Added `comparison_type` and `comparison_label` fields to all datasets
   - Added DE file references for Khanna timepoint comparisons

2. **`app_code/setup.R`**
   - Parse comparison metadata from JSON
   - Store in `dataset_comparison_type` and `dataset_comparison_label` lists

3. **`app_code/ui.R`**
   - Replace static checkbox with dynamic UI output: `uiOutput("by_condition_checkbox")`
   - Add comparison info display: `uiOutput("comparison_info")`

4. **`app_code/server.R`**
   - `output$by_condition_checkbox` - Renders checkbox with appropriate label
   - `output$comparison_info` - Renders information box showing current comparison

### Code Example - Dynamic Checkbox

```r
output$by_condition_checkbox <- renderUI({
  req(sidebar_inputs$study())
  comparison_type <- dataset_comparison_type[[sidebar_inputs$study()]]
  
  label <- if (comparison_type == "timepoint") {
    "Compare by timepoint"
  } else {
    "Compare by disease"
  }
  
  checkboxInput("by_disease", label)
})
```

### Code Example - Comparison Info Display

```r
output$comparison_info <- renderUI({
  req(sidebar_inputs$study())
  
  # Check if comparison data is available
  has_comparison <- # ... check for DE_by_disease files
  
  if (!has_comparison) {
    return(NULL)
  }
  
  comparison_label <- dataset_comparison_label[[sidebar_inputs$study()]]
  
  if (isTruthy(input$by_disease)) {
    # Blue box showing specific comparison
    div(style = "...", comparison_label)
  } else {
    # Orange box showing cluster markers
    div(style = "...", "Cluster markers (all cells)")
  }
})
```

## Adding New Datasets

To add a new dataset with a different comparison type:

1. **Add dataset entry to `dataset_details.json`:**
   ```json
   {
     "new_dataset": {
       "comparison_type": "your_type",
       "comparison_label": "Your Comparison Description",
       "data_levels": { ... },
       "files": {
         "full": {
           "DEGs_auto": "...",
           "DE_by_disease_auto": "...",  // Your comparison files
           ...
         }
       }
     }
   }
   ```

2. **If needed, update the checkbox label logic in `server.R`:**
   ```r
   label <- if (comparison_type == "timepoint") {
     "Compare by timepoint"
   } else if (comparison_type == "your_type") {
     "Compare by your_type"
   } else {
     "Compare by disease"
   }
   ```

3. **Ensure your DE comparison files are in the correct location:**
   - Path: `data/DE_dfs/your_dataset_*_DE_by_*.txt`

## Testing

To test the implementation:

1. Select different datasets and verify:
   - Checkbox label changes appropriately
   - Comparison info box shows correct comparison label
   - Checking/unchecking updates the info box color and text

2. For Khanna dataset specifically:
   - Checkbox should say "Compare by timepoint"
   - When checked, blue box should show "After vs Before Treatment"
   - When unchecked, orange box should show "Cluster markers (all cells)"

3. For disease datasets (Tabib, Gur, Ma, TMKMH):
   - Checkbox should say "Compare by disease"
   - When checked, blue box should show "SSc vs Healthy"
   - When unchecked, orange box should show "Cluster markers (all cells)"

## Future Enhancements

Potential improvements for future versions:

1. **Support for multiple comparison types per dataset**
   - Allow users to select from dropdown if multiple comparisons available

2. **Comparison-specific visualizations**
   - Customize plot labels based on comparison type

3. **Metadata-driven comparison descriptions**
   - Pull comparison details from Seurat object metadata

4. **Comparison validation**
   - Check that referenced DE files exist and show warnings if missing
