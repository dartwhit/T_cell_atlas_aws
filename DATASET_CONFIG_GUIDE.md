# Dataset Configuration Guide

## Overview

This guide explains the improved dataset configuration system that makes `datasets.tsv` the single source of truth for all dataset information, eliminating the need for hardcoded configurations in `setup.R`.

## Key Improvements

### 1. **Single Source of Truth**

-   All dataset metadata is now in `datasets.tsv`
-   Study levels are dynamically generated from the TSV file
-   File paths are intelligently detected using naming patterns
-   No more hardcoded `dataset_files` or `data_level_choices`

### 2. **Flexible Study Levels**

Instead of hardcoding study levels, they are now defined in the TSV using the format:

```         
study_levels: full:Full;fib:Fibroblasts;immune:Immune cells;mye:Myeloid cells
```

This allows each dataset to have its own set of available study levels.

### 3. **Intelligent File Detection**

The system automatically finds files using naming patterns and fallback strategies: - Tries multiple naming conventions - Handles date patterns in filenames - Validates file existence and warns about missing files - Uses alternative patterns if primary patterns don't match

## Enhanced datasets.tsv Format

### Required Columns

| Column | Description | Example |
|--------------------|------------------------------|----------------------|
| `id` | Unique dataset identifier | `tmkmh` |
| `name` | Human-readable name | `TMKMH integrated dataset` |
| `assay` | Assay type | `RNA` |
| `n_cells` | Number of cells | `173299` |
| `image` | Path to dataset image | `imgs/TMKMH_img.png` |
| `desc` | Dataset description | `Integrated dataset from TMKMH study.` |
| `has_scrna` | Has single-cell RNA data | `TRUE` |
| `has_spatial` | Has spatial data | `FALSE` |

### New Columns for Dynamic Configuration

| Column | Description | Example |
|--------------------|------------------------------|----------------------|
| `study_levels` | Available study levels with labels | `full:Full;fib:Fibroblasts;immune:Immune cells` |
| `file_base_pattern` | Base pattern for file names | `TMKMH` |
| `spatial_file` | Spatial data file (if applicable) | `spatial_data.rds` |
| `meta_file` | Metadata file (if applicable) | `TMKMH_metadata.txt` |
| `date_pattern` | Date pattern in filenames | `2025-06-06` |
| `alt_patterns` | Alternative file patterns | `TMKMH_filtered_reduced_*_VAM.rds` |

### Study Levels Format

The `study_levels` column uses a semicolon-separated list of `code:label` pairs:

```         
full:Full;fib:Fibroblasts;immune:Immune cells;mye:Myeloid cells
```

This creates: - `full` level with label "Full" - `fib` level with label "Fibroblasts"\
- `immune` level with label "Immune cells" - `mye` level with label "Myeloid cells"

## File Naming Patterns

The system uses intelligent pattern matching to find files:

### For Full Datasets:

1.  `{pattern}_filtered_reduced_{date}_VAM.rds`
2.  `{pattern}_full_filtered_VAM.rds`
3.  `{pattern}_filtered_VAM.rds`
4.  Alternative patterns from `alt_patterns` column

### For Subsets:

1.  `{pattern}_{level}_filtered_VAM.rds`
2.  `{pattern}_{level}_gene_list.rds`
3.  `{pattern}_{level}_DEGs_sig.txt`
4.  etc.

## Benefits

### 1. **Easy Dataset Updates**

To add a new dataset, simply add a row to `datasets.tsv`. No code changes needed!

### 2. **Flexible Study Levels**

Each dataset can have different combinations of study levels (full, fib, immune, mye, etc.).

### 3. **Automatic File Discovery**

The system finds files automatically using naming conventions and handles variations.

### 4. **Validation & Error Handling**

-   Warns about missing files
-   Validates file existence
-   Graceful fallbacks for different naming patterns

### 5. **Maintainability**

-   No more massive hardcoded structures in `setup.R`
-   Clear separation of data and code
-   Easy to understand and modify

## Migration Guide

### Step 1: Update datasets.tsv

Replace the current `datasets.tsv` with the enhanced format that includes the new columns.

### Step 2: Replace setup.R

Replace the current `setup.R` with the new dynamic version that reads from the TSV.

### Step 3: Test

Verify that all datasets load correctly and study levels appear as expected.

## Example Configuration

``` tsv
id  name    assay   n_cells image   desc    has_scrna   has_spatial study_levels    file_base_pattern   spatial_file    meta_file   date_pattern    alt_patterns
tmkmh   TMKMH integrated dataset    RNA 173299  imgs/TMKMH_img.png  Integrated dataset from TMKMH study.    TRUE    FALSE   full:Full;fib:Fibroblasts;immune:Immune cells   TMKMH   NA  NA  2025-06-06  TMKMH_filtered_reduced_*_VAM.rds
ma  Ma et al. 2024  RNA 73926   imgs/Ma_img.png Ma dataset. TRUE    TRUE    full:Full;mye:Myeloid cells;fib:Fibroblasts Ma  spatial_data.rds    NA  2025-01-28  Ma_filtered_reduced_*_VAM.rds
```

This configuration: - Defines TMKMH with full, fibroblast, and immune cell levels - Defines Ma with full, myeloid, and fibroblast levels\
- Includes spatial data for Ma dataset - Uses intelligent file pattern matching

## Future Enhancements

1.  **File Path Validation**: Add checks for file existence during app startup
2.  **Custom File Mappings**: Allow override of specific file paths if needed
3.  **Metadata Validation**: Validate TSV structure and required fields
4.  **Auto-detection**: Automatically detect available study levels by scanning file names