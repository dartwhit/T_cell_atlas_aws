# Dynamic Dataset Configuration System

## Overview

The T-Cell Atlas AWS application now uses a **dynamic dataset configuration system** that makes `datasets.tsv` the single source of truth for all dataset information. This eliminates the need for hardcoded configurations and makes adding new datasets as simple as adding a row to a TSV file.

## 🎯 Key Benefits

- **Single Source of Truth**: All dataset configuration in one TSV file
- **Dynamic Study Levels**: Each dataset can define its own available study levels (full, fib, immune, mye, etc.)
- **Intelligent File Detection**: Automatically finds data files using naming patterns and fallback strategies
- **Easy Maintenance**: Adding new datasets requires only TSV updates - no code changes
- **Backward Compatible**: Existing modules continue to work without modification
- **Robust Error Handling**: Graceful handling of missing files with detailed warnings

## 📁 File Structure

```
app_code/
├── config/
│   └── datasets.tsv          # ← Single source of truth
├── setup.R                   # ← Dynamic configuration loader
├── data/                     # ← Data files (auto-detected)
├── modules/                  # ← App modules (unchanged)
└── tests/                    # ← Test suite
```

## 📋 Enhanced datasets.tsv Format

### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `id` | Unique dataset identifier | `tmkmh` |
| `name` | Human-readable name | `TMKMH integrated dataset` |
| `assay` | Assay type | `RNA` |
| `n_cells` | Number of cells | `173299` |
| `image` | Path to dataset image | `imgs/TMKMH_img.png` |
| `desc` | Dataset description | `Integrated dataset...` |
| `has_scrna` | Has single-cell RNA data | `TRUE` |
| `has_spatial` | Has spatial data | `FALSE` |

### New Dynamic Configuration Columns

| Column | Description | Example |
|--------|-------------|---------|
| `study_levels` | Available study levels with labels | `full:Full;fib:Fibroblasts;immune:Immune cells` |
| `file_base_pattern` | Base pattern for file names | `TMKMH` |
| `spatial_file` | Spatial data file (if applicable) | `spatial_data.rds` |
| `meta_file` | Metadata file (if applicable) | `TMKMH_metadata.txt` |
| `date_pattern` | Date pattern in filenames | `2025-06-06` |
| `alt_patterns` | Alternative file patterns | `TMKMH_filtered_reduced_*_VAM.rds` |

### Study Levels Format

The `study_levels` column uses semicolon-separated `code:label` pairs:

```
full:Full;fib:Fibroblasts;immune:Immune cells;mye:Myeloid cells
```

This creates:
- `full` level with label "Full"
- `fib` level with label "Fibroblasts"
- `immune` level with label "Immune cells"
- `mye` level with label "Myeloid cells"

## 🔧 How It Works

### 1. Dynamic File Detection

The system uses intelligent pattern matching to find files:

**For Full Datasets:**
1. `{pattern}_filtered_reduced_{date}_VAM.rds`
2. `{pattern}_full_filtered_VAM.rds`
3. `{pattern}_filtered_VAM.rds`
4. Alternative patterns from `alt_patterns` column

**For Subsets:**
1. `{pattern}_{level}_filtered_VAM.rds`
2. `{pattern}_{level}_gene_list.rds`
3. `{pattern}_{level}_DEGs_sig.txt`
4. etc.

### 2. Automatic Structure Generation

The `setup.R` script automatically generates:

- **`dataset_files`**: Complete file path mappings
- **`data_level_choices`**: UI choice lists for each dataset
- **`datasets`**: Metadata object for modules

### 3. Validation and Error Handling

- Validates TSV structure and required columns
- Warns about missing data files
- Provides fallback file patterns
- Graceful degradation for incomplete datasets

## 🚀 Adding New Datasets

To add a new dataset, simply add a row to `datasets.tsv`:

```tsv
new_study	New Study 2024	RNA	50000	imgs/new_study.png	A new exciting study	TRUE	FALSE	full:Full;immune:Immune cells	NewStudy	NA	NA	2024-12-01	NewStudy_*_VAM.rds
```

That's it! No code changes needed.

## 🧪 Testing

### Local Testing

```bash
# Basic configuration test
Rscript test_dynamic_setup.R

# Comprehensive test suite
Rscript tests/test_dataset_config.R

# EC2 deployment simulation
Rscript test_ec2_deployment.R
```

### Automated Testing

GitHub Actions automatically test:
- Configuration loading
- Data structure validation
- File existence checks
- Backward compatibility
- Module integration

## 📊 Test Results

The test suite validates:

✅ **Configuration Loading** - TSV parsing and setup script execution  
✅ **Data Structures** - Proper generation of `dataset_files` and `data_level_choices`  
✅ **File Detection** - Intelligent file pattern matching  
✅ **Study Levels** - Dynamic level configuration per dataset  
✅ **Backward Compatibility** - Existing modules continue to work  
✅ **Error Handling** - Graceful handling of missing files  

## 🔄 Migration from Hardcoded System

### Before (Hardcoded)
```r
# In setup.R - 200+ lines of hardcoded file paths
dataset_files <- list(
  "tmkmh" = list(
    "full" = list("seurat" = "TMKMH_filtered_reduced_2025-06-06_VAM.rds"),
    "fib" = list("seurat" = "TMKMH_fib_filtered_VAM.rds"),
    # ... many more hardcoded paths
  )
)

data_level_choices <- list(
  tmkmh = c("Full" = "full", "Fibroblasts" = "fib", "Immune cells" = "immune")
  # ... hardcoded for each dataset
)
```

### After (Dynamic)
```r
# In setup.R - Automatically generated from TSV
source("setup.R")  # That's it!
# dataset_files and data_level_choices are generated automatically
```

### TSV Configuration
```tsv
id    study_levels                              file_base_pattern
tmkmh full:Full;fib:Fibroblasts;immune:Immune cells  TMKMH
```

## 🛠️ Technical Details

### File Pattern Detection Algorithm

1. **Exact Match**: Check if file exists as specified
2. **Wildcard Expansion**: Replace `*` with regex patterns
3. **Date Substitution**: Use `date_pattern` for dated files
4. **Fallback Patterns**: Try alternative naming conventions
5. **Graceful Degradation**: Warn about missing files but continue

### Performance Optimizations

- **Lazy Loading**: Files are detected at startup, not during runtime
- **Caching**: File paths are cached after first detection
- **Batch Validation**: All files validated in single pass
- **Minimal I/O**: Only checks file existence, doesn't load content

### Error Recovery

- **Missing Files**: Warns but continues with available files
- **Invalid Patterns**: Falls back to standard naming conventions
- **Corrupted TSV**: Provides helpful error messages
- **Partial Data**: Application remains functional with available datasets

## 🔍 Debugging

### Common Issues and Solutions

1. **"Missing seurat paths" Error**
   - Check `file_base_pattern` matches actual file names
   - Verify `date_pattern` if using dated files
   - Add `alt_patterns` for non-standard naming

2. **"Study levels parsing failed"**
   - Ensure format: `code:label;code:label`
   - Check for special characters in labels
   - Validate no extra semicolons

3. **"Data directory not found"**
   - Verify `inDir` variable points to correct location
   - Check file permissions on data directory
   - Ensure relative paths are correct

### Debugging Tools

```r
# In R console
source("app_code/setup.R")

# Check parsed levels
str(dataset_meta$parsed_levels)

# Check detected files
str(dataset_files, max.level = 2)

# Check missing files
# (automatically printed during setup)
```

## 🎯 Future Enhancements

1. **File Validation**: Add MD5 checksums for data integrity
2. **Auto-Discovery**: Automatically detect new datasets by scanning file names
3. **Configuration UI**: Web interface for managing dataset configuration
4. **Version Control**: Track changes to dataset configurations
5. **Cloud Integration**: Support for S3 and other cloud storage backends

## 📖 Example Configuration

```tsv
id	name	assay	n_cells	image	desc	has_scrna	has_spatial	study_levels	file_base_pattern	spatial_file	meta_file	date_pattern	alt_patterns
tmkmh	TMKMH integrated dataset	RNA	173299	imgs/TMKMH_img.png	Integrated dataset from TMKMH study.	TRUE	FALSE	full:Full;fib:Fibroblasts;immune:Immune cells	TMKMH	NA	NA	2025-06-06	TMKMH_filtered_reduced_*_VAM.rds
ma	Ma et al. 2024	RNA	73926	imgs/Ma_img.png	Ma dataset with spatial data.	TRUE	TRUE	full:Full;mye:Myeloid cells;fib:Fibroblasts	Ma	spatial_data.rds	NA	2025-01-28	Ma_filtered_reduced_*_VAM.rds
```

This configuration automatically generates all necessary data structures for the application to run with these datasets, including appropriate study level choices and file path mappings.

---

🎉 **The dynamic configuration system makes dataset management effortless while maintaining full backward compatibility with existing code!**