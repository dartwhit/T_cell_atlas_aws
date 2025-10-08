# SSc skin cell Atlas

This project is a web-based application for exploring single-cell and spatial transcriptomics data related to Systemic Sclerosis (SSc). It's a Shiny app written in R, designed to be hosted on AWS.

## 🚀 Live Application

The application is hosted on AWS and can be accessed [here](http://3.16.156.124/atlas/):

## 🚀 Features

*   **Data Exploration:** Explore and visualize different transcriptomics datasets from published studies.
*   **Interactive Visualizations:**
    *   **UMAP Plots:** Visualize cell clusters and relationships.
    *   **Feature Plots:** Visualize gene expression and pathway activity on UMAP plots.
    *   **Violin Plots, Dot Plots, and Heatmaps:** Compare gene expression across conditions and cell types.
*   **Differential Expression Analysis:** View and download tables of differentially expressed genes (DEGs).
*   **Pathway Analysis:** Analyze gene sets and pathways using the VAM method.
*   **Spatial Data Exploration:** Visualize gene expression in the context of tissue architecture.

## 📊 Dynamic Dataset Configuration

The application uses a **dynamic dataset configuration system** that makes dataset management simple and maintainable:

### Key Benefits
- **Single Source of Truth**: All dataset configuration in `app_code/config/datasets.tsv`
- **Dynamic Study Levels**: Each dataset defines its own available study levels (full, fibroblasts, immune cells, myeloid cells, etc.)
- **Easy Maintenance**: Adding new datasets requires only TSV updates - no code changes
- **Intelligent File Detection**: Automatically finds data files using naming patterns
- **Robust Testing**: Comprehensive GitHub Actions workflows ensure reliability

### Adding New Datasets
To add a new dataset, simply add a row to `app_code/config/datasets.tsv`:

```tsv
new_study	New Study 2024	RNA	50000	imgs/new.png	Description	TRUE	FALSE	full:Full;immune:Immune cells	NewStudy	NA	NA	2024-01-01	NA
```

### Configuration Format
The enhanced `datasets.tsv` includes:
- **Basic Info**: `id`, `name`, `assay`, `n_cells`, `image`, `desc`
- **Capabilities**: `has_scrna`, `has_spatial`
- **Study Levels**: `study_levels` (e.g., `full:Full;fib:Fibroblasts;immune:Immune cells`)
- **File Patterns**: `file_base_pattern`, `date_pattern`, `alt_patterns`
- **Special Files**: `spatial_file`, `meta_file`

## 🛠️ Technologies Used

*   **Backend:** R, Shiny
*   **Bioinformatics:** Seurat, VAM
*   **Frontend:** HTML, CSS, JavaScript
*   **Deployment:** Docker, AWS, GitHub Actions

## 📊 Data

The application uses several publicly available datasets from studies on Systemic Sclerosis and related conditions. The data is in the form of Seurat objects and includes single-cell RNA-seq and spatial transcriptomics data.

The datasets are configured in `app_code/config/datasets.tsv` and include data from the following studies:

*   Tabib et al. 2021
*   Gur et al. 2022
*   Ma et al. 2024
*   Khanna et al. 2022
*   TMKMH integrated dataset

## ☁️ Deployment

This application is designed for deployment on AWS. The `Dockerfile` is used to containerize the application, and the `.github/workflows` directory contains GitHub Actions for automated deployment to a development or production environment.
