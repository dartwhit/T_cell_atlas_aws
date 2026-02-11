# SSc skin cell Atlas

[![Deploy to EC2](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/deploy.yml/badge.svg)](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/deploy.yml)
[![Test PR](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/test-pr.yml/badge.svg)](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/test-pr.yml)

This project is a web-based application for exploring single-cell and spatial transcriptomics data related to Systemic Sclerosis (SSc). It's a Shiny app written in R, designed to be hosted on AWS.

## 🚀 Live Application

The application is hosted on AWS and can be accessed [here](https://ssccellatlas.dartmouth.edu/):


## 🚀 Features

*   **Data Exploration:** Explore and visualize different transcriptomics datasets from published studies.
*   **Interactive Visualizations:**
    *   **UMAP Plots:** Visualize cell clusters and relationships.
    *   **Feature Plots:** Visualize gene expression and pathway activity on UMAP plots.
    *   **Violin Plots, Dot Plots, and Heatmaps:** Compare gene expression across conditions and cell types.
*   **Differential Expression Analysis:** View and download tables of differentially expressed genes (DEGs).
*   **Pathway Analysis:** Analyze gene sets and pathways using the VAM method.
*   **Spatial Data Exploration:** Visualize gene expression in the context of tissue architecture.

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

## 🧪 Development & Testing

### Pull Request Testing

All pull requests to the `main` branch are automatically tested using GitHub Actions. The PR testing workflow (`.github/workflows/test-pr.yml`) performs the following checks:

**Code Quality Checks:**
- ✅ **R Syntax Validation:** Ensures all R files can be parsed without syntax errors
- ℹ️ **R Linting:** Runs `lintr` to identify code style issues (informational only)

**Functional Tests:**
- ✅ **Shiny App Loading:** Verifies the app can initialize without immediate errors
- ✅ **Docker Build:** Ensures the Docker image builds successfully
- ✅ **Container Startup:** Tests that the containerized app can start

**Requirements:**
- Tests run on `ubuntu-latest` with R 4.3
- All dependencies from the Dockerfile are installed
- Tests are designed to be practical and non-blocking for data-related issues

The workflow provides clear feedback on test results and will block merging if critical issues are found (syntax errors, build failures).
