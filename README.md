# SSc skin cell Atlas

[![Deploy to EC2](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/deploy.yml/badge.svg)](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/deploy.yml)
[![Deploy PR](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/deploy-pr.yml/badge.svg)](https://github.com/dartwhit/T_cell_atlas_aws/actions/workflows/deploy-pr.yml)

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

### Pull Request Deployment

All pull requests to the `main` branch are **automatically deployed to EC2** for manual testing. This allows you to test the full UI and backend functionality in a real environment before merging.

**How It Works:**
- When you open or update a PR, GitHub Actions automatically deploys your branch to EC2
- Each PR gets a unique port (3800 + PR number % 100)
- A comment is posted on the PR with the deployment URL
- Deployments typically complete in 2-3 minutes

**Manual Testing:**
1. Open or update your PR
2. Wait for the deployment workflow to complete
3. Click the URL in the automated PR comment
4. Test your changes in the live environment
5. Check UI functionality and backend operations

**Access Information:**
- The deployment URL will be posted as a comment on your PR
- Example: `http://your-ec2-host:38XX/atlas/`
- Each PR maintains its own isolated container and cache

**Container Management:**
- PR containers are automatically cleaned up (keeping only the 5 most recent)
- Containers use isolated cache and log directories
- View logs on EC2: `docker logs atlas-pr-{PR_NUMBER}`

**Note:** This replaces the previous automated testing workflow, which was slow (~1 hour). Manual testing on EC2 provides better validation of real-world functionality.
