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

### Running Locally (Auth Disabled)

For quick local development you can run the Shiny app directly in R — outside Docker — with authentication disabled. Set the `LOCAL_DEV` environment variable to an explicit truthy value (`1`, `true`, `yes`, or `on`) to skip the `shinymanager` login.

**Requirements:**
- R with the app's packages installed (Seurat, VAM, shinymanager, bslib, DT, plotly, shinyjs, etc.)
- A local copy of the data folder at `app_code/data/`, with one subfolder per study matching the dataset IDs in `app_code/config/datasets.tsv` (e.g. `tmkmh/`, `tabib/`, `gur/`, `ma/`, `khanna/`, `li/`)

**Run from a terminal:**
```bash
cd app_code
LOCAL_DEV=1 Rscript -e "shiny::runApp('.')"
```

**Or from the R console / RStudio "Run App":**
```r
Sys.setenv(LOCAL_DEV = "1")
setwd("app_code")
shiny::runApp(".")
```

The app prints its working directory to stderr on startup, so you can confirm it resolves data paths against `app_code/data/`.

> **Note:** `LOCAL_DEV` only disables authentication for local development. When unset (the default in Docker / production), the UI is wrapped with `shinymanager` and the SQLite credential database is required as normal.

### CellChat data

The CellChat tab is enabled only for studies with a complete CellChat entry in
[app_code/config/dataset_details.json](app_code/config/dataset_details.json).
The initial Tabib integration expects these precomputed, stripped CellChat
objects in the study data directory:

```
app_code/data/tabib/cellchat/merged.rds
app_code/data/tabib/cellchat/HC.rds
app_code/data/tabib/cellchat/SSc.rds
```

For deployment, copy the same directory structure to the mounted EC2 data
volume under `tabib/cellchat/`. These objects are data assets and must not be
committed to the repository. The CellChat tab displays a load error if one is
missing or invalid.

### Pull Request Deployment

All pull requests to the `main` branch are **automatically deployed to EC2** for manual testing. This allows you to test the full UI and backend functionality in a real environment before merging.

**How It Works:**
- When you open or update a PR, GitHub Actions automatically deploys your branch to EC2
- Each PR gets a unique port (3800-3809, based on PR number % 10)
  - **Note:** PRs with numbers differing by 10 (e.g., PR #5 and #15) will use the same port. The newer deployment will replace the older one.
- **A comment is automatically posted on your PR** (in the Conversation tab) with the deployment URL
- Deployments typically complete in 2-3 minutes

**Finding Your Deployment:**
- After the deployment workflow completes, check the **Conversation** tab on your PR
- Look for an automated comment titled "🚀 PR Deployed for Manual Testing"
- The comment includes:
  - Direct access URL to your deployment
  - Container name for debugging
  - Testing instructions
- If you don't see the comment, check the Actions tab for any errors in the "Comment on PR with deployment info" step

**EC2 Security Group Setup:**
To access PR deployments, add an inbound rule to your EC2 security group:
- **Type:** Custom TCP
- **Port Range:** 3800-3809
- **Source:** Your IP address or `0.0.0.0/0` (for public access)
- **Description:** PR testing deployments

**Manual Testing:**
1. Open or update your PR
2. Wait for the deployment workflow to complete (~2-3 minutes)
3. Go to your PR's **Conversation** tab
4. Find the automated comment "🚀 PR Deployed for Manual Testing"
5. Click the deployment URL in the comment
6. Test your changes in the live environment
7. Verify UI functionality and backend operations

**Access Information:**
- The deployment URL will be posted as an automated comment on your PR's Conversation tab
- Example URL format: `http://your-ec2-host:38XX/atlas/`
- Each PR maintains its own isolated container and cache
- The comment is updated automatically when you push new commits

**Container Management:**
- PR containers are automatically cleaned up (keeping only the 5 most recent)
- Containers use isolated cache and log directories
- View logs on EC2: `docker logs atlas-pr-{PR_NUMBER}`

**Note:** This replaces the previous automated testing workflow, which was slow (~1 hour). Manual testing on EC2 provides better validation of real-world functionality.
