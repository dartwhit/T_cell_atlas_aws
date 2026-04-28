# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

An R Shiny web application for interactive exploration of single-cell and spatial transcriptomics data from Systemic Sclerosis (SSc) studies. Live at https://ssccellatlas.dartmouth.edu/. Uses Seurat for data analysis, shinymanager for authentication, and deploys to AWS EC2 via Docker.

## Commands

### Local Development (Debug Mode)
```bash
# Build and run a single container with interactive logging
./docker-run-debug.sh
```

### Production (Multi-container)
```bash
docker-compose build      # Build all containers
docker-compose up -d      # Start 3 Shiny + Nginx (port 80)
docker-compose down       # Stop all
```

### Inspect a Running Container
```bash
docker logs atlas-debug --follow
docker exec -it atlas-debug bash
```

No automated test suite. Testing is done manually via PR deployments (see below).

## Architecture

### Entry Points
```
app_code/
├── app.R        # Entry point — sources setup.R, ui.R, server.R
├── setup.R      # Loads config, reads datasets.tsv and dataset_details.json
├── ui.R         # UI layout (bslib Bootstrap 5)
└── server.R     # All reactive server logic (~960 lines)
```

### Shiny Modules
- `dataset_gallery_module.R` — Dataset discovery/selection (Datasets tab)
- `explore_sidebar_module.R` — Controls for study, data level, gene/pathway selection
- `spatial_unit.R` — Spatial transcriptomics viewer (Visium)

### Configuration
- `config/datasets.tsv` — Dataset metadata (name, assay type, display info)
- `config/dataset_details.json` — Per-dataset file paths, data levels, DE file lists, VAM pathway files

### Data Flow
```
datasets.tsv + dataset_details.json
        ↓ (setup.R on startup)
dataset_choices, dataset_files, data_level_choices lists
        ↓
User selects dataset in Gallery → Sidebar updates
        ↓
"Load Data" button → Seurat RDS file loaded reactively
        ↓
Plots (UMAP, feature, violin/box/dot/heatmap), DEG tables, VAM pathways
```

### Key Reactive Values in server.R
- `seurat_obj` — Currently loaded Seurat object
- `gene_list_obj` — Available genes for the dataset
- `VAM_df` — Pathway analysis results (VAM)
- `DEGs_df` — Differential expression table
- `cell_clusters` — Available cluster identities
- `meta_df` — Cell metadata

### UI Tabs
1. **Datasets** — Gallery with search/filter, triggers dataset selection
2. **Explore** — Main analysis: Plots sub-tab (feature plot, UMAP, expression plot), DEGs table, Metadata
3. **Spatial data explorer** — Visium spatial viewer

### Authentication
Shinymanager wraps the entire UI. Credentials stored in SQLite (`users_current.sqlite`), mounted as a shared Docker volume (`users-db`) across all 3 containers.

## Deployment

### PR Deployments (Automatic)
Every PR triggers `.github/workflows/deploy-pr.yml`:
- Deploys to EC2, port `3800 + (PR_number % 10)` (e.g., PR #3 → port 3803)
- Container named `atlas-pr-{PR_NUMBER}`
- GitHub Actions posts the URL as a PR comment
- Note: PRs with the same `PR_number % 10` will collide (e.g., PR #5 and PR #15)

### Production Deployment
Push to `main` triggers `.github/workflows/deploy.yml` — deploys 3-container stack with Nginx load balancer.

### Data Volumes
Data files (Seurat RDS, gene lists, DEG tables, VAM files) are mounted read-only at `/srv/shiny-server/atlas/data/`. They are NOT in the repo. File paths are configured in `config/dataset_details.json`.

## Key Dependencies

| Package | Purpose |
|---------|---------|
| Seurat | scRNA-seq and spatial data analysis |
| VAM | Pathway/gene set scoring |
| shinymanager | Authentication |
| bslib | Bootstrap 5 theming |
| DT | Interactive DEG tables |
| plotly | Interactive plots |
| igraph | Network analysis |
