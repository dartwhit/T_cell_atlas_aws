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

### Run Locally Without Docker (Auth Disabled)
```bash
cd app_code
LOCAL_DEV=1 Rscript -e "shiny::runApp('.')"   # 1/true/yes/on skips shinymanager login
```

### Tests
No app-wide test suite. The CellChat explorer has a standalone R test harness in `tests/cellchat/` (`test_wiring.R`, `test_plots.R`, `helpers.R`) — run per its `tests/cellchat/README.md`. Broader testing is done manually via PR deployments (see below).

## Architecture

### Entry Points
```
app_code/
├── app.R        # Entry point — sources setup.R, ui.R, server.R
├── setup.R      # Loads config, reads datasets.tsv and dataset_details.json; defines read_object()
├── ui.R         # UI layout (bslib Bootstrap 5)
├── server.R     # All reactive server logic (~1200 lines)
├── config/      # datasets.tsv + dataset_details.json
├── modules/     # Shiny modules (see below)
├── www/         # Static assets: styles.css, dataset thumbnails (imgs/)
└── data/        # Seurat/CellChat data (gitignored, not in repo)
```

### Shiny Modules (`app_code/modules/`)
- `dataset_gallery_module.R` — Dataset discovery/selection (Datasets tab)
- `explore_sidebar_module.R` — Controls for study, data level, gene/pathway selection
- `spatial_unit.R` — Spatial transcriptomics viewer (Visium)
- `cellchat_explorer_module.R` — CellChat cell-cell communication explorer (CellChat tab UI + server)
- `cellchat_helpers.R` — Load/validate stripped CellChat objects; network/pathway/heatmap plot helpers

### Data Loading (`read_object()` in setup.R)
Objects are loaded via `read_object(path)`, which prefers a fast `.qs2` sibling
(`qs2::qs_read`, ~7x faster than gzipped RDS) when present, else falls back to
`readRDS()`. So the app works before/without conversion. Convert data with the
scripts in `scripts/` (`convert_rds_to_qs2.R`, `slim_qs2_inplace.R`). Thread count
for reads is controlled by `QS2_READ_THREADS`.

### Configuration
- `config/datasets.tsv` — Dataset metadata (id, name, assay, n_cells, file path, image, tags, has_scrna/has_spatial)
- `config/dataset_details.json` — Per-dataset file paths, data levels, DE file lists, VAM pathway files, and optional `cellchat` block

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
4. **CellChat** — Cell-cell communication explorer; enabled only for studies with a complete `cellchat` block in `dataset_details.json` (currently `tabib` and `tmkmh`)

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
Data files (Seurat RDS/qs2, gene lists, DEG tables, VAM files, CellChat objects) are mounted read-only at `/srv/shiny-server/atlas/data/`. They are NOT in the repo. File paths are configured in `config/dataset_details.json`. Stripped CellChat objects live under `<study>/cellchat/` (`merged.rds`, `HC.rds`, `SSc.rds`).

## Key Dependencies

| Package | Purpose |
|---------|---------|
| Seurat | scRNA-seq and spatial data analysis |
| VAM | Pathway/gene set scoring |
| CellChat | Cell-cell communication (CellChat tab; optional — tab shows an error if the package or data is missing) |
| qs2 | Fast object serialization (optional; falls back to readRDS) |
| shinymanager | Authentication |
| bslib / bsicons | Bootstrap 5 theming and icons |
| shinyWidgets / shinyjs / shinycssloaders | UI controls, JS helpers, loading spinners |
| periscope2 | Downloadable plot/table widgets |
| DT | Interactive DEG tables |
| dplyr / tidyr / stringr | Data wrangling |
| ggplot2 | Plotting |
| jsonlite | Reads dataset_details.json |
