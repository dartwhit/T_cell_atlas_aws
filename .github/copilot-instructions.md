# Copilot Instructions — SSc Cell Atlas Shiny App

## What This Is

R Shiny app for exploration of single-cell and spatial transcriptomics data from SSc studies.
Live at https://ssccellatlas.dartmouth.edu/. Deployed to AWS EC2 via Docker.

## Architecture

```
app_code/
├── app.R          # Entry point: sources setup.R, ui.R, server.R
├── setup.R        # Reads datasets.tsv + dataset_details.json → dataset_files, data_level_choices
├── ui.R           # bslib Bootstrap 5 layout
└── server.R       # All reactive logic (~960 lines)
modules/
├── dataset_gallery_module.R
├── explore_sidebar_module.R
└── spatial_unit.R             # Visium spatial viewer
config/
├── datasets.tsv               # One row per dataset (id, name, has_scrna, has_spatial, …)
└── dataset_details.json       # Per-dataset: data_levels map, file names per level, comparison type
```

Key reactive values in `server.R`: `seurat_obj`, `gene_list_obj`, `VAM_df`, `DEGs_df`, `cell_clusters`, `meta_df`.

## EC2 Data Volume — Critical

**Dataset files are never in the repo.** They live on the EC2 host at `/home/ubuntu/atlas/data/` and are bind-mounted read-only into containers:

```
/home/ubuntu/atlas/data/   →   /srv/shiny-server/atlas/data/   (inside container, read-only)
```

File names in this volume must exactly match the strings in `config/dataset_details.json`.
The `inDir` variable in `server.R` is `"data/"`, so the resolved path inside the container is `data/<filename>`.

### Adding a New Dataset

1. Place all RDS/txt files on the EC2 host under `/home/ubuntu/atlas/data/`.
2. Add a row to `config/datasets.tsv`.
3. Add an entry to `config/dataset_details.json` with matching file names.
4. PR → auto-deploy → verify in the PR container.

## Build & Deploy

### Local debug (single container, interactive logs)
```bash
./docker-run-debug.sh
```

### PR deployment (automatic via GitHub Actions)
Every PR to `main` triggers `.github/workflows/deploy-pr.yml`:
- Container: `atlas-pr-{PR_NUMBER}`, port: `3800 + (PR_NUMBER % 10)`
- Data volume: `/home/ubuntu/atlas/data` (same production data, read-only)
- GitHub Actions posts the access URL as a PR comment

### Production (push to `main`)
```bash
docker-compose build
docker-compose up -d   # 3 Shiny containers + Nginx on port 80
```

## Debugging on EC2

SSH access is via the `EC2_HOST` / `EC2_USER` / `EC2_SSH_KEY` GitHub secrets.

### View live logs
```bash
# PR container
docker logs atlas-pr-<PR_NUMBER> --follow --tail 200

# Production containers
docker-compose logs --follow --tail 200 shiny-app-1
docker logs tcell_atlas_app_1 --follow --tail 200
```

### Shell into a running container
```bash
docker exec -it atlas-pr-<PR_NUMBER> bash
docker exec -it tcell_atlas_app_1 bash
```

### Check data files are present (most common failure cause)
```bash
# From inside the container
ls /srv/shiny-server/atlas/data/

# From the host
ls /home/ubuntu/atlas/data/ | grep -i <dataset_prefix>
```

Missing files cause `readRDS()` in `server.R` line ~344 to throw without a user-visible error —
the reactive observer silently dies. Add `tryCatch` around `readRDS` calls to expose errors in the UI.

### Check which PR containers are running
```bash
docker ps --filter "name=atlas-pr-" --format "table {{.Names}}\t{{.Ports}}\t{{.Status}}"
```

### Inspect Shiny server logs (inside container)
```bash
ls /var/log/shiny-server/
cat /var/log/shiny-server/atlas-shiny-*.log | tail -100
```

### Dataset load failures: diagnostic checklist
1. **File missing**: `ls /home/ubuntu/atlas/data/<filename>` — if absent, scp the file up.
2. **Filename mismatch**: compare `dataset_details.json` entry with actual filenames (case-sensitive).
3. **Wrong data level key**: `data_levels` values in JSON must match the keys in `files` (e.g., `"fib"` maps to `files.fib`).
4. **No `gene_list`**: subset levels need a `gene_list` file; a missing entry causes `readRDS(gene_list_path)` to fail on `NA` path.
5. **VAMcdf assay absent**: `server.R` line ~347 accesses `@assays$VAMcdf` — subset objects must have this assay.

### Authentication (SQLite)
Users DB: `/home/ubuntu/atlas/data/users_current.sqlite` (host) → copied to Docker volume `users-db` on first deploy.
Force re-sync: trigger `deploy.yml` workflow with `force_db_sync: true`.

## Code Conventions

- `setup.R` builds `dataset_files[[id]][[level]][["seurat"]]` from `dataset_details.json`. All server file lookups follow this pattern.
- `data_levels` in JSON is a display-name → key map (e.g., `{"Fibroblasts": "fib"}`); an empty object `{}` means spatial-only (no scRNA selector shown).
- `spatial_seurat: null` in JSON means no spatial data for that dataset; omit the key or set it to a filename to enable spatial.
- DEG files live in `data/DE_dfs/` (the `DE_dir` prefix in `server.R`). Seurat RDS and gene list files are directly in `data/`.
