# To run this test script:
# 1. Make sure you are in the root directory of the T_cell_atlas_aws project.
# 2. Open R.
# 3. Run the command: shiny::runApp("app_code/modules/test_scrna_seq_module.R")
#
# This script creates a self-contained Shiny app to test the scRNA-seq module.
# It generates a lightweight, mock Seurat object and all necessary input files,
# placing them in a temporary directory that is cleaned up after the app closes.

library(shiny)
library(shinyjs)
library(Seurat)
library(Matrix)
library(dplyr)

# --- Setup Temporary Test Environment ---

# Create a temporary directory for mock data. This is cleaned up automatically.
test_data_dir <- file.path(tempdir(), "scrna_seq_module_test")
if (dir.exists(test_data_dir)) {
  unlink(test_data_dir, recursive = TRUE) # Clean up from previous runs
}
dir.create(test_data_dir, recursive = TRUE)
dir.create(file.path(test_data_dir, "DE_dfs"), recursive = TRUE)


# --- Generate Mock Data ---

# 1. Create a Mock Seurat Object
set.seed(42)
n_genes <- 100
n_cells <- 50
mock_counts <- Matrix(
  rpois(n_genes * n_cells, lambda = 2),
  nrow = n_genes,
  ncol = n_cells,
  sparse = TRUE
)
rownames(mock_counts) <- paste0("Gene", 1:n_genes)
colnames(mock_counts) <- paste0("Cell", 1:n_cells)
mock_metadata <- data.frame(
  cell_type = sample(c("T-cell", "B-cell", "Macrophage"), n_cells, replace = TRUE),
  row.names = colnames(mock_counts)
)
seurat_obj <- CreateSeuratObject(
  counts = mock_counts,
  meta.data = mock_metadata
)

# --- Make the Mock Seurat Object more realistic ---

# 1. Add a mock VAMcdf assay for pathway analysis
n_pathways <- 20
mock_vam_data <- matrix(
  runif(n_pathways * n_cells),
  nrow = n_pathways,
  ncol = n_cells
)
rownames(mock_vam_data) <- paste0("HALLMARK-PATHWAY", 1:n_pathways)
colnames(mock_vam_data) <- colnames(seurat_obj)
seurat_obj[["VAMcdf"]] <- CreateAssayObject(counts = mock_vam_data)

# 2. Add dimensional reduction for plotting
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 30, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)


# Save the mock seurat object to the temp directory for the module to access.
mock_seurat_path <- file.path(test_data_dir, "mock_seurat.rds")
saveRDS(seurat_obj, mock_seurat_path)

# 2. Extract Gene List from Mock Data
mock_gene_list <- rownames(seurat_obj)
mock_gene_list_path <- file.path(test_data_dir, "mock_gene_list.rds")
saveRDS(mock_gene_list, mock_gene_list_path)

# 3. Mock DEGs Table (using genes from the mock object)
mock_degs <- data.frame(
  cluster = rep(c("Cluster1", "Cluster2"), each = 5),
  gene = sample(rownames(seurat_obj), 10),
  avg_log2FC = rnorm(10),
  p_val = runif(10, 0, 0.05),
  p_val_adj = runif(10, 0, 0.1)
)
mock_degs_path <- file.path(test_data_dir, "DE_dfs/mock_degs.csv")
write.csv(mock_degs, mock_degs_path, row.names = FALSE)


# --- Mock Global Variables (from setup.R) ---
# These variables are normally created in the main app's setup.R.
# Here, we mock them to point to our temporary test data.

# Mock the 'datasets' data frame that is normally read from datasets.tsv
datasets <- data.frame(
  id = c("test_study"),
  name = c("Test Study")
)

inDir <- paste0(test_data_dir, "/") # The module expects a trailing slash
DE_dir <- "DE_dfs/"

dataset_files <- list(
  "test_study" = list(
    "full" = list(
      "seurat" = "mock_seurat.rds",
      "gene_list" = "mock_gene_list.rds",
      "DEGs_auto" = "mock_degs.csv",
      "DEGs_broad" = "mock_degs.csv",
      "DE_by_disease_auto" = "mock_degs.csv",
      "DE_by_disease_broad" = "mock_degs.csv"
    )
  )
)

data_level_choices <- list(
  "test_study" = c("Full" = "full")
)


# --- Source the Modules ---
# Use the 'here' package to build paths relative to the project root,
# ensuring that the script can be run from any working directory.
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "http://cran.us.r-project.org")
}
library(here)

source(here("app_code", "modules", "explore_sidebar_module.R"))
source(here("app_code", "modules", "scrna_seq_module.R"))


# --- Test App UI ---
ui <- fluidPage(
  useShinyjs(), # Initialize shinyjs
  title = "Test scRNA-seq Module with Mock Data",
  layout_sidebar(
    sidebar = explore_sidebar_UI("explore_sidebar_module"),
    scrna_seq_UI("scrna_seq_module")
  )
)

# --- Test App Server ---
server <- function(input, output, session) {
  # Simulate the gallery page selecting a study.
  mock_selected_study <- reactiveVal(list(id = "test_study"))
  
  # Call the sidebar server module
  sidebar_inputs <- explore_sidebar_server(
    "explore_sidebar_module", 
    selected_study_from_gallery = mock_selected_study,
    dataset_config = datasets
  )
  
  # Call the main scRNA-seq server module
  scrna_data <- scrna_seq_server(
    "scrna_seq_module", 
    sidebar_inputs = sidebar_inputs,
    dataset_files = dataset_files,
    inDir = inDir
  )
  
  # Observer to update the sidebar UI when the gene/pathway lists are loaded
  observe({
    req(scrna_data$gene_list())
    updateSelectizeInput(
      session, 
      "explore_sidebar_module-gene_select", 
      choices = scrna_data$gene_list(), 
      server = TRUE
    )
  })
  
  observe({
    req(scrna_data$pathway_list())
    updateSelectInput(
      session,
      "explore_sidebar_module-pathway_select",
      choices = scrna_data$pathway_list()
    )
  })
  
  # Clean up the temporary directory when the app is closed.
  session$onSessionEnded(function() {
    unlink(test_data_dir, recursive = TRUE)
    cat("Cleaned up temporary test data directory:", test_data_dir, "\n")
  })
}

# --- Run the App ---
shinyApp(ui, server)