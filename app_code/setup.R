inDir <- "data/"
DE_dir <- "DE_dfs/"
cat("Working directory is:", getwd(), "\n", file = stderr())



dataset_meta <- read.delim(
  file.path("config", "datasets.tsv"),
  sep = "\t",
  stringsAsFactors = FALSE
)
dataset_meta$has_scrna <- as.logical(dataset_meta$has_scrna)
dataset_meta$has_spatial <- as.logical(dataset_meta$has_spatial)

# Minimal mapping of dataset ids to available files

dataset_files <- setNames(
  lapply(dataset_meta$file_path, function(fp) {
    list(full = list(seurat = fp))
  }),
  dataset_meta$id
)

# Data level choices limited to full dataset

data_level_choices <- setNames(
  rep(list(c("Full" = "full")), nrow(dataset_meta)),
  dataset_meta$id
)

inDir <- "data/"
DE_dir <- "DE_dfs/"
cat("Working directory is:", getwd(), "\n", file = stderr())



# Define the file paths for your datasets
dataset_files <- list(
  # --------- TMKMH dataset ------------
  "tmkmh" = list(
    # "meta" = "",
    "full" = list(
      "seurat" = "TMKMH_filtered_reduced_2025-06-06_VAM.rds",
      "gene_list" = "TMKMH_gene_list.rds",
      "DEGs_auto" = "TMKMH_auto_DEGs_sig.txt",
      "DEGs_broad" = "TMKMH_broad_DEGs_sig.txt",
      "VAM_df" = "TMKMH_full_filtered_VAM_top.txt",
      "DE_by_disease_broad" = "TMKMH_full_broad_DE_by_disease.txt",
      "DE_by_disease_auto" = "TMKMH_full_auto_DE_by_disease.txt"
    ),
    "fib" = list(
      "seurat" = "TMKMH_fib_filtered_VAM.rds",
      "gene_list" = "TMKMH_fib_gene_list.rds",
      "DEGs" = "TMKMH_fib_DEGs_sig.txt",
      "VAM_df" = "TMKMH_fib_filtered_VAM_top.txt",
      "DE_by_disease" = "TMKMH_fib_DE_by_disease.txt",
      "VAM_by_disease" = "TMKMH_fib_VAMcdf_DE_by_disease.txt"
    ),
    "immune" = list(
      "seurat" = "TMKMH_immune_filtered_VAM.rds",
      "gene_list" = "TMKMH_immune_gene_list.rds",
      "DEGs" = "TMKMH_immune_DEGs_sig.txt",
      "VAM_df" = "TMKMH_immune_filtered_VAM_top.txt",
      "DE_by_disease" = ""
    ),
    "spatial_seurat" = NULL
  )
  ,
  # ----------- Tabib dataset ---------
  "tabib" = list(
    "meta" = "Tabib_metadata.txt",
    "full" = list(
      "seurat" = "Tabib_full_filtered_VAM.rds",
      "gene_list" = "Tabib_gene_list.rds",
      "DEGs_auto" = "Tabib_auto_DEGs.csv",
      "DEGs_broad" = "Tabib_broad_DEGs_sig.txt",
      "VAM_df" = "Tabib_full_filtered_VAM_top.txt",
      "DE_by_disease_broad" = "Tabib_broad_DE_by_disease.txt",
      "DE_by_disease_auto" = "Tabib_auto_DE_by_disease.txt",
      "VAM_by_disease_broad" = "Tabib_broad_VAMcdf_DE_by_disease.txt",
      "VAM_by_disease_auto" = "Tabib_auto_VAMcdf_DE_by_disease.txt"
    ),
    "fib" = list(
      "seurat" = "Tabib_fib_filtered_VAM.rds",
      "gene_list" = "Tabib_fib_gene_list.rds",
      "DEGs" = "Tabib_fib_DEGs_sig.txt",
      "VAM_df" = "Tabib_fib_filtered_VAM_top.txt",
      "DE_by_disease" = "Tabib_fib_DE_by_disease.txt",
      "VAM_by_disease" = "Tabib_fib_VAMcdf_DE_by_disease.txt"
    ),
    "immune" = list(
      "seurat" = "Tabib_immune_filtered_VAM.rds",
      "gene_list" = "Tabib_immune_gene_list.rds",
      "DEGs" = "Tabib_immune_DEGs_sig.txt",
      "VAM_df" = "Tabib_immune_filtered_VAM_top.txt",
      "DE_by_disease" = "Tabib_immune_DE_by_disease.txt",
      "VAM_by_disease" = "Tabib_immune_VAMcdf_DE_by_disease.txt"
    ),
    "spatial_seurat" = NULL
  ),
  # ------------- Gur Dataset ------------
  "gur" = list(
    "meta" = "Gur_metadata.txt",
    "full" = list(
      "seurat" = "Gur_full_filtered_VAM.rds",
      "gene_list" = "Gur_gene_list.rds",
      "DEGs_auto" = "Gur_auto_DEGs_sig.txt",
      "DEGs_broad" = "Gur_broad_DEGs_sig.txt",
      "VAM_df" = "Gur_full_filtered_VAM_top.txt",
      "DE_by_disease_broad" = "Gur_full_broad_DE_by_disease.txt",
      "DE_by_disease_auto" = "Gur_full_auto_DE_by_disease.txt",
      "VAM_by_disease_broad" = "Gur_broad_VAMcdf_DE_by_disease.txt",
      "VAM_by_disease_auto" = "Gur_auto_VAMcdf_DE_by_disease.txt"
      
    ),
    "fib" = list(
      "seurat" = "Gur_fib_filtered_VAM.rds",
      "gene_list" = "Gur_fib_gene_list.rds",
      "DEGs" = "Gur_fib_DEGs_sig.txt",
      "VAM_df" = "Gur_fib_filtered_VAM_top.txt",
      "DE_by_disease" = "Gur_fib_DE_by_disease.txt",
      "VAM_by_disease" = "Gur_fib_VAMcdf_DE_by_disease.txt"
    ),
    "immune" = list(
      "seurat" = "Gur_immune_filtered_VAM.rds",
      "gene_list" = "Gur_immune_gene_list.rds",
      "DEGs" = "Gur_immune_DEGs_sig.txt",
      "VAM_df" = "Gur_immune_filtered_VAM_top.txt",
      "DE_by_disease" = "Gur_immune_DE_by_disease.txt",
      "VAM_by_disease" = "Gur_immune_VAMcdf_DE_by_disease.txt"
      
    ),
    "spatial_seurat" = NULL
  ),
  # -------------- Ma dataset --------------
  "ma" = list(
    "full" = list(
      "seurat" = "Ma_filtered_reduced_2025-01-28_VAM.rds",
      "gene_list" = "Ma_gene_list.rds",
      "DEGs_auto" = "Ma_auto_DEGs_sig.txt",
      "DEGs_broad" = "Ma_broad_DEGs_sig.txt",
      "VAM_df" = "Ma_filtered_reduced_2025-01-28_VAM_top.txt",
      "DE_by_disease_broad" = "Ma_full_broad_DE_by_disease.txt",
      "DE_by_disease_auto" = "Ma_full_auto_DE_by_disease.txt"
    ),
    "mye" = list(
      "seurat" = "Ma_mye_filtered_VAM.rds",
      "gene_list" = "Ma_mye_gene_list.rds",
      "DEGs" = "Ma_mye_DEGs_sig.txt",
      "VAM_df" = "Ma_mye_filtered_VAM_top.txt",
      "DE_by_disease" = "Ma_mye_DE_by_disease.txt"
    ),
    "fib" = list(
      "seurat" = "Ma_fib_filtered_VAM.rds",
      "gene_list" = "Ma_fib_gene_list.rds",
      "DEGs" = "Ma_fib_DEGs_sig.txt",
      "VAM_df" = "Ma_fib_filtered_VAM_top.txt",
      "DE_by_disease" = "Ma_fib_DE_by_disease.txt"
    ),
    "spatial_seurat" = "2025-07-07_MaSSc_Visium_PRECAST_SingleCellPredicted_RegionsNamed_CARD.rds"
  ),
  # ------------- Khanna Dataset -------------
  "khanna" = list(
    "full" = list(
      "seurat" = "Khanna_filtered_reduced_2025-01-28_VAM.rds",
      "gene_list" = "Khanna_gene_list.rds",
      "DEGs_auto" = "Khanna_auto_DEGs_sig.txt",
      "DEGs_broad" = "Khanna_broad_DEGs_sig.txt",
      "VAM_df" = "Khanna_filtered_reduced_2025-01-28_VAM_top.txt"
    ),
    "fib" = list(
      "seurat" = "Khanna_fib_filtered_VAM.rds",
      "gene_list" = "Khanna_fib_gene_list.rds",
      "DEGs" = "Khanna_fib_DEGs_sig.txt",
      "VAM_df" = "Khanna_fib_filtered_VAM_top.txt"
    ),
    "mye" = list(
      "seurat" = "Khanna_mye_filtered_VAM.rds",
      "gene_list" = "Khanna_mye_gene_list.rds",
      "DEGs" = "Khanna_mye_DEGs_sig.txt",
      "VAM_df" = "Khanna_mye_filtered_VAM_top.txt"
      
    ),
    "spatial_seurat" = NULL
  )
)


data_level_choices <- list(
  tabib = c("Full" = "full", 
            "Fibroblasts" = "fib", 
            "Immune cells" = "immune"),
  gur = c("Full" = "full",
          "Fibroblasts" = "fib", 
          "Immune cells" = "immune"),
  ma = c("Full" = "full", 
         "Myeloid cells" = "mye", 
         "Fibroblasts" = "fib"),
  khanna = c("Full" = "full",
             "Myeloid cells" = "mye", 
             "Fibroblasts" = "fib"), # Example: only "Full" data available for Clark
  tmkmh = c("Full" = "full",
            "Fibroblasts" = "fib",
            "Immune cells" = "immune")
)

# Create datasets object for sidebar module (this was missing!)
datasets <- dataset_meta