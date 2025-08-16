inDir <- "/srv/shiny-server/atlas/data/"
DE_dir <- "DE_dfs/"
cat("Working directory is:", getwd(), "\n", file = stderr())


dataset_meta <- read.delim(
  file.path("..", "config", "datasets.tsv"),
  sep = "\t",
  stringsAsFactors = FALSE
)

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
