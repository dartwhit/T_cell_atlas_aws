inDir <- "/srv/shiny-server/atlas/data/"
DE_dir <- "DE_dfs/"
cat("Working directory is:", getwd(), "\n", file = stderr())



# Define the file paths for your datasets
file_map <- read.csv("dataset_files.csv", stringsAsFactors = FALSE)
dataset_files <- split(file_map, file_map$dataset_id) |>
  lapply(function(df) {
    split(df, df$data_level) |>
      lapply(function(x) {
        files <- setNames(x$filename, x$file_type)
        if (length(files) == 1) unname(files) else files
      })
  })

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
