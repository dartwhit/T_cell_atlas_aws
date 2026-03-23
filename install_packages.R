## R package installation script used during Docker build.
## Called via: Rscript install_packages.R
## Stops with a non-zero exit code if any package fails to install,
## so `docker build` fails loudly instead of silently skipping packages.

install_required <- function(pkgs, repos = "https://cloud.r-project.org/") {
  install.packages(pkgs, repos = repos)
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing)) {
    stop("Failed to install the following packages: ", paste(missing, collapse = ", "))
  }
  invisible(NULL)
}

cat("=== Installing core Shiny / utility packages ===\n")
install_required(c(
  "plotly", "igraph", "Matrix", "DT", "ggplot2", "dplyr", "stringr",
  "shiny", "shinymanager", "shinyWidgets", "bslib", "shinycssloaders",
  "bsicons", "periscope2", "shinyjs", "tidyr", "DBI", "RSQLite",
  "jsonlite", "png"
))

cat("=== Installing Seurat ecosystem ===\n")
install_required(c("SeuratObject", "Seurat"))

cat("=== Installing VAM ===\n")
install_required("VAM")

cat("=== All packages installed successfully ===\n")
