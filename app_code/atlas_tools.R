# Read-only query functions for the AI chatbot.
# Each function accepts a seurat_obj explicitly and returns JSON-serialisable data.

query_gene_expression <- function(seurat_obj, genes, group_by = "seurat_clusters") {
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) return(list(error = "None of the requested genes found in dataset."))

  mat <- Seurat::FetchData(seurat_obj, vars = c(genes, group_by))
  result <- mat |>
    dplyr::group_by(.data[[group_by]]) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(genes),
                                   list(mean = mean, pct_expressed = ~ mean(. > 0))),
                     .groups = "drop")
  as.data.frame(result)
}

# NOTE: FindAllMarkers can be slow (minutes) on large objects.
# Use cluster_id to restrict to one cluster when possible.
get_cluster_markers <- function(seurat_obj, cluster_id = NULL, n_genes = 10) {
  if (!is.null(cluster_id)) {
    other_cells <- colnames(seurat_obj)[Seurat::Idents(seurat_obj) != cluster_id]
    markers <- Seurat::FindMarkers(
      seurat_obj,
      ident.1         = cluster_id,
      cells.2         = other_cells,
      only.pos        = TRUE,
      min.pct         = 0.25,
      logfc.threshold = 0.25,
      verbose         = FALSE
    )
    markers$gene    <- rownames(markers)
    markers$cluster <- cluster_id
    head(markers[order(-markers$avg_log2FC), ], n_genes)
  } else {
    markers <- Seurat::FindAllMarkers(
      seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE
    )
    markers |>
      dplyr::group_by(cluster) |>
      dplyr::slice_max(order_by = avg_log2FC, n = n_genes) |>
      as.data.frame()
  }
}

get_dataset_summary <- function(seurat_obj) {
  list(
    n_cells       = ncol(seurat_obj),
    n_genes       = nrow(seurat_obj),
    clusters      = levels(Seurat::Idents(seurat_obj)),
    n_clusters    = length(levels(Seurat::Idents(seurat_obj))),
    metadata_cols = colnames(seurat_obj@meta.data)
  )
}

find_cells_by_metadata <- function(seurat_obj, column, value) {
  if (!column %in% colnames(seurat_obj@meta.data)) {
    return(list(error = paste("Column", column, "not found in metadata.")))
  }
  cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[column]] == value, ])
  list(n_cells = length(cells), sample_cells = head(cells, 10))
}

check_gene_exists <- function(seurat_obj, gene) {
  list(
    gene    = gene,
    exists  = gene %in% rownames(seurat_obj),
    similar = if (!gene %in% rownames(seurat_obj))
                head(agrep(gene, rownames(seurat_obj), value = TRUE, max.distance = 0.2), 5)
              else character(0)
  )
}
