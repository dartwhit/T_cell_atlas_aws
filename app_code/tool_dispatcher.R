dispatch_tool <- function(tool_name, tool_input, seurat_obj) {
  tryCatch({
    result <- switch(tool_name,
      query_gene_expression = query_gene_expression(
        seurat_obj,
        genes    = tool_input$genes,
        group_by = tool_input$group_by %||% "seurat_clusters"
      ),
      get_cluster_markers = get_cluster_markers(
        seurat_obj,
        cluster_id = tool_input$cluster_id,
        n_genes    = tool_input$n_genes %||% 10L
      ),
      get_dataset_summary   = get_dataset_summary(seurat_obj),
      find_cells_by_metadata = find_cells_by_metadata(
        seurat_obj,
        column = tool_input$column,
        value  = tool_input$value
      ),
      check_gene_exists = check_gene_exists(seurat_obj, gene = tool_input$gene),
      list(error = paste("Unknown tool:", tool_name))
    )
    jsonlite::toJSON(result, auto_unbox = TRUE, na = "null")
  }, error = function(e) {
    jsonlite::toJSON(list(error = conditionMessage(e)), auto_unbox = TRUE)
  })
}
