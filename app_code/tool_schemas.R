atlas_tool_schemas <- list(
  list(
    name         = "query_gene_expression",
    description  = "Returns mean expression and percent expressed for one or more genes, grouped by cluster or another metadata column. Use when the user asks about expression levels, which clusters express a gene, or co-expression.",
    input_schema = list(
      type       = "object",
      properties = list(
        genes    = list(type = "array", items = list(type = "string"),
                        description = "Gene symbols to query, e.g. ['CD3E', 'CD8A']."),
        group_by = list(type = "string",
                        description = "Metadata column to group by. Default: 'seurat_clusters'.")
      ),
      required = list("genes")
    )
  ),
  list(
    name         = "get_cluster_markers",
    description  = "Returns top differentially expressed marker genes for one or all clusters. Use when the user asks what defines a cluster or what genes are specific to a cell type. Providing cluster_id is much faster than requesting all clusters at once.",
    input_schema = list(
      type       = "object",
      properties = list(
        cluster_id = list(type = "string",
                          description = "Cluster identifier. Omit to return markers for all clusters (slow)."),
        n_genes    = list(type = "integer",
                          description = "Number of top markers to return per cluster. Default: 10.")
      ),
      required = list()
    )
  ),
  list(
    name         = "get_dataset_summary",
    description  = "Returns basic statistics: cell count, gene count, cluster labels, and available metadata columns. Call this first when the user's question is vague or you need to know what is available.",
    input_schema = list(
      type       = "object",
      properties = list(),
      required   = list()
    )
  ),
  list(
    name         = "find_cells_by_metadata",
    description  = "Counts and identifies cells matching a specific value in a metadata column (e.g. sample == 'Patient_01', Disease == 'SSc').",
    input_schema = list(
      type       = "object",
      properties = list(
        column = list(type = "string", description = "Metadata column name."),
        value  = list(type = "string", description = "Value to filter by.")
      ),
      required = list("column", "value")
    )
  ),
  list(
    name         = "check_gene_exists",
    description  = "Checks whether a gene symbol exists in the dataset and suggests similar names if not. Always call this before querying an unfamiliar gene.",
    input_schema = list(
      type       = "object",
      properties = list(
        gene = list(type = "string", description = "Gene symbol to look up.")
      ),
      required = list("gene")
    )
  )
)
