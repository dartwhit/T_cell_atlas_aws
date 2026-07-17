# Utilities for loading and visualizing precomputed CellChat results.
#
# The app uses stripped CellChat objects that retain network and annotation
# slots while omitting expression matrices. This keeps the interactive views
# light enough to load on demand for an individual study.

cellchat_merged_datasets <- function(obj) {
  if (is.null(obj) || !methods::is(obj, "CellChat")) return(character())
  net <- obj@net
  nm <- names(net)
  if (is.list(net) && length(nm) > 0 && !("count" %in% nm)) nm else character()
}

cellchat_is_merged <- function(obj) length(cellchat_merged_datasets(obj)) > 0

cellchat_celltypes <- function(obj) {
  if (is.null(obj) || !methods::is(obj, "CellChat")) return(character())
  idents <- obj@idents
  if (is.list(idents) && !is.null(idents$joint)) {
    levels(idents$joint)
  } else {
    levels(idents)
  }
}

cellchat_pathways <- function(obj) {
  if (is.null(obj) || !methods::is(obj, "CellChat")) return(character())
  netp <- obj@netP
  if (cellchat_is_merged(obj)) {
    datasets <- cellchat_merged_datasets(obj)
    pathways <- unlist(lapply(datasets, function(dataset) {
      netp[[dataset]]$pathways %||% character()
    }), use.names = FALSE)
    sort(unique(pathways))
  } else {
    sort(unique(netp$pathways %||% character()))
  }
}

cellchat_load_objects <- function(study, config) {
  empty <- list(merged = NULL, conditions = list(), error = NULL)

  if (!requireNamespace("CellChat", quietly = TRUE)) {
    empty$error <- "The CellChat package is not available in this app image."
    return(empty)
  }
  if (!cellchat_config_complete(config)) {
    empty$error <- "This study does not have a complete CellChat configuration."
    return(empty)
  }

  read_configured_object <- function(relative_path) {
    tryCatch(
      read_object(data_path(study, relative_path)),
      error = function(error) error
    )
  }

  merged <- read_configured_object(config$merged)
  if (inherits(merged, "error")) {
    empty$error <- paste("Unable to load the merged CellChat object:", conditionMessage(merged))
    return(empty)
  }
  if (!methods::is(merged, "CellChat") || !cellchat_is_merged(merged)) {
    empty$error <- "The configured merged object is not a valid multi-condition CellChat object."
    return(empty)
  }

  configured_conditions <- config$conditions
  condition_objects <- lapply(configured_conditions, read_configured_object)
  failures <- names(condition_objects)[vapply(condition_objects, inherits, logical(1), what = "error")]
  if (length(failures) > 0) {
    messages <- vapply(condition_objects[failures], conditionMessage, character(1))
    empty$error <- paste0(
      "Unable to load CellChat condition object(s) ",
      paste(failures, collapse = ", "), ": ",
      paste(messages, collapse = "; ")
    )
    return(empty)
  }
  invalid <- names(condition_objects)[!vapply(condition_objects, methods::is, logical(1), class2 = "CellChat")]
  if (length(invalid) > 0) {
    empty$error <- paste("Configured condition object(s) are not CellChat objects:",
                         paste(invalid, collapse = ", "))
    return(empty)
  }

  list(merged = merged, conditions = condition_objects, error = NULL)
}

# ---- Communication table extraction and filtering ---------------------------

cellchat_standardize_comm_columns <- function(data) {
  if (!is.data.frame(data) || nrow(data) == 0) return(tibble::tibble())
  out <- tibble::as_tibble(data)
  candidates <- list(
    source = c("source", "source.labels", "sender"),
    target = c("target", "target.labels", "receiver"),
    ligand = c("ligand", "geneL"),
    receptor = c("receptor", "geneR"),
    pathway = c("pathway_name", "pathway", "signaling"),
    prob = c("prob", "probability", "score"),
    pval = c("pval", "p.value", "p_value", "pval.adj")
  )

  for (name in names(candidates)) {
    match <- candidates[[name]][candidates[[name]] %in% names(out)]
    if (length(match) > 0 && match[[1]] != name) {
      names(out)[names(out) == match[[1]]] <- name
    }
  }
  for (name in names(candidates)) {
    if (!(name %in% names(out))) out[[name]] <- NA
  }

  out <- dplyr::mutate(
    out,
    source = as.character(.data$source),
    target = as.character(.data$target),
    ligand = as.character(.data$ligand),
    receptor = as.character(.data$receptor),
    pathway = as.character(.data$pathway),
    prob = suppressWarnings(as.numeric(.data$prob)),
    pval = suppressWarnings(as.numeric(.data$pval)),
    pair = ifelse(is.na(.data$ligand) | is.na(.data$receptor), NA_character_,
                  paste(.data$ligand, .data$receptor, sep = " -> "))
  )

  preferred <- c("dataset", "source", "target", "ligand", "receptor", "pair",
                 "pathway", "prob", "pval")
  out[, c(intersect(preferred, names(out)), setdiff(names(out), preferred)), drop = FALSE]
}

cellchat_extract_comm <- function(obj, level = c("LR", "pathway")) {
  level <- match.arg(level)
  if (is.null(obj) || !requireNamespace("CellChat", quietly = TRUE)) {
    return(tibble::tibble())
  }
  slot_name <- if (level == "pathway") "netP" else "net"
  raw <- tryCatch(
    CellChat::subsetCommunication(obj, slot.name = slot_name),
    error = function(error) NULL
  )
  if (is.null(raw)) return(tibble::tibble())
  if (is.data.frame(raw)) return(cellchat_standardize_comm_columns(raw))
  if (!is.list(raw)) return(tibble::tibble())

  parts <- lapply(names(raw), function(dataset) {
    out <- cellchat_standardize_comm_columns(raw[[dataset]])
    if (nrow(out) > 0) tibble::add_column(out, dataset = dataset, .before = 1) else out
  })
  dplyr::bind_rows(parts)
}

cellchat_contains_pattern <- function(values, pattern) {
  if (is.null(pattern) || !nzchar(pattern)) return(rep(TRUE, length(values)))
  stringr::str_detect(
    tolower(as.character(values)),
    stringr::fixed(tolower(pattern))
  )
}

cellchat_apply_comm_filters <- function(data, sources = NULL, targets = NULL,
                                        ligand = NULL, receptor = NULL, pathway = NULL,
                                        prob_min = 0, pval_max = 1) {
  if (!is.data.frame(data) || nrow(data) == 0) return(tibble::as_tibble(data))
  out <- data
  matches_selection <- function(values, selected) {
    if (is.null(selected) || length(selected) == 0 || identical(selected, "All")) {
      rep(TRUE, length(values))
    } else {
      values %in% selected
    }
  }

  out <- out[matches_selection(out$source, sources), , drop = FALSE]
  out <- out[matches_selection(out$target, targets), , drop = FALSE]
  out <- out[matches_selection(out$pathway, pathway), , drop = FALSE]
  out <- out[cellchat_contains_pattern(out$ligand, ligand), , drop = FALSE]
  out <- out[cellchat_contains_pattern(out$receptor, receptor), , drop = FALSE]
  out <- out[is.na(out$prob) | out$prob >= prob_min, , drop = FALSE]
  out <- out[is.na(out$pval) | out$pval <= pval_max, , drop = FALSE]
  tibble::as_tibble(out)
}

# ---- Plot wrappers -----------------------------------------------------------

cellchat_selection_or_null <- function(value) {
  if (is.null(value) || length(value) == 0 || all(!nzchar(value))) NULL else value
}

cellchat_plot_circle_comparison <- function(merged, measure = c("weight", "count")) {
  measure <- match.arg(measure)
  datasets <- cellchat_merged_datasets(merged)
  matrices <- lapply(datasets, function(dataset) merged@net[[dataset]][[measure]])
  maximum <- max(vapply(matrices, function(matrix) max(matrix, na.rm = TRUE), numeric(1)))
  old_parameters <- graphics::par(mfrow = c(1, length(datasets)), xpd = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  for (index in seq_along(datasets)) {
    CellChat::netVisual_circle(
      matrices[[index]], weight.scale = TRUE, label.edge = FALSE,
      edge.weight.max = maximum,
      title.name = paste0(datasets[[index]], " (", measure, ")")
    )
  }
  invisible(NULL)
}

cellchat_plot_diff_interaction <- function(merged, measure = c("count", "weight")) {
  CellChat::netVisual_diffInteraction(
    merged, weight.scale = TRUE, measure = match.arg(measure)
  )
  invisible(NULL)
}

cellchat_plot_compare_interactions <- function(merged, measure = c("count", "weight")) {
  plot <- CellChat::compareInteractions(
    merged, show.legend = FALSE, group = seq_along(cellchat_merged_datasets(merged)),
    measure = match.arg(measure)
  )
  print(plot)
  invisible(NULL)
}

cellchat_plot_heatmap_diff <- function(merged, measure = c("count", "weight")) {
  heatmap <- CellChat::netVisual_heatmap(merged, measure = match.arg(measure))
  ComplexHeatmap::draw(heatmap)
  invisible(NULL)
}

cellchat_plot_signaling_role_heatmap <- function(obj, pattern = c("outgoing", "incoming", "all"),
                                                  signaling = NULL) {
  heatmap <- CellChat::netAnalysis_signalingRole_heatmap(
    obj, pattern = match.arg(pattern), signaling = cellchat_selection_or_null(signaling)
  )
  ComplexHeatmap::draw(heatmap)
  invisible(NULL)
}

cellchat_plot_bubble <- function(merged, sources = NULL, targets = NULL, signaling = NULL) {
  plot <- CellChat::netVisual_bubble(
    merged,
    sources.use = cellchat_selection_or_null(sources),
    targets.use = cellchat_selection_or_null(targets),
    signaling = cellchat_selection_or_null(signaling),
    comparison = seq_along(cellchat_merged_datasets(merged)),
    angle.x = 45
  )
  print(plot)
  invisible(NULL)
}

cellchat_plot_ranknet <- function(merged, stacked = TRUE) {
  plot <- CellChat::rankNet(
    merged, mode = "comparison",
    comparison = seq_along(cellchat_merged_datasets(merged)),
    stacked = stacked, do.stat = TRUE
  )
  print(plot)
  invisible(NULL)
}

cellchat_plot_signaling_role_scatter <- function(obj, signaling = NULL) {
  plot <- CellChat::netAnalysis_signalingRole_scatter(
    obj, signaling = cellchat_selection_or_null(signaling)
  )
  print(plot)
  invisible(NULL)
}

cellchat_draw_plot_error <- function(message) {
  graphics::plot.new()
  graphics::text(
    0.5, 0.5,
    labels = paste(strwrap(message, width = 70), collapse = "\n"),
    cex = 1.1
  )
}
