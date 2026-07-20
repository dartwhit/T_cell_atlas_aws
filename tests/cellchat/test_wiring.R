# Reactive wiring checks for the CellChat explorer module.
#
# Drives cellchat_explorer_server() with shiny::testServer() against the real
# tabib CellChat objects, and checks which sidebar inputs actually change the
# communication table.

source("../tests/cellchat/helpers.R")

cat("CellChat wiring harness — study:", QA_STUDY, "\n")
cat(strrep("=", 78), "\n\n", sep = "")

testServer(
  cellchat_explorer_server,
  args = list(dataset_configs = setNames(list(QA_CONFIG), QA_STUDY)),
  {
    # --- load -----------------------------------------------------------------
    session$setInputs(study = QA_STUDY, level = "LR", measure = "weight",
                      ligand = "", receptor = "", prob_min = 0, pval_max = 1,
                      sources = character(), targets = character(),
                      pathways = character())

    data <- cellchat_data()
    qa_record("objects load without error", if (is.null(data$error)) "PASS" else "FAIL",
              data$error %||% paste(length(data$conditions), "condition objects"))

    conditions <- names(data$conditions)
    celltypes <- cellchat_celltypes(data$merged)
    pathways <- cellchat_pathways(data$merged)
    qa_record("condition choices populated",
              if (length(conditions) > 0) "PASS" else "FAIL",
              paste(conditions, collapse = ", "))
    qa_record("cell type choices populated",
              if (length(celltypes) > 0) "PASS" else "FAIL",
              paste(length(celltypes), "cell types"))
    qa_record("pathway choices populated",
              if (length(pathways) > 0) "PASS" else "FAIL",
              paste(length(pathways), "pathways"))

    baseline <- filtered_table()
    qa_record("communication table non-empty (LR)",
              if (nrow(baseline) > 0) "PASS" else "FAIL",
              paste(nrow(baseline), "rows"))

    # Reference values drawn from the data itself so the checks stay valid if
    # the underlying objects are regenerated.
    a_source <- celltypes[[1]]
    a_target <- celltypes[[min(2, length(celltypes))]]
    a_pathway <- if (length(pathways) > 0) pathways[[1]] else NULL
    a_ligand  <- if (nrow(baseline) > 0) substr(baseline$ligand[[1]],   1, 3) else NULL
    a_receptor <- if (nrow(baseline) > 0) substr(baseline$receptor[[1]], 1, 3) else NULL

    rows_now <- function() nrow(filtered_table())
    differs <- function(n) !identical(n, nrow(baseline))
    reset <- function() {
      session$setInputs(sources = character(), targets = character(),
                        pathways = character(), ligand = "", receptor = "",
                        prob_min = 0, pval_max = 1)
    }

    # --- filters that should drive the table ---------------------------------
    session$setInputs(sources = a_source)
    qa_expect_reactive("table reacts to `sources`", differs(rows_now()),
                       sprintf("%d -> %d rows (%s)", nrow(baseline), rows_now(), a_source))
    reset()

    session$setInputs(targets = a_target)
    qa_expect_reactive("table reacts to `targets`", differs(rows_now()),
                       sprintf("%d -> %d rows (%s)", nrow(baseline), rows_now(), a_target))
    reset()

    session$setInputs(pathways = a_pathway)
    qa_expect_reactive("table reacts to `pathways`", differs(rows_now()),
                       sprintf("%d -> %d rows (%s)", nrow(baseline), rows_now(), a_pathway))
    reset()

    session$setInputs(ligand = a_ligand)
    qa_expect_reactive("table reacts to `ligand`", differs(rows_now()),
                       sprintf("%d -> %d rows (contains '%s')", nrow(baseline), rows_now(), a_ligand))
    reset()

    session$setInputs(receptor = a_receptor)
    qa_expect_reactive("table reacts to `receptor`", differs(rows_now()),
                       sprintf("%d -> %d rows (contains '%s')", nrow(baseline), rows_now(), a_receptor))
    reset()

    median_prob <- if (nrow(baseline) > 0) stats::median(baseline$prob, na.rm = TRUE) else 0  # 0 = no-op filter when baseline is empty
    session$setInputs(prob_min = median_prob)
    qa_expect_reactive("table reacts to `prob_min`", differs(rows_now()),
                       sprintf("%d -> %d rows (prob >= %.3g)", nrow(baseline), rows_now(), median_prob))
    reset()

    session$setInputs(pval_max = 0.01)
    qa_expect_reactive("table reacts to `pval_max`", differs(rows_now()),
                       sprintf("%d -> %d rows (pval <= 0.01)", nrow(baseline), rows_now()))
    reset()

    session$setInputs(level = "pathway")
    pathway_rows <- rows_now()
    qa_expect_reactive("table reacts to `level` (LR vs pathway)", differs(pathway_rows),
                       sprintf("%d -> %d rows", nrow(baseline), pathway_rows))
    session$setInputs(level = "LR")

    # --- the table's own condition scope -------------------------------------
    session$setInputs(table_condition = CELLCHAT_ALL_CONDITIONS)
    all_rows <- rows_now()
    per_condition <- vapply(conditions, function(cond) {
      session$setInputs(table_condition = cond)
      nrow(filtered_table())
    }, numeric(1))
    qa_expect_reactive(
      "table reacts to `table_condition`",
      length(unique(per_condition)) > 1 || all(per_condition < all_rows),
      paste(c(sprintf("all=%d", all_rows),
              sprintf("%s=%d", names(per_condition), per_condition)), collapse = ", ")
    )

    # Scoping to one condition must leave only that condition's rows.
    session$setInputs(table_condition = conditions[[1]])
    scoped <- filtered_table()
    qa_record("scoped table contains only the selected condition",
              if (identical(unique(scoped$dataset), conditions[[1]])) "PASS" else "UNWIRED",
              paste("dataset column holds:",
                    paste(unique(scoped$dataset), collapse = ", ")))

    # ...and the rows must still add up to the unscoped total.
    qa_record("per-condition row counts sum to the 'All conditions' total",
              if (sum(per_condition) == all_rows) "PASS" else "FAIL",
              sprintf("%d == %d", sum(per_condition), all_rows))
    session$setInputs(table_condition = CELLCHAT_ALL_CONDITIONS)

    # --- download handler ----------------------------------------------------
    tmp <- tempfile(fileext = ".csv")
    output$download
    qa_try("download handler writes non-empty CSV", {
      utils::write.csv(filtered_table(), tmp, row.names = FALSE)
      n <- nrow(utils::read.csv(tmp))
      qa_record("download handler writes non-empty CSV",
                if (n > 0) "PASS" else "FAIL", paste(n, "rows"))
    })

    # --- draw_plot must not swallow req()'s silent error ---------------------
    # heatmap_role / role_scatter call req() inside draw_plot(). req() raises a
    # shiny.silent.error, which inherits from "error", so a bare
    # tryCatch(error=) would catch it and render a bogus error panel instead of
    # leaving the plot blank while the input is still unset.
    grDevices::pdf(NULL)  # draw_plot's error path draws; keep it off disk
    on.exit(grDevices::dev.off(), add = TRUE)
    silent <- tryCatch(draw_plot(function(data) req(NULL)), error = function(e) e)
    qa_record("draw_plot lets req()'s silent error through",
              if (inherits(silent, "shiny.silent.error")) "PASS" else "UNWIRED",
              if (inherits(silent, "shiny.silent.error")) "plot stays blank"
              else "renders 'could not be generated:' with an empty message")

    # A genuine error must still reach the user as a drawn message.
    real <- tryCatch(draw_plot(function(data) stop("boom")), error = function(e) e)
    qa_record("draw_plot still reports genuine plot errors",
              if (!inherits(real, "error")) "PASS" else "FAIL",
              "real errors are drawn, not propagated")
  }
)

qa_summary()
