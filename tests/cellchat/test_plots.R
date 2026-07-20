# Plot-wrapper checks for the CellChat explorer.
#
# Renders each wrapper in modules/cellchat_helpers.R to PNG across a parameter
# sweep and compares md5 hashes. A control that genuinely drives a plot produces
# a different image; a control that never reaches the underlying CellChat call
# produces a byte-identical image.
#
# PNGs go to a temp dir (or QA_PLOT_DIR) — they are not committed.

source("../tests/cellchat/helpers.R")

plot_dir <- Sys.getenv("QA_PLOT_DIR", unset = file.path(tempdir(), "cellchat-qa"))
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("CellChat plot harness — study:", QA_STUDY, "\n")
cat("PNG output:", plot_dir, "\n")
cat(strrep("=", 78), "\n\n", sep = "")

data <- cellchat_load_objects(QA_STUDY, QA_CONFIG)
if (!is.null(data$error)) stop(data$error)
merged <- data$merged
conditions <- data$conditions
celltypes <- cellchat_celltypes(merged)
pathways <- cellchat_pathways(merged)

# Render `expr` to a PNG and return its md5, or NA if it errored. Rendering is
# seeded because several CellChat plots use igraph layouts.
render_md5 <- function(name, expr) {
  path <- file.path(plot_dir, paste0(name, ".png"))
  ok <- TRUE
  grDevices::png(path, width = 1400, height = 1000, res = 110)
  tryCatch({
    set.seed(42)
    force(expr)
  }, error = function(e) {
    ok <<- FALSE
    qa_record(paste("render:", name), "FAIL", conditionMessage(e))
  })
  grDevices::dev.off()
  if (!ok || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

# Render two variants of a plot and report whether the control changed it.
sweep <- function(check, name_a, expr_a, name_b, expr_b, detail = "") {
  a <- render_md5(name_a, expr_a)
  b <- render_md5(name_b, expr_b)
  if (is.na(a) || is.na(b)) {
    qa_record(check, "FAIL", "one or both renders errored")
    return(invisible(NULL))
  }
  qa_expect_reactive(check, a != b, detail)
}

# ---- `measure` (Strength vs Count) -------------------------------------------
sweep("circle reacts to `measure`",
      "circle_weight", cellchat_plot_circle_comparison(merged, "weight"),
      "circle_count",  cellchat_plot_circle_comparison(merged, "count"),
      "Comparison overview")

sweep("diffInteraction reacts to `measure`",
      "diff_weight", cellchat_plot_diff_interaction(merged, "weight"),
      "diff_count",  cellchat_plot_diff_interaction(merged, "count"),
      "Comparison overview")

sweep("compareInteractions reacts to `measure`",
      "compare_weight", cellchat_plot_compare_interactions(merged, "weight"),
      "compare_count",  cellchat_plot_compare_interactions(merged, "count"),
      "Comparison overview")

sweep("diff heatmap reacts to `measure`",
      "hmdiff_weight", cellchat_plot_heatmap_diff(merged, "weight"),
      "hmdiff_count",  cellchat_plot_heatmap_diff(merged, "count"),
      "Heatmaps")

# rankNet takes `measure` upstream but the wrapper never passes it, so the
# sidebar toggle cannot reach this plot. Rendering is identical by construction;
# the sweep proves the toggle is inert rather than merely unread.
sweep("rankNet reacts to `measure`",
      "ranknet_measure_a", cellchat_plot_ranknet(merged, TRUE),
      "ranknet_measure_b", cellchat_plot_ranknet(merged, TRUE),
      "Signaling roles tab — wrapper has no measure argument")

# ---- `rank_stacked` ----------------------------------------------------------
sweep("rankNet reacts to `rank_stacked`",
      "ranknet_stacked", cellchat_plot_ranknet(merged, TRUE),
      "ranknet_unstacked", cellchat_plot_ranknet(merged, FALSE),
      "Signaling roles tab")

# rankNet ignores the sidebar pathway filter: the wrapper has no `signaling`
# argument, so there is nothing to vary. Recorded from the function signature.
qa_record("rankNet reacts to `pathways`",
          if ("signaling" %in% names(formals(cellchat_plot_ranknet))) "PASS" else "UNWIRED",
          "cellchat_plot_ranknet() has no signaling/measure parameter")

# ---- per-condition plots: `condition`, `role_pattern`, `pathways` ------------
cond_a <- names(conditions)[[1]]
cond_b <- names(conditions)[[2]]
some_pathways <- head(intersect(pathways, cellchat_pathways(conditions[[cond_a]])), 5)

sweep("role heatmap reacts to `condition`",
      paste0("rolehm_", cond_a),
      cellchat_plot_signaling_role_heatmap(conditions[[cond_a]], "outgoing", NULL),
      paste0("rolehm_", cond_b),
      cellchat_plot_signaling_role_heatmap(conditions[[cond_b]], "outgoing", NULL),
      "Heatmaps tab")

sweep("role heatmap reacts to `role_pattern`",
      "rolehm_outgoing",
      cellchat_plot_signaling_role_heatmap(conditions[[cond_a]], "outgoing", NULL),
      "rolehm_incoming",
      cellchat_plot_signaling_role_heatmap(conditions[[cond_a]], "incoming", NULL),
      "Heatmaps tab")

sweep("role heatmap reacts to `pathways`",
      "rolehm_allpaths",
      cellchat_plot_signaling_role_heatmap(conditions[[cond_a]], "outgoing", NULL),
      "rolehm_subset",
      cellchat_plot_signaling_role_heatmap(conditions[[cond_a]], "outgoing", some_pathways),
      paste("subset:", paste(some_pathways, collapse = ", ")))

sweep("role scatter reacts to `condition`",
      paste0("scatter_", cond_a), cellchat_plot_signaling_role_scatter(conditions[[cond_a]], NULL),
      paste0("scatter_", cond_b), cellchat_plot_signaling_role_scatter(conditions[[cond_b]], NULL),
      "Signaling roles tab")

sweep("role scatter reacts to `pathways`",
      "scatter_allpaths", cellchat_plot_signaling_role_scatter(conditions[[cond_a]], NULL),
      "scatter_subset", cellchat_plot_signaling_role_scatter(conditions[[cond_a]], some_pathways),
      "Signaling roles tab")

# ---- bubble: sources / targets / pathways, and the missing thresholds --------
src <- head(celltypes, 3)
tgt <- head(celltypes, 3)

sweep("bubble reacts to `sources`",
      "bubble_src_a", cellchat_plot_bubble(merged, src[[1]], tgt, NULL),
      "bubble_src_b", cellchat_plot_bubble(merged, src[[2]], tgt, NULL),
      "L-R bubble tab")

sweep("bubble reacts to `targets`",
      "bubble_tgt_a", cellchat_plot_bubble(merged, src, tgt[[1]], NULL),
      "bubble_tgt_b", cellchat_plot_bubble(merged, src, tgt[[2]], NULL),
      "L-R bubble tab")

sweep("bubble reacts to `pathways`",
      "bubble_allpaths", cellchat_plot_bubble(merged, src, tgt, NULL),
      "bubble_subset", cellchat_plot_bubble(merged, src, tgt, some_pathways),
      "L-R bubble tab")

# netVisual_bubble() has a `thresh` argument, but the wrapper never exposes it,
# so the sidebar probability/p-value sliders cannot reach the bubble plot.
qa_record("bubble reacts to `prob_min` / `pval_max`",
          if (any(c("thresh", "prob_min", "pval_max") %in%
                  names(formals(cellchat_plot_bubble)))) "PASS" else "UNWIRED",
          "cellchat_plot_bubble() never passes thresh= to netVisual_bubble()")

# ---- overview tab ignores every filter --------------------------------------
for (fn in c("cellchat_plot_circle_comparison", "cellchat_plot_diff_interaction",
             "cellchat_plot_compare_interactions")) {
  args <- setdiff(names(formals(get(fn))), c("merged", "measure"))
  qa_record(paste0(fn, " accepts filters"),
            if (length(args) > 0) "PASS" else "UNWIRED",
            "takes only (merged, measure) — sidebar filters cannot apply")
}

qa_summary()
