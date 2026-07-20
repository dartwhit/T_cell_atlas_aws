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

sweep("rankNet reacts to `measure`",
      "ranknet_weight", cellchat_plot_ranknet(merged, TRUE, measure = "weight"),
      "ranknet_count",  cellchat_plot_ranknet(merged, TRUE, measure = "count"),
      "Signaling roles tab")

# ---- `rank_stacked` ----------------------------------------------------------
sweep("rankNet reacts to `rank_stacked`",
      "ranknet_stacked", cellchat_plot_ranknet(merged, TRUE),
      "ranknet_unstacked", cellchat_plot_ranknet(merged, FALSE),
      "Signaling roles tab")

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

sweep("bubble reacts to `pval_max`",
      "bubble_pval_1", cellchat_plot_bubble(merged, src, tgt, NULL, pval_max = 1),
      "bubble_pval_01", cellchat_plot_bubble(merged, src, tgt, NULL, pval_max = 0.01),
      "L-R bubble tab — pval_max feeds netVisual_bubble(thresh=)")

# `prob_min` has no counterpart in netVisual_bubble(): CellChat exposes only a
# p-value threshold there. It stays a table-only filter, and the sidebar says so.
qa_record("bubble reacts to `prob_min`", "UNWIRED",
          "by design — netVisual_bubble() has no probability cutoff; labelled table-only")

# ---- rankNet now takes the filters ------------------------------------------
sweep("rankNet reacts to `pathways`",
      "ranknet_allpaths", cellchat_plot_ranknet(merged, TRUE),
      "ranknet_subset", cellchat_plot_ranknet(merged, TRUE, signaling = some_pathways),
      "Signaling roles tab")

sweep("rankNet reacts to `sources`",
      "ranknet_src_all", cellchat_plot_ranknet(merged, TRUE),
      "ranknet_src_sub", cellchat_plot_ranknet(merged, TRUE, sources = src),
      "Signaling roles tab")

# rankNet deliberately takes no p-value argument. Guard the reason: CellChat's
# own thresh is inert here, so if a future version starts honouring it this
# check fails and the slider can be wired up for real.
ranknet_thresh_data <- function(th) {
  CellChat::rankNet(merged, mode = "comparison", comparison = 1:2, measure = "weight",
                    thresh = th, stacked = TRUE, do.stat = TRUE,
                    return.data = TRUE)$signaling.contribution
}
qa_record("rankNet exposes no `pval_max` (CellChat thresh is inert)",
          if (identical(ranknet_thresh_data(1), ranknet_thresh_data(0.001))) "PASS" else "FAIL",
          "rankNet(thresh=1) == rankNet(thresh=0.001); slider correctly omitted")

# ---- differential circle plot now takes source/target ------------------------
sweep("diffInteraction reacts to `sources`",
      "diff_src_all", cellchat_plot_diff_interaction(merged, "weight"),
      "diff_src_sub", cellchat_plot_diff_interaction(merged, "weight", sources = src),
      "Comparison overview")

sweep("diffInteraction reacts to `targets`",
      "diff_tgt_all", cellchat_plot_diff_interaction(merged, "weight"),
      "diff_tgt_sub", cellchat_plot_diff_interaction(merged, "weight", targets = tgt),
      "Comparison overview")

# The two aggregate overview plots are deliberately unfiltered: a circle plot of
# every cell type and a total-interaction bar chart are the tab's whole point.
# The sidebar carries a note saying so.
for (fn in c("cellchat_plot_circle_comparison", "cellchat_plot_compare_interactions")) {
  qa_record(paste0(fn, " accepts filters"), "UNWIRED",
            "by design — aggregate overview; sidebar states filters do not apply")
}

qa_summary()
