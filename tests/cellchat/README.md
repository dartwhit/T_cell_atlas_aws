# CellChat explorer — QA harness and findings

Functional QA of the CellChat tab
([cellchat_explorer_module.R](../../app_code/modules/cellchat_explorer_module.R),
[cellchat_helpers.R](../../app_code/modules/cellchat_helpers.R)), added on
`cellchat-integration` and not previously exercised end-to-end.

Run against the real `tabib` objects (25 cell types, 122 pathways, HC + SSc,
32,860 L-R rows). **No app code was changed** — this is a report.

## Running

Both scripts must be run from `app_code/` so `data_path()` resolves:

```bash
cd app_code
Rscript ../tests/cellchat/test_wiring.R
Rscript ../tests/cellchat/test_plots.R
```

Needs `CellChat`, `ComplexHeatmap`, `shiny`, `Seurat` — all already present
locally. No Docker, no browser driver: wiring uses `shiny::testServer()`, plots
use `grDevices::png()`. Set `QA_PLOT_DIR` to keep the rendered PNGs.

Each line prints `PASS`, `UNWIRED`, or `FAIL`:

- **PASS** — the control drives the output, as the UI implies.
- **UNWIRED** — the control demonstrably does *not* reach the output. A finding,
  not a harness failure; does not affect the exit code.
- **FAIL** — the harness itself could not complete the check. Exit code 1.

`test_plots.R` decides `PASS`/`UNWIRED` by rendering two variants and comparing
md5 hashes: a control that really drives a plot produces a different image; one
that never reaches the CellChat call produces a byte-identical image. Where the
wrapper has no parameter at all, the check reads the function signature instead.

## Result: 27 PASS, 9 UNWIRED, 0 FAIL

**Everything renders.** All 12 outputs produce real plots/tables against real
data — no crashes, no blank panels, no missing data files. The problems are all
about *which controls reach which outputs*, plus one cosmetic error-handling bug.

### What works

| Control | Drives | Verified by |
|---|---|---|
| `study` | all outputs | load test — HC + SSc objects, 25 cell types, 122 pathways |
| `measure` | circle, diffInteraction, compareInteractions, diff heatmap | hash sweep |
| `condition` | role heatmap, role scatter | hash sweep (HC vs SSc differ) |
| `role_pattern` | role heatmap | hash sweep (outgoing vs incoming differ) |
| `rank_stacked` | rankNet | hash sweep |
| `sources` / `targets` | bubble, table | hash sweep; 32860 → 393 rows |
| `pathways` | role heatmap, role scatter, bubble, table | hash sweep; 32860 → 21 rows |
| `ligand` / `receptor` | table | 32860 → 214 / 564 rows |
| `prob_min` / `pval_max` | table | 32860 → 16430 / 32668 rows |
| `level` | table, download filename | 32860 → 12970 rows |
| download button | non-empty CSV of the filtered table | write/read round-trip |

### What does not work

**1. `measure` (Strength/Count) is inert on the rankNet plot.**
[helpers.R:253](../../app_code/modules/cellchat_helpers.R#L253) —
`cellchat_plot_ranknet()` has no `measure` parameter, so `CellChat::rankNet()`
is always called with its default. Toggling Strength/Count on the "Signaling
roles / rankNet" tab changes nothing. Renders are byte-identical.

**2. `pathways` is inert on the rankNet plot.**
Same wrapper — no `signaling` parameter, though `rankNet()` accepts one. The
sidebar pathway filter narrows the heatmaps, scatter, bubble and table, but
rankNet keeps showing all 122 pathways.

**3. `prob_min` / `pval_max` are inert on the L-R bubble plot.**
[helpers.R:240](../../app_code/modules/cellchat_helpers.R#L240) —
`cellchat_plot_bubble()` never passes `thresh=` to `netVisual_bubble()`. The two
sliders sit under a heading that reads "Filters" but only affect the table.

**4. `condition` does not filter the communication table.**
[module.R:206-223](../../app_code/modules/cellchat_explorer_module.R#L206-L223) —
`communication_table()` runs on `data$merged` and `filtered_table()` never reads
`input$condition`. Measured: HC = 32,860 rows, SSc = 32,860 rows — identical.
The table always carries both conditions, distinguished only by its `dataset`
column, so each L-R pair appears twice and the control labelled "Condition for
drill-down" cannot drill down. The downloaded CSV has the same problem.

**5. `req()` inside `draw_plot()` renders a bogus error panel.**
[module.R:152-158](../../app_code/modules/cellchat_explorer_module.R#L152-L158) —
`draw_plot()` wraps the plot call in `tryCatch(error = ...)`. `req()` signals a
`shiny.silent.error`, which *inherits from* `error`, so the `req()` calls at
[module.R:184](../../app_code/modules/cellchat_explorer_module.R#L184) and
[module.R:201](../../app_code/modules/cellchat_explorer_module.R#L201) get
caught. Reproduced verbatim — the user sees:

```
CellChat plot could not be generated:
```

with an empty message, where the plot should simply stay blank. This fires on
the Heatmaps and Signaling-roles tabs on first load, before `condition` is
populated.

**6. The whole "Filters" block is inert on the Comparison overview tab.**
`cellchat_plot_circle_comparison()`, `cellchat_plot_diff_interaction()` and
`cellchat_plot_compare_interactions()` all take only `(merged, measure)`. Source,
target, pathway, ligand, receptor and both sliders cannot affect that tab at all.
Arguably by design for an overview, but the sidebar gives no such hint.

**7. Latent: `conditions[[1]] %||% NULL` is not a safe guard.**
[module.R:131](../../app_code/modules/cellchat_explorer_module.R#L131) — on an
empty vector, `conditions[[1]]` throws `subscript out of bounds` *before* `%||%`
is ever reached (confirmed). `cellchat_config_complete()` makes this unreachable
today, so it is a latent trap rather than a live bug.

### Summary

The sidebar's single "Filters" block implies the controls under it apply to
everything on the right. In practice they apply fully only to the Communication
table, partially to the bubble and heatmaps, and not at all to the Comparison
overview or rankNet. The cheapest structural fix is to pass the missing
arguments through the four wrappers (#1, #2, #3), scope the table by condition
(#4), and either move per-tab controls next to their tab or label the sidebar
block with what it governs.
