# CellChat explorer — QA harness

Functional QA of the CellChat tab
([cellchat_explorer_module.R](../../app_code/modules/cellchat_explorer_module.R),
[cellchat_helpers.R](../../app_code/modules/cellchat_helpers.R)), run against the
real `tabib` objects (25 cell types, 122 pathways, HC + SSc, 32,860 L-R rows).

The first pass found nine controls that did not reach the output they appeared
to drive; those are now fixed. The harness stays here so the wiring can be
re-checked whenever the module or the pinned CellChat version changes.

**Current state: 39 PASS, 3 UNWIRED (all documented below as deliberate), 0 FAIL.**

## Running

Both scripts must be run from `app_code/` so `data_path()` resolves:

```bash
cd app_code
Rscript ../tests/cellchat/test_wiring.R
Rscript ../tests/cellchat/test_plots.R
```

Needs `CellChat`, `ComplexHeatmap`, `shiny`, `Seurat` — all already in the app
image. No Docker and no browser driver: wiring uses `shiny::testServer()`, plots
use `grDevices::png()`. Set `QA_PLOT_DIR` to keep the rendered PNGs.

Each line prints `PASS`, `UNWIRED`, or `FAIL`:

- **PASS** — the control drives the output, as the UI implies.
- **UNWIRED** — the control does *not* reach the output. Only three remain, each
  a deliberate, documented limitation; does not affect the exit code.
- **FAIL** — the harness itself could not complete the check. Exit code 1.

`test_plots.R` decides `PASS`/`UNWIRED` by rendering two variants and comparing
md5 hashes: a control that really drives a plot produces a different image; one
that never reaches the CellChat call produces a byte-identical image. This is
what caught the original bugs, and what keeps them from coming back.

## What each control drives

| Control | Drives |
|---|---|
| `study` | everything |
| `measure` (Strength/Count) | circle, diffInteraction, compareInteractions, diff heatmap, **rankNet** |
| `condition` | role heatmap, role scatter |
| `table_condition` | communication table, CSV download |
| `role_pattern` | role heatmap |
| `rank_stacked` | rankNet |
| `sources` / `targets` | **diffInteraction**, bubble, **rankNet**, table |
| `pathways` | role heatmap, role scatter, bubble, **rankNet**, table |
| `ligand` / `receptor` | table |
| `prob_min` | table |
| `pval_max` | **bubble**, table |
| `level` | table, download filename |

Bold entries were wired up in this pass.

## Fixes applied

1. **`measure` reached rankNet.** `cellchat_plot_ranknet()` gained a `measure`
   parameter; the Strength/Count toggle had been silently inert on that tab.
2. **`pathways` and `sources`/`targets` reached rankNet** via `signaling=`,
   `sources.use=` and `targets.use=`.
3. **`pval_max` reached the L-R bubble** via `netVisual_bubble(thresh=)`. Added
   `cellchat_pval_threshold()`, which nudges the cutoff by 1e-9 because CellChat
   filters on `pval < thresh` while the table filters on `pval <= pval_max` —
   without it the two views disagree on interactions sitting exactly on the
   boundary.
4. **The communication table can be scoped to one condition.** It reads from the
   merged object, so it had always returned every condition stacked together
   (HC 14,594 + SSc 18,266 rows) with no way to narrow it. The sidebar
   `condition` control drives the per-condition *plots* and cannot take an "all"
   value without breaking them, so the table got its own `table_condition`
   selector in its tab, defaulting to "All conditions" — matching how `level`,
   `role_pattern` and `rank_stacked` already live beside the output they affect.
5. **`draw_plot()` no longer swallows `req()`.** Its `tryCatch(error=)` caught
   `shiny.silent.error` (which inherits from `error`), so an unset input made the
   Heatmaps and Signaling-roles tabs render `"CellChat plot could not be
   generated:"` with an empty message on first load. Silent errors are now
   re-raised; genuine errors are still drawn, and the harness checks both.
6. **`sources`/`targets` reached the differential circle plot**, which does
   support them. The sidebar now also states which outputs each block governs.
7. **`conditions[[1]] %||% NULL` replaced with a length check.** On an empty
   vector the subscript threw before `%||%` was ever reached.

## Remaining UNWIRED — deliberate

- **`prob_min` does not affect the L-R bubble.** `netVisual_bubble()` exposes
  only a p-value threshold, no probability cutoff. It stays a table-only filter
  and the sidebar says so.
- **`pval_max` is not offered for rankNet.** `rankNet()` accepts `thresh`, but it
  provably does nothing here — identical output from `thresh = 1` down to
  `0.001`. Wiring it would put a slider on screen that changes nothing, so the
  parameter is omitted. `test_plots.R` asserts the inertness, so if a future
  CellChat starts honouring `thresh` the check fails and the slider can be added.
- **The Comparison overview ignores the Filters block.** An aggregate circle plot
  of every cell type and a total-interaction bar chart are that tab's purpose.
  The differential plot beside them *does* now respect source/target.
