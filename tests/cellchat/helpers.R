# Shared bootstrap + result collection for the CellChat QA harness.
#
# Both harness scripts are run from app_code/ so that data_path() resolves
# against the same relative "data/" prefix the app uses:
#
#   cd app_code && Rscript ../tests/cellchat/test_wiring.R

if (!file.exists("setup.R")) {
  stop("Run this harness from app_code/: cd app_code && Rscript ../tests/cellchat/<script>.R")
}

suppressPackageStartupMessages({
  library(shiny)
  source("setup.R")
  source("modules/cellchat_helpers.R")
  source("modules/cellchat_explorer_module.R")
})

QA_STUDY <- "tabib"
QA_CONFIG <- dataset_files[[QA_STUDY]]$cellchat

# ---- result collection -------------------------------------------------------
#
# Three outcomes, deliberately distinguished:
#   PASS     - the control drives the output, as the UI implies
#   UNWIRED  - the control demonstrably does NOT reach the output (a finding,
#              not a harness failure; does not affect the exit code)
#   FAIL     - the harness itself could not complete the check

qa_results <- new.env(parent = emptyenv())
qa_results$rows <- list()

qa_record <- function(check, status, detail = "") {
  qa_results$rows[[length(qa_results$rows) + 1]] <-
    list(check = check, status = status, detail = detail)
  cat(sprintf("%-8s %-52s %s\n", status, check, detail))
  invisible(NULL)
}

# Expect a control to change an output. `changed` is the observed boolean.
qa_expect_reactive <- function(check, changed, detail = "") {
  qa_record(check, if (isTRUE(changed)) "PASS" else "UNWIRED", detail)
}

qa_try <- function(check, expr) {
  tryCatch(expr, error = function(e) {
    qa_record(check, "FAIL", conditionMessage(e))
    NULL
  })
}

qa_summary <- function() {
  statuses <- vapply(qa_results$rows, `[[`, character(1), "status")
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(sprintf("PASS: %d   UNWIRED: %d   FAIL: %d\n",
              sum(statuses == "PASS"), sum(statuses == "UNWIRED"),
              sum(statuses == "FAIL")))
  unwired <- qa_results$rows[statuses == "UNWIRED"]
  if (length(unwired) > 0) {
    cat("\nControls that do not reach their apparent output:\n")
    for (row in unwired) cat("  - ", row$check, "\n", sep = "")
  }
  # Only harness breakage is an error; wiring gaps are the report's content.
  if (any(statuses == "FAIL")) quit(status = 1)
  invisible(NULL)
}
