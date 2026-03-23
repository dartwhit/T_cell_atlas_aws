# ---- Auth toggle ----
# Set to FALSE to skip login (for local testing)
# Set to TRUE for deployed / production use
AUTH_ENABLED <- FALSE

library(shiny)

cat("INFO: app.R starting – AUTH_ENABLED =", AUTH_ENABLED, "\n", file = stderr())

# Wrap source() calls so any init error surfaces in the browser
# rather than producing a blank 'No UI defined' page.
#
# IMPORTANT: use local = environment() so that `ui` and `server` are
# assigned into *this* environment (E_app), not into globalenv().
# Shiny's loadApp() captures the return value of shinyApp() or looks
# for `ui`/`server` in E_app; relying on globalenv can break when
# Shiny-server creates E_app with a non-standard parent chain.
init_error <- NULL

tryCatch({
  cat("INFO: sourcing ui.R\n", file = stderr())
  source("ui.R", local = environment())
  cat("INFO: sourcing server.R\n", file = stderr())
  source("server.R", local = environment())
  cat("INFO: ui.R and server.R loaded successfully\n", file = stderr())
}, error = function(e) {
  msg <- conditionMessage(e)
  call_str <- if (!is.null(conditionCall(e))) paste(deparse(conditionCall(e)), collapse = " ") else "(unknown call)"
  cat("FATAL: app initialisation error in", call_str, "—", msg, "\n", file = stderr())
  init_error <<- msg
})

# Safety net: if source() returned without error but somehow `ui` or
# `server` is still missing, surface a clear message instead of
# letting shinyApp() fail with an opaque "No UI defined" page.
if (is.null(init_error) && (!exists("ui") || !exists("server"))) {
  init_error <- paste(
    "ui or server was not defined after sourcing app files.",
    "This usually means a silent failure in ui.R or server.R.",
    "Check the Shiny Server logs for details."
  )
  cat("FATAL:", init_error, "\n", file = stderr())
}

if (!is.null(init_error)) {
  # Minimal error UI so the problem is visible in the browser
  ui <- fluidPage(
    tags$h2("Application failed to start", style = "color:red"),
    tags$p("Check the Shiny Server logs for details."),
    tags$pre(init_error)
  )
  server <- function(input, output, session) {}
}

cat("INFO: calling shinyApp()\n", file = stderr())
shinyApp(ui = ui, server = server)