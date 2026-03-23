# ---- Auth toggle ----
# Set to FALSE to skip login (for local testing)
# Set to TRUE for deployed / production use
AUTH_ENABLED <- FALSE

library(shiny)

# Wrap source() calls so any init error surfaces in the browser
# rather than producing a blank 'No UI defined' page.
init_error <- NULL

tryCatch({
  source("ui.R")
  source("server.R")
}, error = function(e) {
  msg <- conditionMessage(e)
  cat("FATAL: app initialisation error:", msg, "\n", file = stderr())
  init_error <<- msg
})

if (!is.null(init_error)) {
  # Minimal error UI so the problem is visible in the browser
  ui <- fluidPage(
    tags$h2("Application failed to start", style = "color:red"),
    tags$p("Check the Shiny Server logs for details."),
    tags$pre(init_error)
  )
  server <- function(input, output, session) {}
}

shinyApp(ui = ui, server = server)