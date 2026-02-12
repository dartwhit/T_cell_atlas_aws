library(shiny)

source("ui.R")
source("server.R")

# Wrap UI with shinymanager authentication (unless disabled via environment variable)
disable_auth <- Sys.getenv("DISABLE_AUTH", "false") == "true"
if (disable_auth) {
  cat("⚠️  Authentication DISABLED via DISABLE_AUTH environment variable\n", file = stderr())
  # ui is already defined from ui.R, no wrapping needed
} else {
  cat("✅ Authentication ENABLED - wrapping UI with secure_app()\n", file = stderr())
  library(shinymanager)
  tryCatch({
    ui <- secure_app(ui, enable_admin = TRUE)
    cat("✓ UI wrapped with authentication successfully\n", file = stderr())
  }, error = function(e) {
    cat("❌ ERROR wrapping UI with secure_app():", conditionMessage(e), "\n", file = stderr())
    cat("⚠️  Continuing with unwrapped UI\n", file = stderr())
    # ui remains as the original page_navbar() object
  })
}

# Verify ui is defined
if (!exists("ui")) {
  stop("ERROR: ui object was not found after sourcing ui.R!")
}
cat("✓ UI object ready for shinyApp\n", file = stderr())

shinyApp(ui = ui, server = server)