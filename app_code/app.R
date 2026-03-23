library(shiny)

# ---- Auth toggle ----
# Set to FALSE to skip login (for local testing)
# Set to TRUE for deployed / production use
AUTH_ENABLED <- FALSE

source("ui.R")
source("server.R")
shinyApp(ui = ui, server = server)