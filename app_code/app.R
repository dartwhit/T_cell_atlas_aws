# ---- Auth toggle ----
# Set to FALSE to skip login (for local testing)
# Set to TRUE for deployed / production use
AUTH_ENABLED <- FALSE

library(shiny)

source("ui.R")
source("server.R")
shinyApp(ui = ui, server = server)