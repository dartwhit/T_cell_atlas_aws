library(bslib)
library(shiny)
# ---- Test App ----
test_ui <- bslib::page_fluid(
  theme = bs_theme(bootswatch = "flatly"),
  h3("Spatial module test"),
  spatial_UI("sp1")
)

test_server <- function(input, output, session) {
  spatial_server(
    id = "sp1",
    # Option A: pass a path for quick testing
    rds_path = "../data/2025-07-07_MaSSc_Visium_PRECAST_SingleCellPredicted_RegionsNamed_CARD.rds"

    # Option B (faster at runtime): load once globally and pass `spat_obj = ...`
    # spat_obj = my_spatial  # if you did: my_spatial <- readRDS("...rds")
  )
}

spat_unit_app <- shinyApp(test_ui, test_server)
