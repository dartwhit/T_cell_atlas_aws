library(bslib)
library(shiny)
source("spatial_unit.R")

# ---- Test App ----
test_ui <- bslib::page_fluid(
  theme = bs_theme(bootswatch = "flatly"),
  h3("Spatial module test"),
  layout_sidebar(
    sidebar = sidebar(
      title = "Spatial data options",
      selectInput(NS("sp1", "group_by"),
                  "Group by (metadata)",
                  choices = c("Seurat cluster" = "cluster",
                  "Key regions" =  "Key_Regions"),
                selected = "cluster"),
      selectizeInput(NS("sp1", "feature"), "Feature (Gene/pathway)",
      choices = NULL, options = list(placeholder = "Type to search…")),
      checkboxGroupInput(NS("sp1", "samples"),
                    "Select samples to view"
                  )
    ),
    spatial_UI("sp1")
  )
)

test_server <- function(input, output, session) {
  spatial_server(
    id = "sp1",
    # Option A: pass a path for quick testing
    rds_path = reactive("../data/2025-07-07_MaSSc_Visium_PRECAST_SingleCellPredicted_RegionsNamed_CARD.rds")
  )
}

spat_unit_app <- shinyApp(test_ui, test_server)