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
      selectizeInput(NS("sp1", "feature"), "Feature (Gene/pathway)", multiple = TRUE,
      choices = NULL, options = list(placeholder = "Type to search…")),
      checkboxGroupInput(NS("sp1", "samples"),
                    "Select samples to view"
                  )
    ),
    spatial_UI("sp1")
  )
)

test_server <- function(input, output, session) {
  # Load the spatial object
  spat_obj <- reactive({
    readRDS("../data/spat_Li_RCTD_TMKMHdecon.rds")
  })
  
  # Populate feature choices (genes)
  observe({
    req(spat_obj())
    assay <- if ("SCT" %in% names(spat_obj()@assays)) "SCT" else DefaultAssay(spat_obj())
    genes <- rownames(spat_obj()[[assay]])
    updateSelectizeInput(session, NS("sp1", "feature"), 
                        choices = genes, 
                        server = TRUE)
  })
  
  # Populate sample choices
  observe({
    req(spat_obj())
    samples <- names(spat_obj()@images)
    updateCheckboxGroupInput(session, NS("sp1", "samples"),
                             choices = samples,
                             selected = samples[1])
  })
  
  spatial_server(
    id = "sp1",
    spat_obj = spat_obj
  )
}

spat_unit_app <- shinyApp(test_ui, test_server)