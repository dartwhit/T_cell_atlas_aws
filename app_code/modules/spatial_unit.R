library(shiny)
library(bslib)
library(Seurat)

# ---------- MODULE UI ----------
spatial_UI <- function(id) {
  ns <- NS(id)
  plotOutput(ns("umap"), height = 500)
  layout_sidebar(
    sidebar = sidebar(
      id = ns("side"),
      title = "Spatial data options",
      open = TRUE,
      # div(
      #   class = "d-flex gap-2",
      #   # sidebar_toggle(ns("side"), label = "Toggle sidebar"),
      #   actionButton(ns("reset"),"Reset")
      # ),
      # tags$hr(),
      selectInput(ns("group_by"), 
                  "Group by (metadata)", 
                  choices = c("Seurat cluster" = "cluster",
                  "Key regions" =  "Key_Regions"),
                selected = "cluster"),
      selectizeInput(ns("feature"), "Feature (Gene/pathway)",
      choices = NULL, options = list(placeholder = "Type to search…"))
    ),# End of sidebar
    # ----------- Right side --------------
    layout_columns(
      col_widths = c(6,6),
      card(
        card_header("UMAP (Idents)"),
        plotOutput(ns("umap"))
      ),
      card(
        card_header("Feature plot"),
        plotOutput(ns("featureplot"))
      )
    )
  )
}

# ---------- MODULE SERVER ----------
# You can pass a ready Seurat object via `spat_obj`,
# or pass `rds_path` to load on the fly (useful for testing).
spatial_server <- function(id, spat_obj = NULL, rds_path = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Load object once if not provided
    obj <- if (!is.null(spat_obj)) {
      spat_obj
    } else {
      req(rds_path)
      readRDS(rds_path)
    }

    updateSelectizeInput(session, "feature", choices = rownames(obj[["SCT"]]), server = TRUE)


    

    redn <- reactive({
      if ("umap" %in% names(obj@reductions)) "umap" else
              if ("UMAP" %in% names(obj@reductions)) "UMAP" else NULL
    })

    output$umap <- renderPlot({
      req(redn())
      DimPlot(obj, reduction = redn(), group.by = input$group_by)
    })

    output$featureplot <- renderPlot({
      DefaultAssay(obj) <- "SCT"
      req(input$feature)
            FeaturePlot(
        obj, features = input$feature,
        reduction = redn()
      )
    })
    
  })
}
