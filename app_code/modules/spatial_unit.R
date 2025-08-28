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
      selectizeInput(ns("feature"), "Feature (Gene/pathway)",multiple = TRUE,
      choices = NULL, options = list(placeholder = "Type to search…")),
      # Select sampels to show
      checkboxGroupInput(ns("samples"),
                    "Select samples to view"
                  )
    ),# End of sidebar
    # ----------- Right side --------------
    layout_columns(
      col_widths = c(6,6),
      card(
        card_header("UMAP (Idents)"),
        plotOutput(ns("umap"),height=380)
      ),
      card(
        card_header("Feature plot"),
        plotOutput(ns("featureplot"),height=380)
      )
    ),
      # Dynamic grid of spatial plots (one per checked sample)
    card(
      card_header("Spatial Dimplots"),
      uiOutput(ns("spatial_grid"))
    ),
     # Dynamic UI for spatial featureplots (one per checked sample)
     card(
      card_header("Spatial Featureplots"),
      uiOutput(ns("spatial_feature_grid"))
  )
)}

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

    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

    updateSelectizeInput(session, "feature", choices = rownames(obj[["SCT"]]), server = TRUE)
    updateCheckboxGroupInput(session, "samples", choices = names(obj@images))

    

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
    
    # ---- Dynamic Spatial plots (one per checked sample) ----
    # 1) Build UI placeholders for each selected sample
    output$spatial_grid <- renderUI({
      req(length(input$samples) > 0)
      cards <- lapply(input$samples, function(s) {
        out_id    <- paste0("sp_", safe_id(s))
        slider_id <- paste0("pt_", safe_id(s))
        card(
          card_header(paste("Sample:", s)),
          div(class = "px-3 py-2",
              sliderInput(ns(slider_id), "Point size (spots)",
                          min = 0.1, max = 4, value = 1.6, step = 0.1)
          ),
          plotOutput(ns(out_id), height = 380)
        )
      })
      do.call(layout_columns, c(list(col_widths = rep(4, length(cards))), cards))
    })

     # 2) Attach renderers for each selected sample
    observe({
      # Re-wire outputs whenever the selection changes
      req(length(input$samples) > 0)
      lapply(input$samples, function(s) {
        out_id <- paste0("sp_", safe_id(s))
        slider_id <- paste0("pt_", safe_id(s))
        # local() to capture `s` correctly inside the loop
        local({
          s_local <- s
          out_local <- out_id
          slider_local <- slider_id
          output[[out_local]] <- renderPlot({
            size_val <- input[[slider_local]]
            grp <- if (!is.null(input$group_by)) input$group_by else NULL
            SpatialDimPlot(
              obj,
              images = s_local,         # <-- one plot per sample
              group.by = grp,
              pt.size.factor = size_val
            )
          })
        })
      })
    })


    # ---- Dynamic Spatial Feature plots (one per checked sample) ----
    # 1) Build UI placeholders for each selected sample
    output$spatial_feature_grid <- renderUI({
      req(length(input$samples) > 0, input$feature)
      cards <- lapply(input$samples, function(s) {
        out_id    <- paste0("spfeat_", safe_id(s))
        slider_id <- paste0("ptfeat_", safe_id(s))
        card(
          card_header(paste("Sample:", s)),
          div(class = "px-3 py-2",
              sliderInput(ns(slider_id), "Point size (spots)",
                          min = 0.1, max = 4, value = 1.6, step = 0.1)
          ),
          plotOutput(ns(out_id), height = 380)
        )
      })
      do.call(layout_columns, c(list(col_widths = rep(4, length(cards))), cards))
    })

     # 2) Attach renderers for each selected sample
    observe({
      # Re-wire outputs whenever the selection changes
      req(length(input$samples) > 0, input$feature)
      lapply(input$samples, function(s) {
        out_id <- paste0("spfeat_", safe_id(s))
        slider_id <- paste0("ptfeat_", safe_id(s))
        # local() to capture `s` correctly inside the loop
        local({
          s_local <- s
          out_local <- out_id
          slider_local <- slider_id
          output[[out_local]] <- renderPlot({
            size_val <- input[[slider_local]]
            DefaultAssay(obj) <- "SCT"
            SpatialFeaturePlot(
              obj,
              images = s_local,         # <-- one plot per sample
              features = input$feature,
              pt.size.factor = size_val
            )
          })
        })
      })
  })
}
  )}
