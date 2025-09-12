library(shiny)
library(bslib)
library(Seurat)

# ---------- MODULE UI ----------
spatial_UI <- function(id) {
  ns <- NS(id)
  tagList(
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
    # Sample-centric grid of plots
    conditionalPanel(
      condition = paste0("input['", ns("samples"), "'] && input['", ns("samples"), "'].length > 0"),
      uiOutput(ns("sample_centric_grid"))
    )
  )
}

# ---------- MODULE SERVER ----------
# You can pass a ready Seurat object via `spat_obj`,
# or pass `rds_path` to load on the fly (useful for testing).
spatial_server <- function(id, spat_obj = NULL, rds_path = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    obj <- reactive({
      if (!is.null(spat_obj)) {
        spat_obj
      } else {
        req(rds_path())
        readRDS(rds_path())
      }
    })

    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

    default_pt_sizes <- reactive({
      req(obj())
      sapply(names(obj()@images), function(s) {
        coords <- obj()@images[[s]]@coordinates
        if (nrow(coords) < 2) return(1.6) # Default if less than 2 spots
        dists <- as.matrix(dist(coords))
        diag(dists) <- Inf # Set distance to self to Inf
        min_dists <- apply(dists, 1, min)
        d_mean <- mean(min_dists)
        spot_radius <- obj()@images[[s]]@spot.radius
        scale_factor <- obj()@images[[s]]@scale.factors$spot
        
        # Set diameter to be a fraction of the mean distance.
        pt_size_factor <- (0.1 * d_mean) / (2 * spot_radius * scale_factor)
        pt_size_factor
      })
    })

    observe({
      req(obj())
      updateSelectizeInput(session, "feature", choices = rownames(obj()[["SCT"]]), server = TRUE)
      updateCheckboxGroupInput(session, "samples", choices = names(obj()@images))
    })

    redn <- reactive({
      req(obj())
      if ("umap" %in% names(obj()@reductions)) "umap" else
              if ("UMAP" %in% names(obj()@reductions)) "UMAP" else NULL
    })

    output$umap <- renderPlot({
      req(redn(), obj())
      DimPlot(obj(), reduction = redn(), group.by = input$group_by)
    })

    output$featureplot <- renderPlot({
      req(obj(), input$feature)
      plot_obj <- obj()
      DefaultAssay(plot_obj) <- "SCT"
      FeaturePlot(
        plot_obj, features = input$feature,
        reduction = redn()
      )
    })
    
    # ---- Sample-centric dynamic grid ----
    output$sample_centric_grid <- renderUI({
      req(length(input$samples) > 0)
      
      cards <- lapply(input$samples, function(s) {
        safe_s <- safe_id(s)
        slider_id <- paste0("pt_", safe_s)
        dimplot_id <- paste0("sp_dim_", safe_s)
        featureplot_id <- paste0("sp_feat_", safe_s)
        
        card(
          card_header(paste("Sample:", s)),
          div(class = "px-3 py-2",
              sliderInput(ns(slider_id), "Relative point size",
                          min = 0.2, max = 1.5, value = 1, step = 0.1)
          ),
          layout_columns(
            col_widths = c(6, 6),
            plotOutput(ns(dimplot_id), height = 380),
            conditionalPanel(
              condition = paste0("input['", ns("feature"), "'] && input['", ns("feature"), "'].length > 0"),
              plotOutput(ns(featureplot_id), height = 380)
            )
          )
        )
      })
      
      do.call(layout_columns, c(list(col_widths = rep(6, length(cards))), cards))
    })
    
    observe({
      req(length(input$samples) > 0, default_pt_sizes())
      
      lapply(input$samples, function(s) {
        safe_s <- safe_id(s)
        slider_id <- paste0("pt_", safe_s)
        dimplot_id <- paste0("sp_dim_", safe_s)
        featureplot_id <- paste0("sp_feat_", safe_s)
        
        local({
          s_local <- s
          
          # Render SpatialDimPlot
          output[[dimplot_id]] <- renderPlot({
            size_val_multiplier <- input[[slider_id]]
            if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
            
            default_size <- default_pt_sizes()[s_local]
            final_size <- default_size * size_val_multiplier
            
            grp <- if (!is.null(input$group_by)) input$group_by else NULL
            SpatialDimPlot(
              obj(),
              images = s_local,
              group.by = grp,
              pt.size.factor = final_size
            )
          })
          
          # Render SpatialFeaturePlot
          output[[featureplot_id]] <- renderPlot({
            req(input$feature)
            size_val_multiplier <- input[[slider_id]]
            if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
            
            default_size <- default_pt_sizes()[s_local]
            final_size <- default_size * size_val_multiplier
            
            plot_obj <- obj()
            DefaultAssay(plot_obj) <- "SCT"
            SpatialFeaturePlot(
              plot_obj,
              images = s_local,
              features = input$feature,
              pt.size.factor = final_size
            )
          })
        })
      })
    })
  })
}
