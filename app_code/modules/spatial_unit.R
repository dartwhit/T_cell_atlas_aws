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
        cat("Loading spatial object from spat_obj parameter\n", file = stderr())
        # Update object to ensure compatibility
        tryCatch({
          updated_obj <- UpdateSeuratObject(spat_obj)
          cat("✓ Seurat object updated\n", file = stderr())
          updated_obj
        }, error = function(e) {
          cat("Warning: Could not update Seurat object:", e$message, "\n", file = stderr())
          spat_obj
        })
      } else {
        req(rds_path())
        cat("Loading spatial object from path:", rds_path(), "\n", file = stderr())
        tryCatch({
          loaded_obj <- readRDS(rds_path())
          cat("✓ Spatial object loaded successfully\n", file = stderr())
          cat("  - Images:", paste(names(loaded_obj@images), collapse = ", "), "\n", file = stderr())
          cat("  - Assays:", paste(names(loaded_obj@assays), collapse = ", "), "\n", file = stderr())
          
          # Update object to ensure compatibility with current Seurat version
          cat("Updating Seurat object...\n", file = stderr())
          updated_obj <- UpdateSeuratObject(loaded_obj)
          cat("✓ Seurat object updated for compatibility\n", file = stderr())
          updated_obj
        }, error = function(e) {
          cat("✗ Error loading/updating spatial object:", e$message, "\n", file = stderr())
          showNotification(paste("Error loading spatial data:", e$message), type = "error", duration = NULL)
          NULL
        })
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
      # Determine which assay to use for feature selection
      assay_to_use <- if ("SCT" %in% names(obj()@assays)) {
        "SCT"
      } else if ("Spatial" %in% names(obj()@assays)) {
        "Spatial"
      } else {
        "RNA"
      }
      
      tryCatch({
        updateSelectizeInput(session, "feature", choices = rownames(obj()[[assay_to_use]]), server = TRUE)
      }, error = function(e) {
        cat("Error updating feature choices:", e$message, "\n", file = stderr())
      })
      
      updateCheckboxGroupInput(session, "samples", choices = names(obj()@images))
    })

    redn <- reactive({
      req(obj())
      if ("umap" %in% names(obj()@reductions)) "umap" else
              if ("UMAP" %in% names(obj()@reductions)) "UMAP" else NULL
    })

    output$umap <- renderPlot({
      req(redn(), obj())
      tryCatch({
        DimPlot(obj(), reduction = redn(), group.by = input$group_by)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error generating UMAP:\n", e$message), cex = 1.2, col = "red")
      })
    })

    output$featureplot <- renderPlot({
      req(obj(), input$feature)
      tryCatch({
        plot_obj <- obj()
        # Determine which assay to use
        if ("SCT" %in% names(plot_obj@assays)) {
          DefaultAssay(plot_obj) <- "SCT"
        } else if ("Spatial" %in% names(plot_obj@assays)) {
          DefaultAssay(plot_obj) <- "Spatial"
        } else {
          DefaultAssay(plot_obj) <- "RNA"
        }
        FeaturePlot(
          plot_obj, features = input$feature,
          reduction = redn()
        )
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error generating feature plot:\n", e$message), cex = 1.2, col = "red")
      })
    })
    
    # ---- Sample-centric dynamic grid ----
    output$sample_centric_grid <- renderUI({
      req(length(input$samples) > 0)
      
      cards <- lapply(input$samples, function(s) {
        safe_s <- safe_id(s)
        slider_id <- paste0("pt_", safe_s)
        dimplot_id <- paste0("sp_dim_", safe_s)
        featureplot_id <- paste0("sp_feat_", safe_s)
        dld_dim_id <- paste0("dld_dim_", safe_s)
        dld_feat_id <- paste0("dld_feat_", safe_s)
        
        card(
          card_header(paste("Sample:", s)),
          div(class = "px-3 py-2",
              sliderInput(ns(slider_id), "Relative point size",
                          min = 0.2, max = 1.5, value = 1, step = 0.05)
          ),
          layout_columns(
            col_widths = c(6, 6),
            plotOutput(ns(dimplot_id), height = 380),
            conditionalPanel(
              condition = paste0("input['", ns("feature"), "'] && input['", ns("feature"), "'].length > 0"),
              plotOutput(ns(featureplot_id), height = 380)
            )
          ),
          card_footer(
            downloadButton(ns(dld_dim_id), "Download Clusters"),
            conditionalPanel(
              condition = paste0("input['", ns("feature"), "'] && input['", ns("feature"), "'].length > 0"),
              downloadButton(ns(dld_feat_id), "Download Features")
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
        dld_dim_id <- paste0("dld_dim_", safe_s)
        dld_feat_id <- paste0("dld_feat_", safe_s)
        
        local({
          s_local <- s
          
          # Reactive for DimPlot
          dim_plot_reactive <- reactive({
            tryCatch({
              size_val_multiplier <- input[[slider_id]]
              if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
              default_size <- default_pt_sizes()[s_local]
              final_size <- default_size * size_val_multiplier
              grp <- if (!is.null(input$group_by)) input$group_by else NULL
              SpatialDimPlot(obj(), images = s_local, group.by = grp, pt.size.factor = final_size)
            }, error = function(e) {
              cat("Error in SpatialDimPlot for sample", s_local, ":", e$message, "\n")
              plot.new()
              text(0.5, 0.5, paste("Error:\n", e$message), cex = 1, col = "red")
            })
          })
          
          # Reactive for FeaturePlot
          feat_plot_reactive <- reactive({
            req(input$feature)
            tryCatch({
              size_val_multiplier <- input[[slider_id]]
              if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
              default_size <- default_pt_sizes()[s_local]
              final_size <- default_size * size_val_multiplier
              plot_obj <- obj()
              # Determine which assay to use
              if ("SCT" %in% names(plot_obj@assays)) {
                DefaultAssay(plot_obj) <- "SCT"
              } else if ("Spatial" %in% names(plot_obj@assays)) {
                DefaultAssay(plot_obj) <- "Spatial"
              } else {
                DefaultAssay(plot_obj) <- "RNA"
              }
              SpatialFeaturePlot(plot_obj, images = s_local, features = input$feature, pt.size.factor = final_size)
            }, error = function(e) {
              cat("Error in SpatialFeaturePlot for sample", s_local, ":", e$message, "\n", file = stderr())
              plot.new()
              text(0.5, 0.5, paste("Error:\n", e$message), cex = 1, col = "red")
            })
          })

          # Render plots
          output[[dimplot_id]] <- renderPlot(dim_plot_reactive())
          output[[featureplot_id]] <- renderPlot(feat_plot_reactive())
          
          # Download handlers
          output[[dld_dim_id]] <- downloadHandler(
            filename = function() { paste0(s_local, "_clusters.png") },
            content = function(file) {
              png(file, width = 7, height = 7, units = "in", res = 150)
              print(dim_plot_reactive())
              dev.off()
            }
          )
          
          output[[dld_feat_id]] <- downloadHandler(
            filename = function() { paste0(s_local, "_features.png") },
            content = function(file) {
              png(file, width = 7, height = 7, units = "in", res = 150)
              print(feat_plot_reactive())
              dev.off()
            }
          )
        })
      })
    })
  })
}
