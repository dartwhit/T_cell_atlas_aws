library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)

# ---------- MODULE UI ----------
spatial_UI <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(
      col_widths = c(6,6),
      card(
        card_header("UMAP"),
        plotOutput(ns("umap"),height=380)
      ),
      card(
        card_header("Feature Plot"),
        plotOutput(ns("featureplot"),height=380)
      )
    ),
    # Sample-centric grid of plots
    uiOutput(ns("sample_centric_grid"))
  )
}

# ---------- SPATIAL SIDEBAR UI ----------
spatial_sidebar_UI <- function(id) {
  ns <- NS(id)
  sidebar(
    title = "Spatial Analysis Controls",
    position = "left",
    selectInput(ns("spatial_study_selector"), 
               "Select study", 
               choices = NULL),
    selectizeInput(ns("feature"), 
                  "Select feature to plot", 
                  choices = NULL, 
                  multiple = FALSE),
    selectInput(ns("group_by"), 
               "Color cells by", 
               choices = c("Default" = "", 
                         "Clusters" = "cluster",
                         "Original identity" = "orig.ident",
                         "Sample ID" = "SampleID",
                         "Predicted labels" = "predicted_CARD_labels"),
               selected = ""),
    checkboxGroupInput(ns("samples"), 
                      "Select samples to display", 
                      choices = NULL)
  )
}

# ---------- SPATIAL SIDEBAR SERVER ----------
spatial_sidebar_server <- function(id, studies_with_spatial) {
  moduleServer(id, function(input, output, session) {
    
    # Update the spatial study selector when studies change
    observe({
      req(studies_with_spatial())
      updateSelectInput(session, "spatial_study_selector", choices = studies_with_spatial())
    })
    
    return(
      list(
        spatial_study_selector = reactive(input$spatial_study_selector),
        feature = reactive(input$feature),
        group_by = reactive(input$group_by),
        samples = reactive(input$samples)
      )
    )
  })
}

# ---------- MODULE SERVER ----------
# You can pass a ready Seurat object via `spat_obj`,
# or pass `rds_path` to load on the fly (useful for testing).
spatial_server <- function(id, spat_obj = NULL, rds_path = NULL, sidebar_inputs = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    obj <- reactive({
      if (!is.null(spat_obj)) {
        if (is.reactive(spat_obj)) {
          spat_obj()
        } else {
          spat_obj
        }
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

    redn <- reactive({
      req(obj())
      if ("UMAP" %in% names(obj()@reductions)) "UMAP" else
              if ("umap" %in% names(obj()@reductions)) "umap" else NULL
    })

    output$umap <- renderPlot({
      req(redn(), obj())
      group_by_val <- if (!is.null(sidebar_inputs)) sidebar_inputs$group_by() else input$group_by
      
      # Handle empty or NULL group_by values
      if (is.null(group_by_val) || group_by_val == "") {
        p <- DimPlot(obj(), reduction = redn())
        title_text <- "UMAP"
      } else {
        # Check if the group_by column exists in metadata
        if (group_by_val %in% colnames(obj()@meta.data)) {
          p <- DimPlot(obj(), reduction = redn(), group.by = group_by_val)
          title_text <- paste("UMAP colored by", group_by_val)
        } else {
          # Fallback to no grouping if column doesn't exist
          p <- DimPlot(obj(), reduction = redn())
          title_text <- paste("UMAP (", group_by_val, "not found)")
        }
      }
      p + labs(title = title_text)
    })

    output$featureplot <- renderPlot({
      req(obj())
      feature_val <- if (!is.null(sidebar_inputs)) sidebar_inputs$feature() else input$feature
      req(feature_val)
      
      plot_obj <- obj()
      
      # Use SCT assay if available, otherwise use default assay
      if ("SCT" %in% names(plot_obj@assays)) {
        DefaultAssay(plot_obj) <- "SCT"
      }
      
      # Check if feature exists in the current assay
      current_assay <- DefaultAssay(plot_obj)
      available_features <- rownames(plot_obj[[current_assay]])
      
      if (!feature_val %in% available_features) {
        # Try to find the feature in other assays
        for (assay_name in names(plot_obj@assays)) {
          if (feature_val %in% rownames(plot_obj[[assay_name]])) {
            DefaultAssay(plot_obj) <- assay_name
            current_assay <- assay_name
            break
          }
        }
      }
      
      if (feature_val %in% rownames(plot_obj[[current_assay]])) {
        p <- FeaturePlot(
          plot_obj, features = feature_val,
          reduction = redn()
        )
        p + labs(title = paste("Feature:", feature_val, "on UMAP"))
      } else {
        # Create an empty plot if feature not found
        ggplot() + 
          labs(title = paste("Feature not found:", feature_val)) +
          theme_void()
      }
    })
    
    # ---- Sample-centric dynamic grid ----
    output$sample_centric_grid <- renderUI({
      samples_val <- if (!is.null(sidebar_inputs)) sidebar_inputs$samples() else input$samples
      req(length(samples_val) > 0)
      
      cards <- lapply(samples_val, function(s) {
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
              condition = paste0("input['", if (!is.null(sidebar_inputs)) "spatial_sidebar-feature" else ns("feature"), "'] && input['", if (!is.null(sidebar_inputs)) "spatial_sidebar-feature" else ns("feature"), "'].length > 0"),
              plotOutput(ns(featureplot_id), height = 380)
            )
          ),
          card_footer(
            downloadButton(ns(dld_dim_id), "Download Clusters"),
            conditionalPanel(
              condition = paste0("input['", if (!is.null(sidebar_inputs)) "spatial_sidebar-feature" else ns("feature"), "'] && input['", if (!is.null(sidebar_inputs)) "spatial_sidebar-feature" else ns("feature"), "'].length > 0"),
              downloadButton(ns(dld_feat_id), "Download Features")
            )
          )
        )
      })
      
      do.call(layout_columns, c(list(col_widths = rep(6, length(cards))), cards))
    })
    
    observe({
      samples_val <- if (!is.null(sidebar_inputs)) sidebar_inputs$samples() else input$samples
      req(length(samples_val) > 0, default_pt_sizes())
      
      lapply(samples_val, function(s) {
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
            size_val_multiplier <- input[[slider_id]]
            if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
            default_size <- default_pt_sizes()[s_local]
            final_size <- default_size * size_val_multiplier
            grp <- if (!is.null(sidebar_inputs)) sidebar_inputs$group_by() else input$group_by
            
            # Handle empty or NULL group_by values for spatial plots
            if (is.null(grp) || grp == "" || !grp %in% colnames(obj()@meta.data)) {
              p <- SpatialDimPlot(obj(), images = s_local, pt.size.factor = final_size)
              title_text <- paste0("Sample ", s_local)
            } else {
              p <- SpatialDimPlot(obj(), images = s_local, group.by = grp, pt.size.factor = final_size)
              title_text <- paste0("Sample ", s_local, " (", grp, ")")
            }
            
            p + labs(title = title_text) +
              theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
          })
          
          # Reactive for FeaturePlot
          feat_plot_reactive <- reactive({
            feature_val <- if (!is.null(sidebar_inputs)) sidebar_inputs$feature() else input$feature
            req(feature_val)
            size_val_multiplier <- input[[slider_id]]
            if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
            default_size <- default_pt_sizes()[s_local]
            final_size <- default_size * size_val_multiplier
            plot_obj <- obj()
            
            # Use SCT assay if available, otherwise use default assay
            if ("SCT" %in% names(plot_obj@assays)) {
              DefaultAssay(plot_obj) <- "SCT"
            }
            
            # Check if feature exists in the current assay
            current_assay <- DefaultAssay(plot_obj)
            available_features <- rownames(plot_obj[[current_assay]])
            
            if (!feature_val %in% available_features) {
              # Try to find the feature in other assays
              for (assay_name in names(plot_obj@assays)) {
                if (feature_val %in% rownames(plot_obj[[assay_name]])) {
                  DefaultAssay(plot_obj) <- assay_name
                  current_assay <- assay_name
                  break
                }
              }
            }
            
            if (feature_val %in% rownames(plot_obj[[current_assay]])) {
              p <- SpatialFeaturePlot(plot_obj, images = s_local, features = feature_val, pt.size.factor = final_size)
              p + labs(title = paste0("Expression level of ", feature_val, " in Sample ", s_local)) +
                theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
            } else {
              # Create an empty plot if feature not found
              ggplot() + 
                labs(title = paste0("Feature not found: ", feature_val, " in Sample ", s_local)) +
                theme_void() +
                theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
            }
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
