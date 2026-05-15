library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(tidyr)

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
    ),
    # Expression by cluster/region
    conditionalPanel(
      condition = paste0("input['", ns("feature"), "'] && input['", ns("feature"), "'].length > 0"),
      card(
        card_header(
          div(
            class = "d-flex align-items-center justify-content-between w-100",
            span("Expression by cluster / region"),
            div(
              class = "d-flex align-items-center gap-2",
              radioButtons(
                ns("expr_plot_type"), label = NULL,
                choices = c("Violin" = "violin", "Box" = "box"),
                selected = "violin", inline = TRUE
              ),
              downloadButton(ns("dld_expr"), "Download", class = "btn-sm")
            )
          )
        ),
        plotOutput(ns("expr_by_cluster"), height = 420)
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
    
    # ---- Expression by cluster/region plot ----
    group_by_col <- reactive({
      # Use the raw input$group_by value directly — it matches the metadata column
      # name (the same value passed to SpatialDimPlot/DimPlot group.by).
      # Fall back to the first metadata column if it somehow doesn't exist.
      grp <- input$group_by
      meta_cols <- colnames(obj()@meta.data)
      if (!grp %in% meta_cols) {
        # Try common alternatives
        alts <- c("seurat_clusters", "ident", "orig.ident")
        grp <- alts[alts %in% meta_cols][1]
      }
      grp
    })

    expr_plot_data <- reactive({
      req(obj(), length(input$feature) > 0)
      genes <- input$feature
      grp   <- group_by_col()
      req(!is.na(grp), grp %in% colnames(obj()@meta.data))
      df <- FetchData(obj(), vars = c(grp, genes))
      colnames(df)[1] <- "group"
      df$group <- factor(df$group)
      if (length(genes) == 1) {
        df$gene <- genes
        colnames(df)[2] <- "expression"
      } else {
        df <- tidyr::pivot_longer(df, cols = -group, names_to = "gene", values_to = "expression")
      }
      df
    })

    expr_plot_reactive <- reactive({
      req(expr_plot_data())
      df   <- expr_plot_data()
      grp_label <- switch(input$group_by,
        "cluster"     = "Seurat cluster",
        "Key_Regions" = "Key region",
        input$group_by
      )
      p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = expression, fill = group))
      if (isTRUE(input$expr_plot_type == "box")) {
        p <- p + ggplot2::geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4)
      } else {
        p <- p +
          ggplot2::geom_violin(scale = "width", trim = TRUE) +
          ggplot2::geom_jitter(width = 0.15, size = 0.3, alpha = 0.25, show.legend = FALSE)
      }
      p <- p +
        ggplot2::facet_wrap(~gene, scales = "free_y") +
        ggplot2::labs(x = grp_label, y = "Normalized expression", fill = grp_label) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.background = ggplot2::element_rect(fill = "#f0f0f0")
        )
      p
    })

    output$expr_by_cluster <- renderPlot({ expr_plot_reactive() })

    output$dld_expr <- downloadHandler(
      filename = function() {
        genes <- paste(input$feature, collapse = "-")
        paste0("expr_by_", input$group_by, "_", genes, ".png")
      },
      content = function(file) {
        n_genes <- length(input$feature)
        w <- max(7, 3.5 * ceiling(sqrt(n_genes)))
        h <- max(5, 3   * ceiling(n_genes / ceiling(sqrt(n_genes))))
        ggplot2::ggsave(file, expr_plot_reactive(), width = w, height = h, dpi = 150)
      }
    )

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
            size_val_multiplier <- input[[slider_id]]
            if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
            default_size <- default_pt_sizes()[s_local]
            final_size <- default_size * size_val_multiplier
            grp <- if (!is.null(input$group_by)) input$group_by else NULL
            SpatialDimPlot(obj(), images = s_local, group.by = grp, pt.size.factor = final_size)
          })
          
          # Reactive for FeaturePlot
          feat_plot_reactive <- reactive({
            req(input$feature)
            size_val_multiplier <- input[[slider_id]]
            if (is.null(size_val_multiplier)) size_val_multiplier <- 1.0
            default_size <- default_pt_sizes()[s_local]
            final_size <- default_size * size_val_multiplier
            plot_obj <- obj()
            DefaultAssay(plot_obj) <- "SCT"
            SpatialFeaturePlot(plot_obj, images = s_local, features = input$feature, pt.size.factor = final_size)
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
