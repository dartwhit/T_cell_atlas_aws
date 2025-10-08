cat("✅ server.R initialized at", Sys.time(), "\n", file = stderr())

library(shiny)
library(Seurat)
library(DT)
library(ggplot2)
library(dplyr)
library(stringr)
library(VAM)
library(shinyjs)
library(tidyr)
source("setup.R")
source("modules/dataset_gallery_module.R")
source("modules/explore_sidebar_module.R")
source("modules/spatial_unit.R")
source("modules/scrna_seq_module.R")

options(shiny.trace = TRUE)


# Define server logic required to draw a histogram
server <- function(input, output,session) {
  # Get studies with spatial data
  studies_with_spatial <- reactive({
    names(dataset_files)[sapply(dataset_files, function(x) !is.null(x$spatial_seurat))]
  })
  
  # Initialize spatial sidebar
  spatial_sidebar_inputs <- spatial_sidebar_server("spatial_sidebar", 
                                                  studies_with_spatial = studies_with_spatial)
  
  # Get the path to the selected spatial data
  spatial_data_path <- reactive({
    req(spatial_sidebar_inputs$spatial_study_selector())
    paste0(inDir, dataset_files[[spatial_sidebar_inputs$spatial_study_selector()]][["spatial_seurat"]])
  })

  # Load spatial object for sharing between sidebar and main module
  spatial_obj <- reactive({
    req(spatial_data_path())
    readRDS(spatial_data_path())
  })
  
  # Update spatial sidebar with the spatial object when it's loaded
  observe({
    req(spatial_obj())
    updateSelectizeInput(
      session, 
      "spatial_sidebar-feature", 
      choices = rownames(spatial_obj()[["SCT"]]), 
      server = TRUE
    )
    updateCheckboxGroupInput(
      session,
      "spatial_sidebar-samples", 
      choices = names(spatial_obj()@images)
    )
  })

  spatial_server(
    id = "sp1",
    spat_obj = spatial_obj,
    sidebar_inputs = spatial_sidebar_inputs
  )

  options(shiny.trace = FALSE, shiny.fullstacktrace = FALSE, shiny.sanitize.errors = TRUE)

  selected_study_from_gallery <- dataset_gallery_server("gallery_module")

  sidebar_inputs <- explore_sidebar_server("explore_sidebar_module", 
                                           selected_study_from_gallery = selected_study_from_gallery,
                                           dataset_config = datasets)

  observeEvent(selected_study_from_gallery(), {
    study_info <- selected_study_from_gallery()
    req(study_info)

    if (is.list(study_info) && !is.null(study_info$view)) {
      if (study_info$view == "scrna") {
        # Navigate to the scRNA-seq explorer page
        nav_select("nav_page", selected = "Explore")
        # The explore_sidebar_module will automatically update its own study selector
      } else if (study_info$view == "spatial") {
        # Navigate to the spatial explorer page
        nav_select("nav_page", selected = "Spatial")
        # Update the study selector on the spatial page (now properly namespaced)
        updateSelectInput(session, "spatial_sidebar-spatial_study_selector", selected = study_info$id)
      }
    } else {
      # Fallback for old behavior or unexpected data
      nav_select("nav_page", selected = "Explore")
    }
  }, ignoreInit = TRUE)

  
  # # ########################## Starting page #######################
  # output$Tabib_img <- renderImage({

  #   screen_width <- input$dimension[1]
  #   image_width <- (screen_width / 2)*0.9

  #   list(src = "imgs/Tabib_img.png",
  #        width = paste0(image_width, "px"))
  # }, deleteFile = FALSE)

  # output$Gur_img <- renderImage({

  #   screen_width <- input$dimension[1]
  #   image_width <- (screen_width / 2)*0.9

  #   list(src = "imgs/Gur_img.png",
  #        width = paste0(image_width, "px"))
  # }, deleteFile = FALSE)

  # output$Ma_img <- renderImage({

  #   screen_width <- input$dimension[1]
  #   image_width <- (screen_width / 2)*0.9

  #   list(src = "imgs/Ma_img.png",
  #        width = paste0(image_width, "px"))
  # }, deleteFile = FALSE)

  # output$Khanna_img <- renderImage({

  #   screen_width <- input$dimension[1]
  #   image_width <- (screen_width / 2)*0.9

  #   list(src = "imgs/Khanna_img.png",
  #        width = paste0(image_width, "px"))
  # }, deleteFile = FALSE)
  
  # output$TMKMH_img <- renderImage({
  #   image_width <- input$dimentions[1]*0.9
  #   list(src = "imgs/TMKMH_img.png",
  #     width = paste0(image_width,"px")
  #   )
  # },
  # deleteFile = FALSE)
  
  


  ################################## Explore page #################################
  scrna_data <- scrna_seq_server(
    "scrna_seq_module", 
    sidebar_inputs = sidebar_inputs,
    dataset_files = dataset_files,
    inDir = inDir
  )
  
  # Observer to update the sidebar UI when the gene/pathway lists are loaded
  observe({
    req(scrna_data$gene_list())
    updateSelectizeInput(
      session, 
      "explore_sidebar_module-gene_select", 
      choices = scrna_data$gene_list(), 
      server = TRUE
    )
  })
  
  observe({
    req(scrna_data$pathway_list())
    updateSelectInput(
      session,
      "explore_sidebar_module-pathway_select",
      choices = scrna_data$pathway_list()
    )
  })

}
