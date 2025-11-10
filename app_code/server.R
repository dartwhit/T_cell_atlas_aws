cat("✅ server.R initialized at", Sys.time(), "\n", file = stderr())

library(shiny)
library(shinymanager)
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

options(shiny.trace = TRUE)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  db_path <- "/srv/shiny-server/atlas/data/users_lab.db"  # absolute! not relative
  cat("Shinymanager DB (app):", db_path, "\n")
  
  # Connect to the credentials db
  credentials <- read_db_decrypt(
    conn = db_path,
    name = 'credentials'
  )
  
  res_auth <- secure_server(
    check_credentials = credentials
  )

  output$welcome <- renderText(paste("Welcome,", res_auth$user))



  # Get studies with spatial data
  studies_with_spatial <- reactive({
    names(dataset_files)[sapply(dataset_files, function(x) !is.null(x$spatial_seurat))]
  })
  
  # Update the spatial study selector
  observe({
    updateSelectInput(session, "spatial_study_selector", choices = studies_with_spatial())
  })
  
  # Get the path to the selected spatial data
  spatial_data_path <- reactive({
    req(input$spatial_study_selector)
    paste0(inDir, dataset_files[[input$spatial_study_selector]][["spatial_seurat"]])
  })

    spatial_server(
    id = "sp1",
    rds_path = spatial_data_path
  )

  options(shiny.trace = FALSE, shiny.fullstacktrace = FALSE, shiny.sanitize.errors = TRUE)

  selected_study_from_gallery <- dataset_gallery_server("gallery_module")

  sidebar_inputs <- explore_sidebar_server("explore_sidebar_module", 
                                           selected_study_from_gallery = selected_study_from_gallery)

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
        nav_select("nav_page", selected = "spatial")
        # Update the study selector on the spatial page
        updateSelectInput(session, "spatial_study_selector", selected = study_info$id)
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
  
  
  observeEvent(input$explore_Tabib,{
    updateSelectInput(session, "study", selected = "tabib")
  })
  
  observeEvent(input$explore_Gur,{
    updateSelectInput(session, "study", selected = "gur")
  })
  
  observeEvent(input$explore_Ma,{
    updateSelectInput(session, "study", selected = "ma")
  })
  
  observeEvent(input$explore_Khanna,{
    updateSelectInput(session, "study", selected = "khanna")
  })
  
  observeEvent(input$explore_tmkmh,{
    updateSelectInput(session, "study", selected = "tmkmh")
  })

  ################################## Explore page #################################
  
  
  reset_trigger <- reactiveVal(FALSE)
  # seurat object to use
  seurat_obj <- reactiveVal(NULL)
  VAM_df <- reactiveVal()
  DEGs_df <- reactiveVal()
  gene_list_obj <- reactiveVal()
  cell_clusters <- reactiveVal(NULL)
  pathway_list <- reactiveVal()
  meta_df <- reactiveVal()
  
  DEG_path <- reactive({
    input$by_disease
    # req(sidebar_inputs$study(), sidebar_inputs$data_level(), input$by_disease, sidebar_inputs$anno())
    
    if (input$by_disease == FALSE) {
      if (sidebar_inputs$data_level() == "full") {
        if (sidebar_inputs$anno() == TRUE) {
          paste0(inDir, DE_dir, dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["DEGs_auto"]])
        } else {
          paste0(inDir, DE_dir, dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["DEGs_broad"]])
        }
      } else {
        paste0(inDir, DE_dir, dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["DEGs"]])
      }
    } else {
      if (sidebar_inputs$data_level() == "full") {
        if (sidebar_inputs$anno() == TRUE) {
          paste0(inDir, DE_dir, dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["DE_by_disease_auto"]])
        } else {
          paste0(inDir, DE_dir, dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["DE_by_disease_broad"]])
        }
      } else {
        paste0(inDir, DE_dir, dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["DE_by_disease"]])
      }
    }
  })
  
  VAM_path <- reactive({
    input$by_disease
    
    if(sidebar_inputs$data_level() == "full") {
      if(sidebar_inputs$anno() ==TRUE) {
        paste0(inDir,DE_dir,dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["VAM_by_disease_auto"]])
      } else {
        paste0(inDir,DE_dir,dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["VAM_by_disease_broad"]])
      }
    } else {
      paste0(inDir, DE_dir, dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["VAM_by_disease"]])
      
    }
  })
  
  
  # ## Refresh DE table when one of the input changes: study, data_level, anno, feature_type
  # observeEvent({sidebar_inputs$study(), sidebar_inputs$anno(), sidebar_inputs$data_level(), sidebar_inputs$feature_type()},
  #              {
  #                
  #                
  #              })
  #Utilize reset trigger
  

  ## Reload seurat object and DE tables when load button is clicked or annotation changed
  observeEvent({sidebar_inputs$load_btn()
    sidebar_inputs$anno()},
    {
    # Clear previous loaded dataset
    reset_trigger(TRUE)
    invalidateLater(0,session)
    reset_trigger(FALSE)
    seurat_obj(NULL)
    gc()
    
    req(sidebar_inputs$study(), sidebar_inputs$data_level())
    
    shinyjs::show(selector = "#explore_sidebar_module-hidden_menu")
    
    ################## get files paths of the selected dataset #############
    seurat_path <- paste0(inDir,dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["seurat"]])
    gene_list_path <- paste0(inDir,dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["gene_list"]])
    
    
    req(DEG_path())
    print(paste("Loading DEGs from", DEG_path()))
     deg_df <- if (endsWith(DEG_path(), ".csv")) {
      read.csv(DEG_path())
    } else {
      read.delim(DEG_path(), sep = "\t")
    }
    
    DEGs_df(deg_df)
    cell_clusters(levels(as.factor(deg_df$cluster)))
    updateSelectInput(session, "cell_cluster", choices = cell_clusters())
    
    meta_file <- dataset_files[[sidebar_inputs$study()]][["meta"]]
    metadata_path <- paste0(inDir,meta_file)

    ################## Load data ##################
    # Seurat object and gene lists
    print(seurat_path)
    seurat_obj(readRDS(seurat_path))
    gene_list_obj(readRDS(gene_list_path))
    
    vam_names <- rownames(seurat_obj()@assays$VAMcdf)
    pathway_names_for_display <- str_remove(vam_names,"HALLMARK-")
    pathway_choices <- setNames(vam_names, pathway_names_for_display)
    pathway_list(pathway_choices)
    
    
    # Try loading metadata,
    meta_data <- tryCatch(
      {
        if (!is.null(meta_file)) {
          read.delim(metadata_path)
        } else {
          NULL
        }
      },
      error = function(e) {
        NULL
      }
    )
    meta_df(meta_data)
    


    
    # Try loading VAM df
    # VAM_file <- dataset_files[[sidebar_inputs$study()]][[sidebar_inputs$data_level()]][["VAM_df"]]
    # VAM_path <- paste0(inDir,DE_dir,VAM_file)
    # print(VAM_path)
    
    VAM_data <- tryCatch(
      {
        if (!is.null(VAM_file)) {
          read.delim(VAM_path)
        } else {
          NULL
        }
      },
      error = function(e) {
        NULL
      }
    )
    # print(VAM_data)
    VAM_df(VAM_data)
    
    
    # update gene, pathway and cell type list
    # cell_clusters(levels(as.factor(data$cluster)))
    
    updateSelectizeInput(session, "explore_sidebar_module-gene_select", choices = gene_list_obj(), server = TRUE)
    updateSelectInput(session,"explore_sidebar_module-pathway_select", choices = pathway_list())
    
    
    print(paste("Loaded",sidebar_inputs$data_level(),"from dataset:", sidebar_inputs$study()))  
    showNotification(paste("Loading data for study:", sidebar_inputs$study(), 
                           "and data level:", sidebar_inputs$data_level()), type = "message")

    })
  
  observeEvent(input$by_disease, {
    req(sidebar_inputs$load_btn() > 0)
    req(DEG_path())
    
    print(paste("Reloading DEGs from", DEG_path()))
     deg_df <- if (endsWith(DEG_path(), ".csv")) {
      read.csv(DEG_path())
    } else {
      read.delim(DEG_path(), sep = "\t")
    }
    
    DEGs_df(deg_df)
    cell_clusters(levels(as.factor(deg_df$cluster)))
    updateSelectInput(session, "cell_cluster", choices = cell_clusters())
  })
  
  # trigger reset when changing study or dataset level
  observeEvent(list(sidebar_inputs$study(),sidebar_inputs$data_level()),{
    reset_trigger(TRUE)
    invalidateLater(0, session)
    reset_trigger(FALSE)
    
    # Hide the sidebar options
    shinyjs::hide(selector = "#explore_sidebar_module-hidden_menu")
    
    # Reset reactive values
    seurat_obj(NULL)
    gene_list_obj(NULL)
    VAM_df(NULL)
    DEGs_df(NULL)
    cell_clusters(NULL)
    pathway_list(NULL)
    meta_df(NULL)
    
    #Clear gene/pathway input fields
    updateSelectizeInput(session, "explore_sidebar_module-gene_select", choices = NULL, selected = character(0))
    updateTextInput(session, "explore_sidebar_module-gene_input",value = "")
    
    
    updateSelectInput(session, "explore_sidebar_module-pathway_select", choices = NULL, selected = character(0))
    
    
    updateSelectInput(session, "cell_cluster", choices = NULL, selected = character(0))

  })
  
  
  
  # Update default assay of the seurat object
  observe({
    obj <- seurat_obj()
    
    # Ensure the object is not NULL
    if (!is.null(obj)) {
      # Change the default assay based on input
      if (sidebar_inputs$feature_type() == "Genes") {
        DefaultAssay(obj) <- "RNA"
      } else if (sidebar_inputs$feature_type() == "Pathways") {
        DefaultAssay(obj) <- "VAMcdf"
      }
      
      # Update the reactiveVal or reactive expression with the modified object
      seurat_obj(obj)
    }
  })
  
  
  # ------------ UMAP of the full dataset ---------------
  full_umap <- reactive({
    req(seurat_obj())
    if(sidebar_inputs$anno() == TRUE){
      DimPlot(seurat_obj(), group.by = "seurat_clusters")
    }else{
      DimPlot(seurat_obj())
    }
    
  })
  
  download_umap <- function(input, output){
    full_umap()
  }
  
  downloadablePlot(id = "umap_dld",
                   filenameroot = "umap",
                   downloadfxns = list(png = download_umap),
                   visibleplot = download_umap
  )

  
  # ---------- Featureplot of the queried gene ----------
  gene_queried <- reactive({
    
    # Check reset trigger
    if (reset_trigger()) {
      return(NULL)
    }
    if (input$update_gene_queried){
      new_gene_queried()
    }
    
    if (!sidebar_inputs$use_textinput() && length(sidebar_inputs$gene_select()) > 0) {
      return(sidebar_inputs$gene_select())
    } else if (sidebar_inputs$use_textinput() && nchar(sidebar_inputs$gene_input()) > 0) {
      # Splitting by comma, tab, space, or newline
      genes <- unlist(strsplit(sidebar_inputs$gene_input(), "[,\n]+"))
      genes <- genes[genes != ""]  # Remove any empty strings
      return(genes)
    } else {
      return(NULL)
    }
  })
  
  pathway_queried <- reactive({    # Check reset trigger
    if (reset_trigger()) {
      return(NULL)
    }
    if (input$update_gene_queried){
      if(sidebar_inputs$feature_type() == "Pathways"){
        new_gene_queried()
        
      }
    }
    if (!sidebar_inputs$use_textinput_VAM() && length(sidebar_inputs$pathway_select()) > 0) {
      return(sidebar_inputs$pathway_select())
    } else if (sidebar_inputs$use_textinput_VAM() && nchar(sidebar_inputs$pathway_input()) > 0) {
      # Splitting by comma, tab, space, or newline
      pathways <- unlist(strsplit(sidebar_inputs$pathway_input(), "[,\n]+"))
      pathways <- pathways[pathways != ""]  # Remove any empty strings
      return(pathways)
    } else {
      return(NULL)
    }
  })
  
  ### Append user-defined pathway
  observeEvent(sidebar_inputs$submit_geneset(),{
    curr_obj <- seurat_obj()
    if (sidebar_inputs$geneset_name() != "" && nchar(sidebar_inputs$VAM_geneset()>0)){
      new_geneset_name <- sidebar_inputs$geneset_name()
      # Parse gene list
      new_geneset <- unlist(strsplit(sidebar_inputs$VAM_geneset(), "[,\n]+"))
      new_geneset <- new_geneset[new_geneset != ""]  # Remove any empty strings
      print(new_geneset)  # Debug print to console
      
      # Create collection for VAM
      new_collection <- list()
      new_collection[[new_geneset_name]] <- new_geneset
      
      all_genes<-rownames(curr_obj@assays$RNA@counts)
      new_collection <- createGeneSetCollection(gene.ids=all_genes, gene.set.collection = new_collection)
      
      vam_res <- vamForCollection(t(curr_obj[["RNA"]]@counts), new_collection)
      
      # Add VAM scores of new pathway to object
      new_VAM_dist <- CreateAssayObject(rbind(curr_obj[["VAMdist"]]@counts, t(vam_res$distance.sq)))
      new_VAM_score <- CreateAssayObject(rbind(curr_obj[["VAMcdf"]]@counts, t(vam_res$cdf.value)))
      curr_obj[["VAMdist"]] <- new_VAM_dist
      curr_obj[["VAMcdf"]] <-new_VAM_score
      
      # Update seurat object
      seurat_obj(curr_obj)
      
      vam_names <- rownames(curr_obj@assays$VAMcdf)
      pathway_names_for_display <- str_remove(vam_names,"HALLMARK-")
      pathway_choices <- setNames(vam_names, pathway_names_for_display)
      pathway_list(pathway_choices)
      
      updateSelectInput(session,"explore_sidebar_module-pathway_select", choices = pathway_list(), selected = new_geneset_name)
    }
  })
  
  output$show_switch <- reactive({
     (sidebar_inputs$feature_type() == "Genes") && (length(gene_queried()) > 3)
  })
  outputOptions(output, "show_switch", suspendWhenHidden= FALSE)
  

  
  featureplot_plot_gene <- reactive({
    if (!is.null(gene_queried())) {
      selected_genes <- gene_queried()
      
      if(length(selected_genes) == 2) {
        tryCatch({
          FeaturePlot(
            object = seurat_obj(),
            features = selected_genes,
            # cols = c("blue", "red"),
            blend = TRUE,
            pt.size = 0.3
          )
        }, error = function(e) {
          print(paste("Error in FeaturePlot with blending:", e$message))
        })
      } else {
        tryCatch({
          FeaturePlot(
            object = seurat_obj(),
            features = selected_genes,
            pt.size = 0.3
          )
        }, error = function(e) {
          print(paste("Error in FeaturePlot:", e$message))
        })
      }
    }
  })
  
  
  featureplot_plot_pathway <- reactive({
    if (!is.null(pathway_queried())) {
      selected_pathways <- pathway_queried()
      
      if(length(selected_pathways) == 2) {
        tryCatch({
          FeaturePlot(
            object = seurat_obj(),
            features = selected_pathways,
            # cols = c("blue", "red"),
            blend = TRUE,
            pt.size = 0.3
          )
        }, error = function(e) {
          print(paste("Error in FeaturePlot with blending:", e$message))
        })
      } else {
        tryCatch({
          FeaturePlot(
            object = seurat_obj(),
            features = selected_pathways,
            pt.size = 0.3
          )
        }, error = function(e) {
          print(paste("Error in FeaturePlot:", e$message))
        })
      }
    }
  })
  
  
  output$featurePlot <- renderPlot({
    if(sidebar_inputs$feature_type() == "Genes"){
      featureplot_plot_gene()
    }else{
      featureplot_plot_pathway()
    }
    
  })
  
  output$featurePlot_msg <- renderText({
    "Please select a gene to display the plot."
  })
  
  output$featurePlot_UI <- renderUI({
    if (!is.null(gene_queried()) || !is.null(pathway_queried())) {
      tagList(
        plotOutput("featurePlot"),
        downloadButton("dld_feature_plot", "Download"),
        selectInput("file_type_1", "Plot format", choices = c("png", "pdf", "jpg"))
      )
    } else {
      textOutput("featurePlot_msg")
    }
  })
  
  #### Gene expression levels comparison ######


  
  output$geneExpPlot <- renderPlot({
    # Initialize variables
    feature_names <- NULL
    curr_obj <- seurat_obj()
    f_type <- sidebar_inputs$feature_type()
    
    # Determine feature names based on input
    if (f_type == "Genes") {
      feature_names <- gene_queried()  # Assuming this returns a vector of genes
    } else {
      feature_names <- pathway_queried()  # Assuming this returns a vector of pathways
    }
    
    # Ensure required inputs are present
    req(curr_obj, feature_names)
    
    # Check if all features exist in the Seurat object
    if (!all(feature_names %in% rownames(curr_obj))) {
      missing_features <- feature_names[!feature_names %in% rownames(curr_obj)]
      stop(paste("Features queried not found in seurat_obj:", paste(missing_features, collapse = ", ")))
    }
    
    # Determine the type of plot to render
    if (f_type == "Pathways" && length(feature_names) > 3) {
      # Show heatmap for >3 pathways
      DoHeatmap(curr_obj, 
                features = feature_names, 
                assay = "VAMcdf", 
                slot = "data")
    } else if (f_type == "Genes" && length(feature_names) > 3 && input$heatmap) {
      # Show heatmap for >3 genes if heatmap option is selected
      DoHeatmap(curr_obj, 
                features = feature_names, 
                assay = "RNA", 
                slot = "data")
    } else {
      # Show default dot plot or violin plot
      if (length(feature_names) <= 3) {
        # Show a violin plot for <= 3 features
        if (is.null(input$plot_type) || input$plot_type) { # Default to VlnPlot
          assay_to_use <- if (f_type == "Genes") {
            if ("SCT" %in% Assays(curr_obj)) "SCT" else "RNA"
          } else {
            "VAMcdf"
          }
          VlnPlot(curr_obj, 
                  features = feature_names, 
                  assay = assay_to_use,
                  split.by = "Disease", 
                  split.plot = TRUE,
                  pt.size = 0,
                  cols = c("HC" = "#cb07a4ff", "SSc" = "#09b646ff"))
        } else { # BoxPlot
          assay_to_use <- if (f_type == "Genes") {
            if ("SCT" %in% Assays(curr_obj)) "SCT" else "RNA"
          } else {
            "VAMcdf"
          }
          
          plot_data <- FetchData(curr_obj, vars = c(feature_names, "Disease"), assay = assay_to_use) %>%
            pivot_longer(cols = -Disease, names_to = "feature", values_to = "expression")
          
          ggplot(plot_data, aes(x = Disease, y = expression, fill = Disease)) +
            geom_boxplot(outlier.shape = NA) +
            scale_fill_manual(values = c("HC" = "#cb07a4ff", "SSc" = "#09b646ff")) +
            facet_wrap(~feature, scales = "free_y") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      } else {
        # Show a dot plot for > 3 features by default
        DotPlot(curr_obj, 
                features = feature_names, 
                cols = c("blue", "red"), 
                split.by = "Disease")
      }
    }
  })
  
  

  output$geneExpPlot_msg <- renderText({
    "Please choose a gene."
  })
  
  output$geneExpPlot_Khanna_msg <- renderText({
    "Not applicable"
  })

  output$exp_plot_UI <- renderUI({
    if (!is.null(gene_queried()) || !is.null(pathway_queried())) {
      if(sidebar_inputs$study() == "khanna"){
        textOutput("geneExpPlot_Khanna_msg")
      }else{
        plotOutput("geneExpPlot")
      }
      
    } else {
      textOutput("geneExpPlot_msg")
    }
  })

  output$dld_feature_plot <- downloadHandler(
    filename = function(){
      generated_filename <- paste0("featureplot.",input$file_type_1)
      print(paste("Generated filename:", generated_filename))
      return(generated_filename)  # Ensure the filename is returned
      "featureplot.pdf"
    },
    content = function(file){
      filetype <- tools::file_ext(file)
      if (sidebar_inputs$feature_type() == "Genes") {
        p <- plot(featureplot_plot_gene())
      } else {
        p <- plot(featureplot_plot_pathway())
      }
      
      ggsave(file,plot = p, device = filetype)
      
    
      # 
      # if (filetype == "png"){
      #   png(file)
      #   print("here")
      #   p
      #   dev.off()
      # }else if (filetype == "pdf"){
      #   pdf(file)
      # }else if (filetype == "jpg"){
      #   jpeg(file)
      # }
      png(file)
      print("here")
      p
      dev.off()

      print(paste("File type:", filetype))
    }
  )

  
  # Show DEGs of the selected cell cluster
  cell_cluster_selected <- reactive({
    input$cell_cluster
  })


  DEGs_df_show <- reactive({
    
    req(sidebar_inputs$study(), sidebar_inputs$feature_type())
    
    if(sidebar_inputs$feature_type() == "Genes"){
      req(DEGs_df())
      
      DEGs_df() %>%
        dplyr::filter(cluster == cell_cluster_selected()) %>%
        dplyr::select(gene, avg_log2FC, p_val, p_val_adj) %>%
        dplyr::mutate(avg_log2FC = round(avg_log2FC,3)
        )
    }else if(sidebar_inputs$feature_type() == "Pathways"){
      req(VAM_df())
      VAM_df() %>%
        dplyr::filter(cluster == cell_cluster_selected()) %>%
        dplyr::select(gene, avg_log2FC, p_val, p_val_adj) %>%
        dplyr::mutate(avg_log2FC = round(avg_log2FC, 3))
    }
    
  })

  # output$DEGs_df<- renderDT({
  #   req(DEGs_df(), cancelOutput = TRUE) # prevent execution if data is null
  #   DEGs_df_show()
  # }
  # )
  
  
  output$marker_tbl <- renderDT({
    DEGs_df_show()
  })
  
  output$DEGs_table_ui <- renderUI({
    if (!is.null(DEGs_df_show())) {
       tagList(
         DT::dataTableOutput("marker_tbl"),
         br(),
         downloadButton("downloadDEGs", "Download Marker List")
       )
    } else {
      NULL
    }
  })
  
  output$downloadDEGs <- downloadHandler(
    filename = function(){
      paste0(sidebar_inputs$study(),"_DEGs_",sidebar_inputs$data_level(),".txt")
    },
    content = function(file){
      write.table(DEGs_df_show(), file)
    }
  )
  
  output$DEGs_msg <- renderUI({
    if (is.null(DEGs_df()) && sidebar_inputs$feature_type() == "Genes") {
      div(style = "color: red; font-weight: bold; text-align: center;",
          "Please click the 'Load Data' button to display the Genes table.")
    } else if (is.null(VAM_df()) && sidebar_inputs$feature_type() == "Pathways") {
      div(style = "color: red; font-weight: bold; text-align: center;",
          "Please click the 'Load Data' button to display the Pathways table.")
    }
  })
  
  
  # Update featureplot using the selected DEG and jump to UMAP explorer
  new_gene_queried <- reactive({
    idx <- input$marker_tbl_rows_selected
    print(idx)
    
    new_gene <- DEGs_df_show()[idx,]$gene
    # if(sidebar_inputs$feature_type() == "Genes"){
    #   new_gene <- DEGs_df()[idx,]$gene
    # }else if(sidebar_inputs$feature_type() == "Pathways"){
    #   new_gene <- VAM_df()[idx,]$gene
    #   print(paste0("New pathway: ", new_gene))
    # }
    new_gene
  })

  observeEvent(input$update_gene_queried,{
    
    nav_select("explore_tabs", selected = "plots")
    
    # Create branches based on feature type
    feature_type_selected <- sidebar_inputs$feature_type()
    print(feature_type_selected)
    updateAwesomeRadio(session,"feature_type", selected = feature_type_selected)
    
    
    input_id <- if (feature_type_selected == "Genes") "explore_sidebar_module-gene_select" else "explore_sidebar_module-pathway_select"
    if (feature_type_selected == "Genes") {
      choices <- gene_list_obj()
      selected_value <- as.character(new_gene_queried())
    } else {
      choices <- pathway_list()
      selected_name <- as.character(new_gene_queried())
      selected_value <- choices[names(choices) == selected_name]
    }
    
    if (is.null(choices)) {
      choices <- character(0)
    }
    
    updateSelectizeInput(session,
                         # "gene_select",
                         input_id,
                         selected = selected_value,
                         choices = choices,
                         server = TRUE
    )
    
})
  


  output$meta_df <- renderUI({
    if (is.null(meta_df())) {
      textOutput("meta_data_msg")
    } else {
      DTOutput("meta_table")
    }
  })
  
  output$meta_table <- renderDT({
    meta_df()
  })
  
  output$meta_data_msg <- renderText({
    "Metadata not available for the selected dataset."
  })


}