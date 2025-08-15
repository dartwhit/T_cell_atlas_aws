cat("✅ server.R initialized at", Sys.time(), "\n", file = stderr())

library(shiny)
library(Seurat)
library(DT)
library(ggplot2)
library(dplyr)
library(stringr)
library(VAM)
library(shinyjs)
source("setup.R")

options(shiny.trace = TRUE)


# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  
  # ########################## Starting page #######################
  lapply(dataset_info$id, function(id) {
    output[[paste0(id, "_img")]] <- renderImage({
      screen_width <- input$dimension[1]
      img_path <- file.path("imgs", dataset_info$image[dataset_info$id == id])
      is_featured <- dataset_info$featured[dataset_info$id == id]
      image_width <- if (is_featured) screen_width * 0.9 else (screen_width / 2) * 0.9
      list(src = img_path, width = paste0(image_width, "px"))
    }, deleteFile = FALSE)

    observeEvent(input[[paste0("explore_", id)]], {
      updateSelectInput(session, "study", selected = id)
    })
  })
  
  ################################## Explore page #################################
  

  
  output$data_level_ui <- renderUI({
    selectInput("data_level",
                "Select data to visualize",
                choices = data_level_choices[[input$study]])
  })
  
  
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
    # req(input$study, input$data_level, input$by_disease, input$anno)
    
    if (input$by_disease == FALSE) {
      if (input$data_level == "full") {
        if (input$anno == TRUE) {
          paste0(inDir, DE_dir, dataset_files[[input$study]][[input$data_level]][["DEGs_auto"]])
        } else {
          paste0(inDir, DE_dir, dataset_files[[input$study]][[input$data_level]][["DEGs_broad"]])
        }
      } else {
        paste0(inDir, DE_dir, dataset_files[[input$study]][[input$data_level]][["DEGs"]])
      }
    } else {
      if (input$data_level == "full") {
        if (input$anno == TRUE) {
          paste0(inDir, DE_dir, dataset_files[[input$study]][[input$data_level]][["DE_by_disease_auto"]])
        } else {
          paste0(inDir, DE_dir, dataset_files[[input$study]][[input$data_level]][["DE_by_disease_broad"]])
        }
      } else {
        paste0(inDir, DE_dir, dataset_files[[input$study]][[input$data_level]][["DE_by_disease"]])
      }
    }
  })
  
  VAM_path <- reactive({
    input$by_disease
    
    if(input$data_level == "full") {
      if(input$anno ==TRUE) {
        paste0(inDir,DE_dir,dataset_files[[input$study]][[input$data_level]][["VAM_by_disease_auto"]])
      } else {
        paste0(inDir,DE_dir,dataset_files[[input$study]][[input$data_level]][["VAM_by_disease_broad"]])
      }
    } else {
      paste0(inDir, DE_dir, dataset_files[[input$study]][[input$data_level]][["VAM_by_disease"]])
      
    }
  })
  
  
  # ## Refresh DE table when one of the input changes: study, data_level, anno, feature_type
  # observeEvent({input$study, input$anno, input$data_level, input$feature_type},
  #              {
  #                
  #                
  #              })
  #Utilize reset trigger
  

  ## Reload seurat object and DE tables when load button is clicked or annotation changed
  observeEvent({input$load_btn
    input$anno},{
    # Clear previous loaded dataset
    reset_trigger(TRUE)
    invalidateLater(0,session)
    reset_trigger(FALSE)
    seurat_obj(NULL)
    gc()
    
    req(input$study, input$data_level)
    
    shinyjs::show("hidden_menu")
    
    ################## get files paths of the selected dataset #############
    seurat_path <- paste0(inDir,dataset_files[[input$study]][[input$data_level]][["seurat"]])
    gene_list_path <- paste0(inDir,dataset_files[[input$study]][[input$data_level]][["gene_list"]])
    
    
    req(DEG_path())
    print(paste("Loading DEGs from", DEG_path()))
    deg_df <- if (grepl("\\.csv$", DEG_path())) {
      read.csv(DEG_path())
    } else {
      read.delim(DEG_path(), sep = "\t")
    }
    
    DEGs_df(deg_df)
    cell_clusters(levels(as.factor(deg_df$cluster)))
    updateSelectInput(session, "cell_cluster", choices = cell_clusters())
    
    meta_file <- dataset_files[[input$study]][["meta"]]
    metadata_path <- paste0(inDir,meta_file)

    ################## Load data ##################
    # Seurat object and gene lists
    seurat_obj(readRDS(seurat_path))
    gene_list_obj(readRDS(gene_list_path))
    
    vam_names <- rownames(seurat_obj()@assays$VAMcdf)
    vam_names <- str_remove(vam_names,"HALLMARK-")
    pathway_list()
    
    
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
    # VAM_file <- dataset_files[[input$study]][[input$data_level]][["VAM_df"]]
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
    
    updateSelectizeInput(session, "gene_select", choices = gene_list_obj(), server = TRUE)
    updateSelectInput(session,"pathway_select", choices = rownames(seurat_obj()@assays$VAMdist))
    # updateSelectInput(session, "cell_cluster", choices = cell_clusters())
    updateSelectInput(session,"pathway_select", choices = pathway_list())
    
    
    print(paste("Loaded",input$data_level,"from dataset:", input$study))  
    showNotification(paste("Loading data for study:", input$study, 
                           "and data level:", input$data_level), type = "message")

  })
  
  observeEvent(input$by_disease, {
    req(input$load_btn > 0)
    req(DEG_path())
    
    print(paste("Reloading DEGs from", DEG_path()))
    deg_df <- if (grepl("\\.csv$", DEG_path())) {
      read.csv(DEG_path())
    } else {
      read.delim(DEG_path(), sep = "\t")
    }
    
    DEGs_df(deg_df)
    cell_clusters(levels(as.factor(deg_df$cluster)))
    updateSelectInput(session, "cell_cluster", choices = cell_clusters())
  })
  
  # trigger reset when changing study or dataset level
  observeEvent(list(input$study,input$data_level),{
    reset_trigger(TRUE)
    invalidateLater(0, session)
    reset_trigger(FALSE)
    
    #Clear gene/pathway input fields
    updateSelectizeInput(session, "gene_select", choices = NULL, selected = character(0))
    updateTextInput(session, "gene_input",value = "")
    
    
    updateSelectInput(session, "pathway_select", choices = NULL, selected = character(0))
    updateTextInput(session, "pathway_input", value = "")
    
    
    updateSelectInput(session, "cell_cluster", choices = NULL, selected = character(0))

  })
  
  
  
  # Update default assay of the seurat object
  observe({
    obj <- seurat_obj()
    
    # Ensure the object is not NULL
    if (!is.null(obj)) {
      # Change the default assay based on input
      if (input$feature_type == "Genes") {
        DefaultAssay(obj) <- "RNA"
      } else if (input$feature_type == "Pathways") {
        DefaultAssay(obj) <- "VAMcdf"
      }
      
      # Update the reactiveVal or reactive expression with the modified object
      seurat_obj(obj)
    }
  })
  
  
  # ------------ UMAP of the full dataset ---------------
  full_umap <- reactive({
    req(seurat_obj())
    if(input$anno == TRUE){
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
    
    if (!input$use_textinput && length(input$gene_select) > 0) {
      return(input$gene_select)
    } else if (input$use_textinput && nchar(input$gene_input) > 0) {
      # Splitting by comma, tab, space, or newline
      genes <- unlist(strsplit(input$gene_input, "[,\n]+"))
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
      if(input$feature_type == "Pathways"){
        new_gene_queried()
        
      }
    }
    if (!input$use_textinput_VAM && length(input$pathway_select) > 0) {
      return(input$pathway_select)
    } else if (input$use_textinput_VAM && nchar(input$pathway_input) > 0) {
      # Splitting by comma, tab, space, or newline
      pathways <- unlist(strsplit(input$pathway_input, "[,\n]+"))
      pathways <- pathways[pathways != ""]  # Remove any empty strings
      return(pathways)
    } else {
      return(NULL)
    }
  })
  
  ### Append user-defined pathway
  observeEvent(input$submit_geneset,{
    curr_obj <- seurat_obj()
    if (input$geneset_name != "" && nchar(input$VAM_geneset>0)){
      new_geneset_name <- input$geneset_name
      # Parse gene list
      new_geneset <- unlist(strsplit(input$VAM_geneset, "[,\n]+"))
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
      pathway_list(rownames(curr_obj@assays$VAMcdf))
      updateSelectInput(session,"pathway_select", choices = pathway_list())
    }
  })
  
  output$show_switch <- reactive({
     (input$feature_type == "Genes") && (length(gene_queried()) > 3)
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
    if(input$feature_type == "Genes"){
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
    f_type <- input$feature_type
    
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
        VlnPlot(curr_obj, 
                features = feature_names, 
                # assay = "SCT",
                split.by = "Disease", 
                split.plot = TRUE)
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
      if(input$study == "khanna"){
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
      if (input$feature_type == "Genes") {
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
    
    req(input$study, input$feature_type)
    
    if(input$feature_type == "Genes"){
      req(DEGs_df())
      
      DEGs_df() %>%
        dplyr::filter(cluster == cell_cluster_selected()) %>%
        dplyr::select(gene, avg_log2FC, p_val, p_val_adj) %>%
        dplyr::mutate(avg_log2FC = round(avg_log2FC,3)
        )
    }else if(input$feature_type == "Pathways"){
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
      paste0(input$study,"_DEGs_",input$data_level,".txt")
    },
    content = function(file){
      write.table(DEGs_df_show(), file)
    }
  )
  
  output$DEGs_msg <- renderUI({
    if (is.null(DEGs_df()) && input$feature_type == "Genes") {
      div(style = "color: red; font-weight: bold; text-align: center;",
          "Please click the 'Load Data' button to display the Genes table.")
    } else if (is.null(VAM_df()) && input$feature_type == "Pathways") {
      div(style = "color: red; font-weight: bold; text-align: center;",
          "Please click the 'Load Data' button to display the Pathways table.")
    }
  })
  
  
  # Update featureplot using the selected DEG and jump to UMAP explorer
  new_gene_queried <- reactive({
    idx <- input$marker_tbl_rows_selected
    print(idx)
    
    new_gene <- DEGs_df_show()[idx,]$gene
    # if(input$feature_type == "Genes"){
    #   new_gene <- DEGs_df()[idx,]$gene
    # }else if(input$feature_type == "Pathways"){
    #   new_gene <- VAM_df()[idx,]$gene
    #   print(paste0("New pathway: ", new_gene))
    # }
    new_gene
  })

  observeEvent(input$update_gene_queried,{
    
    nav_select("explore_tabs", selected = "plots")
    
    # Create branches based on feature type
    feature_type_selected <- input$feature_type
    print(feature_type_selected)
    updateAwesomeRadio(session,"feature_type", selected = feature_type_selected)
    
    
    input_id <- ifelse(feature_type_selected == "Genes", "gene_select", "pathway_select") # Get correct input id for genes and pathways
    if (feature_type_selected == "Genes") {
      choices <- gene_list_obj()
    } else {
      choices <- pathway_list()
    }
    
    if (is.null(choices)) {
      choices <- character(0)
    }
    
    updateSelectizeInput(session,
                         # "gene_select",
                         input_id,
                         selected = as.character(new_gene_queried()),
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


