#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#test change

library(shiny)
library(shinyWidgets)
library(bslib)
library(shinycssloaders)
library(DT)
library(bsicons)
library(Seurat)
library(periscope2)
library(shinyjs)
source("setup.R")


# Feature plot of queried genes
card_feature <-card("Featureplot of selected genes",
                          withSpinner(
                            uiOutput("featurePlot_UI")
                          )
)
# UMAP
card_umap <- card("UMAP of the dataset",
                        # withSpinner(
                        #   plotOutput("full_umap")
                        # )

                    downloadablePlotUI(id                 = "umap_dld", 
                                       downloadtypes      = c("png"), 
                                       download_hovertext = "Download the umap here!",
                                       height             = "500px", 
                                       btn_halign         = "left")
                  
                  
)

# Expression of queried genes
card_gene_plot <-card("Plot of selected genes",
                            conditionalPanel("output.show_switch == true",
                                             switchInput("heatmap","Show heatmap")
                                             ),
                            
                            withSpinner(
                              
                              uiOutput("exp_plot_UI")
                              )
)


# DEGs table 
card_tabib_DEGs <- card(
  layout_sidebar(
    sidebar = sidebar(
      title = "Differential Expression setup",
      position = "left",
      selectInput("cell_cluster", "Choose a cell cluster",choices = NULL),
      checkboxInput("by_disease","Compare by disease"),
      actionButton("update_gene_queried","Show gene in UMAP space")
    ),
    uiOutput("DEGs_msg"),
    uiOutput("DEGs_table_ui")  )
  
)




# ########################################### Define UI ########################################
ui <- page_navbar(
  useShinyjs(), # Enable shinyjs

  tags$head(
    tags$script('
    $(document).on("shiny:connected", function(e) {
      dimension = [window.innerWidth, window.innerHeight];
      Shiny.onInputChange("dimension", dimension);
      $(window).resize(function(e) {
        dimension = [window.innerWidth, window.innerHeight];
        Shiny.onInputChange("dimension", dimension);
      });
      // Switch to the Explore tab when any explore button is clicked
      $("button[id^='explore_']").on("click", function() {
        $("a[data-value=\'Explore\']").tab("show");
      });
    });
  ')
  ),
  theme = bs_theme(bootswatch = "flatly"),
  # Page title
  title = "SSc cell atlas",
  id = "nav_page",
  nav_panel("Get started",
            {
              featured <- dataset_info[dataset_info$featured, ][1,]
              card(
                full_screen = TRUE,
                card_header(featured$title),
                actionButton(paste0("explore_", featured$id), "Explore"),
                imageOutput(paste0(featured$id, "_img"))
              )
            },
            do.call(
              layout_column_wrap,
              c(
                list(width = 1/4),
                lapply(
                  split(dataset_info[!dataset_info$featured, ],
                        seq_len(nrow(dataset_info[!dataset_info$featured, ]))),
                  function(ds) {
                    card(
                      full_screen = TRUE,
                      card_header(ds$title),
                      actionButton(paste0("explore_", ds$id), "Explore"),
                      imageOutput(paste0(ds$id, "_img")),
                      if (!is.na(ds$link) && nzchar(ds$link))
                        a(
                          href = ds$link,
                          p(ds$citation),
                          style = "color:grey", target = "_blank"
                        )
                    )
                  }
                )
              )
            )
  ),
  
  
  
  # Explore page
    nav_panel("Explore",
            layout_sidebar(
              sidebar = sidebar(
                title = "Select dataset and genes",
                position = "left",
                selectInput(
                  "study", "Select study to explore",
                  choices = setNames(dataset_info$id, dataset_info$short_name)
                ),

                
                conditionalPanel(
                  condition="input.data_level == 'full'",
                  switchInput("anno",
                              "Original seurat clusters",
                              value = FALSE)
                ),
                # Select object to visualize
                # selectInput("data_level","Select data to visualize",
                #             choices = c("Full" = "full",
                #                         "Fibroblasts" = "fib",
                #                         "Immune cells" = "immune",
                #                         "Myeloid cells" = "mye")),
                
                # Data level selection
                uiOutput("data_level_ui"),
                actionBttn(
                  "load_btn","Load Data"
                ),
                
                # Content to be hidden
                hidden(
                  div(id = "hidden_menu",

                      # Select to visualize genes or pathways
                      awesomeRadio("feature_type",
                                   "Select a type of features to visualize",
                                   choices = c("Genes", "Pathways"),
                                   selected = "Genes"
                                   
                      ),
                      ##### Branch based on selected feature type
                      # -------- Gene mode --------
                      conditionalPanel(
                        condition = "input.feature_type == 'Genes'",
                        # Load gene list
                        checkboxInput("use_textinput", "Use my own gene list", value = FALSE),
                        
                        conditionalPanel(
                          condition = "input.use_textinput == 0",
                          selectizeInput("gene_select","Select a gene", choices = NULL, multiple = TRUE)
                        ),
                        
                        conditionalPanel(
                          condition = "input.use_textinput == 1",
                          textAreaInput("gene_input", "Paste gene list:")
                        )
                      ),
                      # --------- Pathway mode ------
                      conditionalPanel(
                        condition = "input.feature_type == 'Pathways'",
                        # Load gene list
                        
                        # # Load pathway list according to pathway categories
                        # checkboxGroupButtons("pathway_cat","Pathway category",
                        #                      choices = c("Hallmark","Reactome")),
                        
                        ### Subset: rownames(fibs@assays$VAMcdf) to get a specific category of pathways
                        checkboxInput("use_textinput_VAM", "Add my own gene set", value = FALSE),
                        selectInput("pathway_select"," Select pathways", choices = NULL, multiple = TRUE),
                        conditionalPanel(
                          condition = "input.use_textinput_VAM ==1",
                          textInput("geneset_name","Name of your gene set"),
                          textAreaInput("VAM_geneset", "Please paste a list of genes for VAM"),
                          actionBttn("submit_geneset","Add geneset")
                        )
                      )
                      
                      )# End of div
                ),
                

                


                
              
                
              ),
              navset_tab(
                id = "explore_tabs",
                nav_panel(title = "Plots",
                          value = "plots",
                          layout_columns(
                            col_widths = c(6,6,12,12),
                            card_feature, # Feature plot of genes
                            card_umap, # UMAP of dataset
                            card_gene_plot # Violin plot / heatmap / dot plot of genes
                            
                          )
                ),
                nav_panel(title = "DEGs table",
                          value = "DEGs_table",
                          card_tabib_DEGs),
                nav_panel(title = "Metadata",
                          value = "meta_table",
                          uiOutput("meta_df"))
              ),
              
              
              
            ),
  ),
  
  
  
  nav_spacer(),
  nav_panel(
    icon(
      "calendar",
      align = "right"
    ) 
  )
  

)


