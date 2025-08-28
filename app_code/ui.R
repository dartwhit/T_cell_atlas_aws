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

## Source UI modules
source("modules/spatial_unit.R")


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
    tags$script(HTML('
      $(document).on("shiny:connected", function() {
        function updateDimensions() {
          Shiny.onInputChange("dimension", [window.innerWidth, window.innerHeight]);
        }
        updateDimensions();
        $(window).resize(updateDimensions);

        // Any element with id starting with "explore_" switches to Explore tab
        $("[id^=\'explore_\']").on("click", function() {
          $("a[data-value=\'Explore\']").tab("show");
        });
      });
    ')),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")


  ),
  theme = bs_theme(bootswatch = "flatly"),
  # Page title
  title = "SSc cell atlas",
  id = "nav_page",
  # nav_panel("Datasets",
  #           card(full_screen = TRUE,
  #                card_header("TMKMH integrated dataset"),
  #                actionButton("explore_tmkmh","Explore"),
  #                imageOutput("TMKMH_img")
                 
  #                )
  #           ,
  #           layout_column_wrap(
  #             width = 1/4,

  #             #Tabib
  #             card(full_screen = TRUE,
  #                  card_header("Tabib et al. 2021"),
  #                  actionButton("explore_Tabib","Explore"),
  #                  imageOutput("Tabib_img"),
  #                  a(href = "https://www.nature.com/articles/s41467-021-24607-6",
  #                    p(strong("Tabib, T., Huang, M., Morse, N., Papazoglou, A., Behera, R., Jia, M., Bulik, M., Monier, D. E., Benos, P. v., Chen, W., Domsic, R., & Lafyatis, R. (2021)."),"Myofibroblast transcriptome indicates SFRP2hi fibroblast progenitors in systemic sclerosis skin. Nature Communications, 12(1), 4384."),
  #                    style = "color:grey", target="_blank")
  #             ),
  #             # Gur
  #             card(full_screen = TRUE,
  #                  card_header("Gur et al. 2022"),
  #                  actionButton("explore_Gur","Explore"),
  #                  imageOutput("Gur_img")),

  #             # Ma
  #             card(full_screen = TRUE,
  #                  card_header("Ma et al. 2024"),
  #                  actionButton("explore_Ma","Explore"),
  #                  imageOutput("Ma_img")),

  #             # Clark
  #             card(full_screen = TRUE,
  #                  card_header("Khanna et al. 2022"),
  #                  actionButton("explore_Khanna","Explore"),
  #                  imageOutput("Khanna_img")),

  #           )
  # ),
    nav_panel(
    "Datasets",
    layout_sidebar(
      sidebar = sidebar(
        textInput("search", "Search"),
        selectInput("assay", "Assay", choices = c("All", unique(dataset_meta$assay)))
        # sliderInput(
        #   "cell_count",
        #   "Cell count",
        #   min = 0,
        #   max = max(dataset_meta$n_cells),
        #   value = c(0, max(dataset_meta$n_cells))
        # )
      ),
      uiOutput("gallery")
    )
  ),
  
  
  
  # Explore page
    nav_panel("Explore",
            layout_sidebar(
              sidebar = sidebar(
                title = "Select dataset and genes",
                position = "left",
                selectInput("study","Select study to explore",
                            choices = c(
                              "TMKMH" = "tmkmh",
                              "Tabib et al." = "tabib",
                              "Gur et al." = "gur",
                              "Ma et al." = "ma",
                              "Khanna et al" = "khanna"
                            )),

                
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
  
  nav_panel("Spatial data explorer",
value = "spatial",
bslib::page_fluid(
  # theme = bs_theme(bootswatch = "flatly"),
  h3("Spatial module test"),
  spatial_UI("sp1")
)
  ),
  
  nav_spacer(),
  nav_panel(
    icon(
      "calendar",
      align = "right"
    ) 
  )
  

)


