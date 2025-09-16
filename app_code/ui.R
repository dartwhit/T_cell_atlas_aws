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
source("modules/dataset_gallery_module.R")
source("modules/explore_sidebar_module.R")


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
                            switchInput("plot_type", "Plot Type", value = TRUE, onLabel = "Vln", offLabel = "Box"),
                            
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
    dataset_gallery_UI("gallery_module")
  ),
  
  
  
  # Explore page
    nav_panel("Explore",
            layout_sidebar(
              sidebar = explore_sidebar_UI("explore_sidebar_module"),
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
    layout_sidebar(
      sidebar = sidebar(
        title = "Spatial data options",
        selectInput("spatial_study_selector", "Select Study", choices = NULL),
        selectInput(NS("sp1", "group_by"), 
                    "Group by (metadata)", 
                    choices = c("Seurat cluster" = "cluster",
                    "Key regions" =  "Key_Regions"),
                  selected = "cluster"),
        selectizeInput(NS("sp1", "feature"), "Feature (Gene/pathway)",multiple = TRUE,
        choices = NULL, options = list(placeholder = "Type to search…")),
        checkboxGroupInput(NS("sp1", "samples"),
                      "Select samples to view"
                    )
      ),
      spatial_UI("sp1")
    )
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


