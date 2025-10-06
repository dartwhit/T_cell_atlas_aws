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
source("modules/scrna_seq_module.R")


# ########################################### Define UI ########################################
ui <- page_navbar(
  useShinyjs(), # Enable shinyjs

  tags$head(
    tags$script(HTML('\n      $(document).on("shiny:connected", function() {\n        function updateDimensions() {\n          Shiny.onInputChange("dimension", [window.innerWidth, window.innerHeight]);\n        }\n        updateDimensions();\n        $(window).resize(updateDimensions);\n\n        // Any element with id starting with "explore_" switches to Explore tab\n        $("[id^=\'explore_\']").on("click", function() {\n          $("a[data-value=\'Explore\']").tab("show");\n        });\n      });\n    '))
    ,
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
  #                
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
              scrna_seq_UI("scrna_seq_module")
            )
  )
)


