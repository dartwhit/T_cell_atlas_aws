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
  ),
  
  # Spatial page with bslib page_fluid
  nav_panel("Spatial",
            page_fluid(
              layout_columns(
                col_widths = c(3, 9),
                card(
                  card_header("Spatial Analysis Controls"),
                  card_body(
                    spatial_sidebar_UI("spatial_sidebar")
                  )
                ),
                card(
                  full_screen = TRUE,
                  card_body(
                    class = "p-0",
                    spatial_UI("sp1")
                  )
                )
              )
            )
  )
)


