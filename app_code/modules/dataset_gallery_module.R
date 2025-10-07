library(shiny)
library(shinyWidgets)

# ---------- MODULE UI ----------
dataset_gallery_UI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      textInput(ns("search"), "Search"),
      radioGroupButtons(
        ns("assay"), 
        label = "Assay", 
        choices = c("All", "scRNA-seq", "Spatial"),
        selected = "All"
      )
    ),
    uiOutput(ns("gallery"))
  )
}

# ---------- MODULE SERVER ----------
dataset_gallery_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    filtered_datasets <- reactive({
      data <- dataset_meta
      
      # Filter by assay type
      if (input$assay == "scRNA-seq") {
        data <- data[data$has_scrna, ]
      } else if (input$assay == "Spatial") {
        data <- data[data$has_spatial, ]
      }

      # Filter by search term
      if (nzchar(input$search)) {
        term <- tolower(input$search)
        data <- data[
          grepl(term, tolower(data$name)) |
            grepl(term, tolower(data$tags)) |
            grepl(term, tolower(data$desc)),
        ]
      }
      data
    })

    output$gallery <- renderUI({
      data <- filtered_datasets()
      if (nrow(data) == 0) {
        return(p("No datasets found."))
      }
      cards <- lapply(seq_len(nrow(data)), function(i) {
        row <- data[i, ]
        
        # Create badges for available assay types
        badges <- list()
        if (row$has_scrna) {
          badges <- append(badges, list(tags$span(class = "badge", "scRNA-seq")))
        }
        if (row$has_spatial) {
          badges <- append(badges, list(tags$span(class = "badge", "Spatial")))
        }

        # Create action buttons based on available data
        action_buttons <- list()
        has_both <- row$has_scrna && row$has_spatial
        
        if (has_both) {
          action_buttons <- list(
            actionButton(ns(paste0("explore_scrna_", row$id)), "View scRNA-seq"),
            actionButton(ns(paste0("explore_spatial_", row$id)), "View Spatial")
          )
        } else if (row$has_scrna) {
          action_buttons <- list(actionButton(ns(paste0("explore_scrna_", row$id)), "View scRNA-seq"))
        } else if (row$has_spatial) {
          action_buttons <- list(actionButton(ns(paste0("explore_spatial_", row$id)), "View Spatial"))
        }

        tags$div(
          class = "dataset-card",
          tags$img(src = row$image, class = "dataset-img", alt = row$name),
          h4(row$name),
          tags$div(class = "badges", badges),
          tags$p(sprintf("%s cells", row$n_cells)),
          tags$div(class = "card-actions", action_buttons)
        )
      })
      div(class = "dataset-gallery", cards)
    })

    # Create a reactive value to store the selected study and view
    selected_study <- reactiveVal(NULL)

    # Observe button clicks for all datasets
    lapply(dataset_meta$id, function(id) {
      # Observer for scRNA-seq button
      observeEvent(input[[paste0("explore_scrna_", id)]], {
        # Add timestamp to ensure value always changes
        selected_study(list(id = id, view = "scrna", timestamp = Sys.time()))
      })
      # Observer for Spatial button
      observeEvent(input[[paste0("explore_spatial_", id)]], {
        # Add timestamp to ensure value always changes
        selected_study(list(id = id, view = "spatial", timestamp = Sys.time()))
      })
    })

    # Return the reactive value
    return(selected_study)
  })
}
