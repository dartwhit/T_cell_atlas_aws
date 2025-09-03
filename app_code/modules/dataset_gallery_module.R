library(shiny)

# ---------- MODULE UI ----------
dataset_gallery_UI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      textInput(ns("search"), "Search"),
      selectInput(ns("assay"), "Assay", choices = c("All", unique(dataset_meta$assay)))
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
      if (input$assay != "All") {
        data <- data[data$assay == input$assay, ]
      }
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
        badges <- lapply(c(row$assay, strsplit(row$tags, ",")[[1]]), function(tg) {
          tags$span(class = "badge", tg)
        })
        tags$div(
          class = "dataset-card",
          tags$img(src = row$image, class = "dataset-img", alt = row$name),
          h4(row$name),
          badges,
          tags$p(sprintf("%s cells", row$n_cells)),
          actionButton(ns(paste0("explore_", row$id)), "Analyze")
        )
      })
      div(class = "dataset-gallery", cards)
    })

    # Create a reactive value to store the selected study
    selected_study <- reactiveVal(NULL)

    # Observe the button clicks
    lapply(dataset_meta$id, function(id) {
      observeEvent(input[[paste0("explore_", id)]], {
        selected_study(id)
      })
    })

    # Return the reactive value
    return(selected_study)
  })
}
