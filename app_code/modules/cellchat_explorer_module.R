library(shiny)
library(bslib)
library(DT)

# ---------- MODULE UI ----------
cellchat_explorer_UI <- function(id, choices) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      title = "CellChat options",
      selectInput(ns("study"), "Select study", choices = choices),
      radioButtons(
        ns("measure"), "Aggregate measure",
        choices = c("Strength" = "weight", "Count" = "count"),
        selected = "weight", inline = TRUE
      ),
      selectInput(ns("condition"), "Condition for drill-down", choices = NULL),
      tags$hr(),
      tags$strong("Filters"),
      selectizeInput(
        ns("sources"), "Source cell type(s)", choices = NULL, multiple = TRUE,
        options = list(placeholder = "All")
      ),
      selectizeInput(
        ns("targets"), "Target cell type(s)", choices = NULL, multiple = TRUE,
        options = list(placeholder = "All")
      ),
      selectizeInput(
        ns("pathways"), "Pathway(s)", choices = NULL, multiple = TRUE,
        options = list(placeholder = "All")
      ),
      textInput(ns("ligand"), "Ligand contains", ""),
      textInput(ns("receptor"), "Receptor contains", ""),
      sliderInput(ns("prob_min"), "Minimum probability", 0, 1, 0, step = 0.01),
      sliderInput(ns("pval_max"), "Maximum p-value", 0, 1, 1, step = 0.01)
    ),
    tagList(
      uiOutput(ns("status")),
      navset_tab(
        nav_panel(
          "Comparison overview",
          p("Aggregated interactions per condition, differential interactions, and total comparison."),
          plotOutput(ns("circle"), height = "460px"),
          layout_columns(
            col_widths = c(6, 6),
            card(card_header("Differential interactions"), plotOutput(ns("diff"), height = "400px")),
            card(card_header("Total interactions"), plotOutput(ns("compare"), height = "400px"))
          )
        ),
        nav_panel(
          "Heatmaps",
          radioButtons(
            ns("role_pattern"), "Signaling role pattern",
            choices = c("Outgoing" = "outgoing", "Incoming" = "incoming", "Overall" = "all"),
            selected = "outgoing", inline = TRUE
          ),
          plotOutput(ns("heatmap_diff"), height = "460px"),
          plotOutput(ns("heatmap_role"), height = "560px")
        ),
        nav_panel(
          "L-R bubble",
          p("Ligand-receptor pairs across selected source and target cell types."),
          plotOutput(ns("bubble"), height = "720px")
        ),
        nav_panel(
          "Signaling roles / rankNet",
          checkboxInput(ns("rank_stacked"), "Stacked (absolute) rankNet", value = TRUE),
          plotOutput(ns("ranknet"), height = "640px"),
          plotOutput(ns("role_scatter"), height = "480px")
        ),
        nav_panel(
          "Communication table",
          radioButtons(
            ns("level"), "Resolution",
            choices = c("Ligand-receptor" = "LR", "Pathway" = "pathway"),
            selected = "LR", inline = TRUE
          ),
          downloadButton(ns("download"), "Download filtered CSV"),
          br(), br(),
          DTOutput(ns("table"))
        )
      )
    )
  )
}

# ---------- MODULE SERVER ----------
cellchat_explorer_server <- function(id, dataset_configs, gallery_selection = NULL) {
  moduleServer(id, function(input, output, session) {
    selected_study <- reactive({ input$study })

    if (!is.null(gallery_selection)) {
      observeEvent(gallery_selection(), {
        selection <- gallery_selection()
        if (is.list(selection) && identical(selection$view, "cellchat") &&
            selection$id %in% names(dataset_configs)) {
          updateSelectInput(session, "study", selected = selection$id)
        }
      }, ignoreInit = TRUE)
    }

    cellchat_data <- reactive({
      study <- selected_study()
      req(study, study %in% names(dataset_configs))
      cellchat_load_objects(study, dataset_configs[[study]])
    })

    output$status <- renderUI({
      data <- cellchat_data()
      if (!is.null(data$error)) {
        div(class = "alert alert-danger", role = "alert", data$error)
      } else {
        NULL
      }
    })

    observeEvent(cellchat_data(), {
      data <- cellchat_data()
      if (!is.null(data$error)) {
        updateSelectInput(session, "condition", choices = character())
        updateSelectizeInput(session, "sources", choices = character(), selected = character())
        updateSelectizeInput(session, "targets", choices = character(), selected = character())
        updateSelectizeInput(session, "pathways", choices = character(), selected = character())
        return()
      }

      conditions <- names(data$conditions)
      updateSelectInput(
        session, "condition", choices = conditions,
        selected = conditions[[1]] %||% NULL
      )
      celltypes <- cellchat_celltypes(data$merged)
      pathways <- cellchat_pathways(data$merged)
      updateSelectizeInput(session, "sources", choices = celltypes, selected = character())
      updateSelectizeInput(session, "targets", choices = celltypes, selected = character())
      updateSelectizeInput(session, "pathways", choices = pathways, selected = character())
    }, ignoreInit = FALSE)

    selected_condition_object <- reactive({
      data <- cellchat_data()
      req(is.null(data$error), input$condition)
      data$conditions[[input$condition]]
    })

    draw_plot <- function(plot_function) {
      data <- cellchat_data()
      if (!is.null(data$error)) {
        cellchat_draw_plot_error(data$error)
        return(invisible(NULL))
      }
      tryCatch(
        plot_function(data),
        error = function(error) {
          cellchat_draw_plot_error(paste("CellChat plot could not be generated:",
                                         conditionMessage(error)))
        }
      )
    }

    output$circle <- renderPlot({
      draw_plot(function(data) {
        cellchat_plot_circle_comparison(data$merged, input$measure)
      })
    })
    output$diff <- renderPlot({
      draw_plot(function(data) {
        cellchat_plot_diff_interaction(data$merged, input$measure)
      })
    })
    output$compare <- renderPlot({
      draw_plot(function(data) {
        cellchat_plot_compare_interactions(data$merged, input$measure)
      })
    })
    output$heatmap_diff <- renderPlot({
      draw_plot(function(data) {
        cellchat_plot_heatmap_diff(data$merged, input$measure)
      })
    })
    output$heatmap_role <- renderPlot({
      draw_plot(function(data) {
        object <- selected_condition_object()
        req(object)
        cellchat_plot_signaling_role_heatmap(object, input$role_pattern, input$pathways)
      })
    })
    output$bubble <- renderPlot({
      draw_plot(function(data) {
        cellchat_plot_bubble(data$merged, input$sources, input$targets, input$pathways)
      })
    })
    output$ranknet <- renderPlot({
      draw_plot(function(data) {
        cellchat_plot_ranknet(data$merged, input$rank_stacked)
      })
    })
    output$role_scatter <- renderPlot({
      draw_plot(function(data) {
        object <- selected_condition_object()
        req(object)
        cellchat_plot_signaling_role_scatter(object, input$pathways)
      })
    })

    communication_table <- reactive({
      data <- cellchat_data()
      if (!is.null(data$error)) return(tibble::tibble())
      cellchat_extract_comm(data$merged, input$level)
    })

    filtered_table <- reactive({
      cellchat_apply_comm_filters(
        communication_table(),
        sources = input$sources,
        targets = input$targets,
        ligand = input$ligand,
        receptor = input$receptor,
        pathway = input$pathways,
        prob_min = input$prob_min,
        pval_max = input$pval_max
      )
    })

    output$table <- renderDT({
      datatable(
        filtered_table(), filter = "top",
        options = list(pageLength = 25, scrollX = TRUE)
      )
    })

    output$download <- downloadHandler(
      filename = function() {
        paste0("cellchat_", selected_study(), "_", input$level, "_filtered.csv")
      },
      content = function(file) {
        utils::write.csv(filtered_table(), file, row.names = FALSE)
      }
    )
  })
}
