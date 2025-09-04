library(shiny)
library(shinyWidgets)

# ---------- MODULE UI ----------
explore_sidebar_UI <- function(id) {
  ns <- NS(id)
  sidebar(
    title = "Select dataset and genes",
    position = "left",
    selectInput(ns("study"),"Select study to explore",
                choices = c(
                  "TMKMH" = "tmkmh",
                  "Tabib et al." = "tabib",
                  "Gur et al." = "gur",
                  "Ma et al." = "ma",
                  "Khanna et al" = "khanna"
                )),

    
    conditionalPanel(
      condition="input.data_level == 'full'",
      ns = ns,
      switchInput(ns("anno"),
                  "Original seurat clusters",
                  value = FALSE)
    ),
    
    uiOutput(ns("data_level_ui")),
    actionBttn(
      ns("load_btn"),"Load Data"
    ),
    
    hidden(
      div(id = ns("hidden_menu"),

          awesomeRadio(ns("feature_type"),
                       "Select a type of features to visualize",
                       choices = c("Genes", "Pathways"),
                       selected = "Genes"
                       
          ),
          
          conditionalPanel(
            condition = "input.feature_type == 'Genes'",
            ns = ns,
            checkboxInput(ns("use_textinput"), "Use my own gene list", value = FALSE),
            
            conditionalPanel(
              condition = "input.use_textinput == 0",
              ns = ns,
              selectizeInput(ns("gene_select"),"Select a gene", choices = NULL, multiple = TRUE)
            ),
            
            conditionalPanel(
              condition = "input.use_textinput == 1",
              ns = ns,
              textAreaInput(ns("gene_input"), "Paste gene list:")
            )
          ),
          
          conditionalPanel(
            condition = "input.feature_type == 'Pathways'",
            ns = ns,
            checkboxInput(ns("use_textinput_VAM"), "Add my own gene set", value = FALSE),
            selectInput(ns("pathway_select")," Select pathways", choices = NULL, multiple = TRUE),
            conditionalPanel(
              condition = "input.use_textinput_VAM ==1",
              ns = ns,
              textInput(ns("geneset_name"),"Name of your gene set"),
              textAreaInput(ns("VAM_geneset"), "Please paste a list of genes for VAM"),
              actionBttn(ns("submit_geneset"),"Add geneset")
            )
          )
          
      )
    )
  )
}

# ---------- MODULE SERVER ----------
explore_sidebar_server <- function(id, selected_study_from_gallery) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$data_level_ui <- renderUI({
      selectInput(ns("data_level"),
                  "Select data to visualize",
                  choices = data_level_choices[[input$study]])
    })

    observeEvent(selected_study_from_gallery(), {
      updateSelectInput(session, "study", selected = selected_study_from_gallery())
    })

    return(
      list(
        study = reactive(input$study),
        data_level = reactive(input$data_level),
        anno = reactive(input$anno),
        load_btn = reactive(input$load_btn),
        feature_type = reactive(input$feature_type),
        use_textinput = reactive(input$use_textinput),
        gene_select = reactive(input$gene_select),
        gene_input = reactive(input$gene_input),
        use_textinput_VAM = reactive(input$use_textinput_VAM),
        pathway_select = reactive(input$pathway_select),
        geneset_name = reactive(input$geneset_name),
        VAM_geneset = reactive(input$VAM_geneset),
        submit_geneset = reactive(input$submit_geneset)
      )
    )
  })
}
