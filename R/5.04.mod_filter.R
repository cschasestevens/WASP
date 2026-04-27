#' Filter Module UI
#'
#' @param id A string. The module namespace id.
#' @export
mod_filter_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Both selectors are rendered dynamically based on the uploaded data
    uiOutput(ns("column_select")),
    uiOutput(ns("value_select"))
  )
}

#' Filter Module Server
#'
#' @param id A string. The module namespace id.
#' @param data A reactive Seurat object from mod_data_import_server.
#' @return A reactive filtered Seurat object.
#' @export
mod_filter_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {

    # Dynamically build column selector from Seurat metadata
    output$column_select <- renderUI({
      req(data())
      ns <- session$ns
      selectInput(
        ns("column"),
        label   = "Filter by metadata column:",
        choices = names(data()@meta.data)
      )
    })

    # Dynamically build value selector based on the chosen column
    output$value_select <- renderUI({
      req(data(), input$column)
      ns <- session$ns
      selectInput(
        ns("value"),
        label    = "Select value:",
        choices  = c("All", unique(data()@meta.data[[input$column]])),
        selected = "All"
      )
    })

    # Return the filtered Seurat object
    reactive({
      req(data())
      if (is.null(input$value) || input$value == "All") {
        data()
      } else {
        subset(data(), subset = .data[[input$column]] == input$value)
      }
    })

  })
}
