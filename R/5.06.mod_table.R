#' Summary Table Module UI
#'
#' @param id A string. The module namespace id.
#' @export
mod_summary_table_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tableOutput(ns("table"))
  )
}

#' Summary Table Module Server
#'
#' @param id A string. The module namespace id.
#' @param data A reactive Seurat object from mod_filter_server.
#' @export
mod_summary_table_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {

    output$table <- renderTable({
      req(data())
      head(data()@meta.data, 10)
    })

  })
}
