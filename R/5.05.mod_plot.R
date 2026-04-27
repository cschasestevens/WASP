#' Plot Module UI
#'
#' @param id A string. The module namespace id.
#' @export
mod_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(
      ns("scatter")
    )
  )
}

#' Plot Module Server
#'
#' @param id A string. The module namespace id.
#' @param data A reactive Seurat object from mod_filter_server.
#' @export
mod_plot_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {

    output$scatter <- renderPlot({
        req(data())
        Seurat::DimPlot(data())
      })

  })
}
