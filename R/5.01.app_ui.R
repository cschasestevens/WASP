#' WASP GUI - UI
#'
#' UI for WASP GUI app.
#'
#' @export
app_ui <- function() {
  fluidPage(
    titlePanel("My Modular App"),
    sidebarLayout(
      sidebarPanel(
        mod_data_import_ui("import"),  # data import module
        hr(),
        mod_filter_ui("filter")        # filter module
      ),
      mainPanel(
        mod_plot_ui("plot"),            # plot module
        mod_summary_table_ui("table")  # summary table module
      )
    )
  )
}
