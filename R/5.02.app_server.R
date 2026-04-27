#' Run WASP GUI
#'
#' Main server function for WASP GUI.
#'
#' @param input Shiny app input.
#' @param output Shiny app output.
#' @param session Shiny app session information
#' @export
app_server <- function(input, output) {
  # ---- 1. Import RDS File ----
  imported_data <- mod_data_import_server("import")
  # 2. Filter — takes the Seurat object, returns a filtered reactive Seurat object
  filtered_data <- mod_filter_server("filter", data = imported_data)

  # 3. Downstream modules consume the filtered Seurat object
  mod_plot_server("plot", data = filtered_data)
  mod_summary_table_server("table", data = filtered_data)

}
