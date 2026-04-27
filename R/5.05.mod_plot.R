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

    output$scatter <- renderPlot(
      expr = {
        req(data())
        WASP::sc_umap(
          so = data(),
          md_var = "CellType",
          slot1 = "wnn.umap"
        )
      },
      # width and height accept functions that are called reactively.
      # "auto" tells Shiny to use the CSS width of the output container.
      width  = "auto",
      height = function() {
        # Build the clientData key for this output's width.
        # Inside a module, session$ns() adds the module prefix to the id,
        # e.g. "plot-plot". The [[ ]] accessor is needed because the key
        # contains a hyphen.
        width_key <- paste0("output_", session$ns("plot"), "_width")
        w <- session$clientData[[width_key]]
        # Return a fallback height if the width is not yet available
        if (is.null(w) || w == 0) return(400)
        round(w * 1)
      },
      # res controls resolution (pixels per inch).
      # 96 matches standard screen resolution.
      res = 96
    )

  })
}
