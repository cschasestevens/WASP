#' Data Import Module UI
#'
#' @param id A string. The module namespace id.
#' @export
mod_data_import_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(
      ns("file"),
      label       = "Upload a Seurat object (.rds):",
      accept      = ".rds",
      placeholder = "No file selected"
    ),
    # Rendered dynamically after a file is successfully loaded
    uiOutput(ns("import_status"))
  )
}

#' Data Import Module Server
#'
#' @param id A string. The module namespace id.
#' @return A reactive containing the imported Seurat object.
#' @export
mod_data_import_server <- function(id, session = shiny::getDefaultReactiveDomain()) {
  moduleServer(id, function(input, output, session) {

    imported_data <- reactive({
      # req() stops execution silently until a file is uploaded
      req(input$file)

      # Validate the file extension
      ext <- tools::file_ext(input$file$name)
      validate(
        need(ext == "rds", "Please upload a .rds file.")
      )

      # Attempt to read the file
      obj <- tryCatch(
        readRDS(input$file$datapath),
        error = function(e) {
          validate(need(FALSE, paste("Could not read file:", e$message)))
        }
      )

      # Validate that the object is a Seurat object
      validate(
        need(
          inherits(obj, "Seurat"),
          "The uploaded .rds file does not contain a Seurat object."
        )
      )

      obj
    })

    # Render a status message showing cell and feature counts on success
    output$import_status <- renderUI({
      req(imported_data())
      obj <- imported_data()
      tags$p(
        style = "color: green;",
        sprintf(
          "\u2713 Loaded Seurat object: %d cells, %d features.",
          ncol(obj),
          nrow(obj)
        )
      )
    })

    return(imported_data)

  })
}
