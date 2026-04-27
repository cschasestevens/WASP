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
    uiOutput(ns("import_status"))
  )
}

#' Data Import Module Server
#'
#' @param id A string. The module namespace id.
#' @return A reactive containing the imported Seurat object.
#' @export
mod_data_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    imported_data <- reactive({
      req(input$file)

      # Validate the file extension
      ext <- tools::file_ext(input$file$name)
      validate(
        need(ext == "rds", "Please upload a .rds file.")
      )

      # Attempt to read the file
      obj <- tryCatch(
        {
          d1 <- readRDS(input$file$datapath)
          return(d1)
        },
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

      # Return object for filtering
      return(obj)

    })

    output$import_status <- shiny::renderUI({
      shiny::req(imported_data())
      shiny::tags$p(
        style = "color: green;",
        sprintf(
          "\u2713 Loaded Seurat object containing: %d cells, %d assays.",
          ncol(imported_data()),
          length(SeuratObject::Assays(imported_data()))
        )
      )
    })

    return(imported_data)

  })
}
