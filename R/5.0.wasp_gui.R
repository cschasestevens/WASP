#' Run WASP GUI
#'
#' Launches the WASP GUI. Useful for interactive analysis
#' of several commonly used WASP functions.
#'
#' @param ... Arguments passed to the WASP GUI.
#' @export
run_wasp_gui <- function(...) {
  app_dir <- system.file("shiny", "app", package = "WASP")
  if (app_dir == "") {
    stop("Could not find the app file. Try reinstalling the package...",
         call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
