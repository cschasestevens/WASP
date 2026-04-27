#' Launch the Shiny application
#'
#' @param ... Arguments passed to \code{\link[shiny]{runApp}}, such as
#'   \code{port} or \code{launch.browser}
#'
#' @export
run_app <- function(...) {
  # Set max size for RDS object
  options(shiny.maxRequestSize = 10 * 1024^3)
  # Select app file to run
  app_dir <- system.file("shiny", "app.r", package = "WASP")
  # Return error if app files are broken
  if (app_dir == "") {
    stop("Could not find the app directory. Try reinstalling the package.",
         call. = FALSE)
  }
  # run app file
  shiny::runApp(app_dir, ...)
}
