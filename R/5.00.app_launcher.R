#' Launch the Shiny application
#'
#' @param ... Arguments passed to \code{\link[shiny]{runApp}}, such as
#'   \code{port} or \code{launch.browser}
#'
#' @export
run_app <- function(...) {

  # Increase upload limit to 2 GB — adjust as needed
  options(shiny.maxRequestSize = 10 * 1024^3)

  app_dir <- system.file("shiny", "app.r", package = "WASP")

  # system.file() returns "" if the directory cannot be found
  if (app_dir == "") {
    stop("Could not find the app directory. Try reinstalling the package.",
         call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
