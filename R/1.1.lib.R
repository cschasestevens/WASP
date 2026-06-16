#' Universal Color Palette
#'
#' Combines the npg, aaas, and lancet ggsci palettes for use with datasets
#' containing up to 36 groups.
#'
#' @return Vector of colors to replace default discrete color scale.
#' @import ggsci
#' @examples
#'
#'  # col_univ()
#'
#' @export
col_univ <- function() {
  cat("Using default COMb color scheme (54 colors)", "\n")
  c(
    ggsci::pal_npg("nrc")(10),
    ggsci::pal_aaas("default")(10),
    ggsci::pal_lancet("lanonc")(9),
    ggsci::pal_frontiers("default")(7),
    ggsci::pal_nejm("default")(8),
    ggsci::pal_jco("default")(10)
  )
}

#' Gradient Color Palette
#'
#' Returns 12 colors from the viridis color palette.
#'
#' @param scm Gradient scheme to use
#' (1 = viridis, 2 = yellow/brown, 3 = blue/red-a, 4 = blue/red-b)
#' @return Vector of colors to replace default gradient color scale.
#' @importFrom viridis viridis
#' @importFrom RColorBrewer brewer.pal
#' @examples
#'
#' col_grad()
#'
#' @export
col_grad <- function(
  scm = 1
) {
  if (scm == 1) {
    cat("Using COMb gradient", "#1:", "viridis", "\n")
    c1 <- viridis::viridis(12)
  }
  if (scm == 2) {
    cat("Using COMb gradient", "#2:", "YlOrBr", "\n")
    c1 <- RColorBrewer::brewer.pal(name = "YlOrBr", n = 9)
  }
  if (scm == 3) {
    cat("Using COMb gradient", "#3:", "RWB heatmap B", "\n")
    c1 <- c("#2e86c1", "white", "#f5b7b1", "#e74c3c")
  }
  if (scm == 4) {
    cat("Using COMb gradient", "#4:", "Zero/positive scale", "\n")
    c1 <- c("lightblue", "red", "darkred")
  }
  if (scm == 5) {
    cat("Using COMb gradient", "#5:", "RWB heatmap A", "\n")
    c1 <- c("#2e86c1", "white", "#e74c3c")
  }
  if (scm == 6) {
    cat("Using COMb gradient", "#6:", "RWB heatmap C", "\n")
    c1 <- c("dodgerblue4", "#2e86c1", "white", "#e74c3c")
  }
  if (scm == 7) {
    cat("Using COMb gradient", "#7:", "RWB heatmap D", "\n")
    c1 <- c("dodgerblue4", "#2e86c1", "white", "#e74c3c", "darkred") # nolint
  }
  return(c1) # nolint
}

#' Generic Plot Theme
#'
#' General plotting theme.
#'
#' @param thm_type COMb theme type (either "default" or "extra_axis_detail").
#' @param txt_legtitle Legend title font size.
#' @param txt_legtext Legend text font size.
#' @param txt_plottitle Plot title font size.
#' @param txt_strip Strip font size.
#' @param txt_axes Axis title and text font size.
#' @param size_leg Legend key size (in cm).
#' @param leg Plot legend coordinates.
#'
#' @return ggplot2 theme parameters to replace default plot theme.
#' @import ggplot2
#' @examples
#'
#' # ms_theme()
#'
#' @export
sc_theme1 <- function(
  thm_type = "default",
  txt_legtitle = 14,
  txt_legtext = 10,
  txt_plottitle = 14,
  txt_strip = 12,
  txt_axes = 14,
  size_leg = 0.2,
  leg = c(0.95, 0.95)
) {
  cat("Adding default COMb plot theme", "\n")
  #---- Plot legend ----
  thm_leg_main <- ggplot2::theme(
    legend.title = ggplot2::element_text(
      size = txt_legtitle,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(size = txt_legtext),
    legend.key.size = ggplot2::unit(size_leg, "cm"),
    legend.key = ggplot2::element_blank(),
    legend.position.inside = leg
  )
  #---- All Plots ----
  thm_gen <- ggplot2::theme(
    # Plot title
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = txt_plottitle
    ),
    # Strip
    strip.background = ggplot2::element_rect(fill = "slategray2"),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = txt_strip
    ),
    # Margins
    plot.margin = ggplot2::unit(
      c(0.5, 0.25, 0.5, 0.25),
      "cm"
    ),
    # Axis titles and text
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = txt_axes,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = txt_axes
    ),
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = txt_axes
    ),
    axis.title.y = ggplot2::element_text(
      face = "bold",
      size = txt_axes
    )
  )
  #---- For default plots ----
  if (thm_type == "default") {
    thm_gen <- thm_gen +
      ggplot2::theme(
        # Panel
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey85"),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        # y-axis ticks
        axis.ticks.y = ggplot2::element_blank()
      )
  }
  if (thm_type == "extra_axis_detail") {
    thm_gen <- thm_gen +
      ggplot2::theme(
        # Panel
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey85"),
        panel.grid.major.x = ggplot2::element_line(colour = "grey85")
      )
  }
  # Output
  thm_out <- thm_gen +
    thm_leg_main
  return(thm_out) # nolint
}
