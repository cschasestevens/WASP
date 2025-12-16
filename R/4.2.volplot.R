#' Volcano Plot
#'
#' Generates a volcano plot from a DEG results object or data frame for
#' a chosen comparison and/or cell type.
#'
#' @param l_deg A list of DGEA results returned by sc.DGEA().
#' @param ct A character string pattern for matching to specific cell types.
#' @param comp_name Comparison name, provided as a character string.
#' @param comp_filt Filter data by comparison name.
#' @param gene_name Column name containing gene or motif names.
#' @param filt_lab (optional) Select labels to include on plot.
#' @param diff_col Column name containing effect size.
#' @param p_col Column name containing adjust p-values.
#' @param p_cut Numeric value for p-value cutoff to indicate significance.
#' @param f_cut Fold change cutoff indicating high effect size.
#' @param f_lim Numeric value for fold change limits on the x-axis.
#' @param y_limit Upper-bound -log10 p-value to display on the y-axis.
#' @param x_title X-axis legend title, provided as a character string.
#' @return A volcano plot for the chosen cell type and treatment comparison.
#' @examples
#'
#' # p_vol <- sc_volcano(
#' #   l_deg = diff.output.activity,
#' #   ct = "1.Secretory",
#' #   comp_name = "KO vs. Control",
#' #   gene_name = "gene",
#' #   filt_lab = grepl("NKX|FOXA", diff.output.activity[["gene"]]),
#' #   diff_col = "avg_diff",
#' #   p_col = "p_val_adj",
#' #   p_cut = 0.005,
#' #   f_cut = 0.25,
#' #   f_lim = 6,
#' #   y_limit = 100,
#' #   x_title = "x-axis title"
#' # )
#'
#' @export
sc_volcano <- function( # nolint
  l_deg,
  ct,
  comp_name = NULL,
  comp_filt = FALSE,
  gene_name = "GENE",
  filt_lab = NULL,
  diff_col = "log2FC",
  p_col = "H.qval",
  p_cut = 0.05,
  f_cut = 0.25,
  f_lim = 5,
  y_limit = 50,
  x_title
) {
  if(is.null(filt_lab) == TRUE) { # nolint
    if(missing(ct) == TRUE && class(l_deg) == "data.frame") { # nolint
      ld <- l_deg
      # Plot function
      v <- EnhancedVolcano::EnhancedVolcano(
        ld,
        lab = ld[[gene_name]],
        subtitle = ggplot2::element_blank(),
        caption = ggplot2::element_blank(),
        x = diff_col,
        y = p_col,
        pCutoff = p_cut,
        FCcutoff = f_cut,
        cutoffLineType = "twodash",
        legendLabels = c("NS", "Fold Change",
                         "p-value", "FC+P"),
        legendLabSize = 12,
        labFace = "bold",
        col = ggsci::pal_npg("nrc")(10)[c(4, 3, 5, 8)],
        colAlpha = 0.6,
        legendIconSize = 4,
        pointSize = 2,
        border = "full",
        borderWidth = 1.5,
        legendPosition = "'right'",
        labSize = 4,
        drawConnectors = TRUE,
        typeConnectors = "open",
        widthConnectors = 0.35,
        colConnectors = "grey75",
        max.overlaps = 18,
        min.segment.length = ggplot2::unit(
          1,
          "mm"
        )
      ) +
        sc_theme1() + # nolint
        ggplot2::labs(
          color = "Key",
          x = x_title,
          y = paste("-log10 P-value (FDR)")
        ) +
        ggplot2::ggtitle(
          paste(
            comp_name,
            sep = " "
          )
        ) +
        ggplot2::theme(
          plot.margin = ggplot2::unit(c(.2, .2, .2, .2), "cm")
        ) +
        ggplot2::coord_cartesian(xlim = c(-f_lim, f_lim),
                                 ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(
          breaks = seq(-f_lim, f_lim, (f_lim / 4))
        )
    }
    if(missing(ct) == FALSE && class(l_deg) == "data.frame") { # nolint
      # Subset Input data
      ld <- l_deg
      if (comp_filt == TRUE) {
        ld <- ld[ld[["Comparison"]] == comp_name, ]
      }
      ld <- ld[ld[["CellType"]] == ct, ]
      # Plot function
      v <- EnhancedVolcano::EnhancedVolcano(
        ld,
        lab = ld[[gene_name]],
        title = ggplot2::element_blank(),
        subtitle = ggplot2::element_blank(),
        caption = ggplot2::element_blank(),
        x = diff_col,
        y = p_col,
        pCutoff = p_cut,
        FCcutoff = f_cut,
        cutoffLineType = "twodash",
        legendLabels = c("NS", "Fold Change",
                         "p-value", "FC+P"),
        legendLabSize = 12,
        labFace = "bold",
        col = ggsci::pal_npg("nrc")(10)[c(4, 3, 5, 8)],
        colAlpha = 0.7,
        legendIconSize = 4,
        pointSize = 2,
        border = "full",
        borderWidth = 1.5,
        legendPosition = "'right'",
        labSize = 4,
        drawConnectors = TRUE,
        typeConnectors = "open",
        min.segment.length = ggplot2::unit(
          1,
          "mm"
        )
      ) +
        sc_theme1() + # nolint
        ggplot2::labs(color = "Key") +
        ggplot2::labs(
          color = "Key",
          x = x_title,
          y = paste("-log10 P-value (FDR)")
        ) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(.2, .2, .2, .2), "cm")) +
        ggplot2::coord_cartesian(xlim = c(-f_lim, f_lim),
                                 ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(breaks = seq(-f_lim, f_lim, (f_lim / 4)))
    }
    if(missing(ct) == FALSE && class(l_deg) == "list") { # nolint
      # Subset Input data
      ld <- l_deg[[1]]
      if (comp_filt == TRUE) {
        ld <- ld[ld[["Comparison"]] == comp_name, ]
      }
      ld <- ld[ld[["CellType"]] == ct, ]
      # Plot function
      v <- EnhancedVolcano::EnhancedVolcano(
        ld,
        lab = ld[[gene_name]],
        title = ggplot2::element_blank(),
        subtitle = ggplot2::element_blank(),
        caption = ggplot2::element_blank(),
        x = diff_col,
        y = p_col,
        pCutoff = p_cut,
        FCcutoff = f_cut,
        cutoffLineType = "twodash",
        legendLabels = c("NS", "Fold Change",
                         "p-value", "FC+P"),
        legendLabSize = 12,
        labFace = "bold",
        col = ggsci::pal_npg("nrc")(10)[c(4, 3, 5, 8)],
        colAlpha = 0.7,
        legendIconSize = 4,
        pointSize = 2,
        border = "full",
        borderWidth = 1.5,
        legendPosition = "'right'",
        labSize = 4,
        drawConnectors = TRUE,
        typeConnectors = "open",
        max.overlaps = 18,
        min.segment.length = ggplot2::unit(
          1,
          "mm"
        )
      ) +
        sc_theme1() + # nolint
        ggplot2::labs(color = "Key") +
        ggplot2::labs(
          color = "Key",
          x = x_title,
          y = paste("-log10 P-value (FDR)")
        ) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(.2, .2, .2, .2), "cm")) +
        ggplot2::coord_cartesian(xlim = c(-f_lim, f_lim),
                                 ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(breaks = seq(-f_lim, f_lim, (f_lim / 4)))
    }
  }
  if(is.null(filt_lab) == FALSE) { # nolint
    if(missing(ct) == TRUE && class(l_deg) == "data.frame") { # nolint
      ld <- l_deg
      # Plot function
      v <- EnhancedVolcano::EnhancedVolcano(
        ld,
        lab = ld[[gene_name]],
        selectLab = filt_lab,
        subtitle = ggplot2::element_blank(),
        caption = ggplot2::element_blank(),
        x = diff_col,
        y = p_col,
        pCutoff = p_cut,
        FCcutoff = f_cut,
        cutoffLineType = "twodash",
        legendLabels = c("NS", "Fold Change",
                         "p-value", "FC+P"),
        legendLabSize = 12,
        labFace = "bold",
        col = ggsci::pal_npg("nrc")(10)[c(4, 3, 5, 8)],
        colAlpha = 0.6,
        legendIconSize = 4,
        pointSize = 2,
        border = "full",
        borderWidth = 1.5,
        legendPosition = "'right'",
        labSize = 4,
        drawConnectors = TRUE,
        typeConnectors = "open",
        widthConnectors = 0.35,
        colConnectors = "grey75",
        max.overlaps = 18,
        min.segment.length = ggplot2::unit(
          1,
          "mm"
        )
      ) +
        sc_theme1() + # nolint
        ggplot2::labs(
          color = "Key",
          x = x_title,
          y = paste("-log10 P-value (FDR)")
        ) +
        ggplot2::ggtitle(
          paste(
            comp_name,
            sep = " "
          )
        ) +
        ggplot2::theme(
          plot.margin = ggplot2::unit(c(.2, .2, .2, .2), "cm")
        ) +
        ggplot2::coord_cartesian(xlim = c(-f_lim, f_lim),
                                 ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(
          breaks = seq(-f_lim, f_lim, (f_lim / 4))
        )
    }
    if(missing(ct) == FALSE && class(l_deg) == "data.frame") { # nolint
      # Subset Input data
      ld <- l_deg
      if (comp_filt == TRUE) {
        ld <- ld[ld[["Comparison"]] == comp_name, ]
      }
      ld <- ld[ld[["CellType"]] == ct, ]
      # Plot function
      v <- EnhancedVolcano::EnhancedVolcano(
        ld,
        lab = ld[[gene_name]],
        selectLab = filt_lab,
        title = ggplot2::element_blank(),
        subtitle = ggplot2::element_blank(),
        caption = ggplot2::element_blank(),
        x = diff_col,
        y = p_col,
        pCutoff = p_cut,
        FCcutoff = f_cut,
        cutoffLineType = "twodash",
        legendLabels = c("NS", "Fold Change",
                         "p-value", "FC+P"),
        legendLabSize = 12,
        labFace = "bold",
        col = ggsci::pal_npg("nrc")(10)[c(4, 3, 5, 8)],
        colAlpha = 0.7,
        legendIconSize = 4,
        pointSize = 2,
        border = "full",
        borderWidth = 1.5,
        legendPosition = "'right'",
        labSize = 4,
        drawConnectors = TRUE,
        typeConnectors = "open",
        min.segment.length = ggplot2::unit(
          1,
          "mm"
        )
      ) +
        sc_theme1() + # nolint
        ggplot2::labs(color = "Key") +
        ggplot2::labs(
          color = "Key",
          x = x_title,
          y = paste("-log10 P-value (FDR)")
        ) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(.2, .2, .2, .2), "cm")) +
        ggplot2::coord_cartesian(xlim = c(-f_lim, f_lim),
                                 ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(breaks = seq(-f_lim, f_lim, (f_lim / 4)))
    }
    if(missing(ct) == FALSE && class(l_deg) == "list") { # nolint
      # Subset Input data
      ld <- l_deg[[1]]
      if (comp_filt == TRUE) {
        ld <- ld[ld[["Comparison"]] == comp_name, ]
      }
      ld <- ld[ld[["CellType"]] == ct, ]
      # Plot function
      v <- EnhancedVolcano::EnhancedVolcano(
        ld,
        lab = ld[[gene_name]],
        selectLab = filt_lab,
        title = ggplot2::element_blank(),
        subtitle = ggplot2::element_blank(),
        caption = ggplot2::element_blank(),
        x = diff_col,
        y = p_col,
        pCutoff = p_cut,
        FCcutoff = f_cut,
        cutoffLineType = "twodash",
        legendLabels = c("NS", "Fold Change",
                         "p-value", "FC+P"),
        legendLabSize = 12,
        labFace = "bold",
        col = ggsci::pal_npg("nrc")(10)[c(4, 3, 5, 8)],
        colAlpha = 0.7,
        legendIconSize = 4,
        pointSize = 2,
        border = "full",
        borderWidth = 1.5,
        legendPosition = "'right'",
        labSize = 4,
        drawConnectors = TRUE,
        typeConnectors = "open",
        max.overlaps = 18,
        min.segment.length = ggplot2::unit(
          1,
          "mm"
        )
      ) +
        sc_theme1() + # nolint
        ggplot2::labs(color = "Key") +
        ggplot2::labs(
          color = "Key",
          x = x_title,
          y = paste("-log10 P-value (FDR)")
        ) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(.2, .2, .2, .2), "cm")) +
        ggplot2::coord_cartesian(xlim = c(-f_lim, f_lim),
                                 ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(breaks = seq(-f_lim, f_lim, (f_lim / 4)))
    }
  }
  return(v)
}
