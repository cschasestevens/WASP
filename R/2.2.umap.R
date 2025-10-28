#' Standard GEX/ATAC UMAP
#'
#' Generates a single UMAP with a specific
#' metadata overlay for a Seurat object
#' containing gene expression and/or chromatin accessibility data.
#'
#' @param so An input Seurat object.
#' @param md_var Name of the clustering column or overlay to plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @param dims1 Should the data be plotted in 2 or 3 dimensions?
#' @param col1 Color scheme to use for plot overlay.
#' @param pos_leg Legend position (specify x and y coordinates or type "none")
#' (either "2D" or "3D").
#' @param show_lab Should labels be shown on the UMAP? Set to FALSE if plotting
#' intensity data.
#' @return A UMAP with points grouped by a specific metadata column.
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @importFrom plotly add_markers plot_ly layout
#' @importFrom htmlwidgets saveWidget
#' @import ggrepel
#' @examples
#'
#' # p_umap <- sc_umap(
#' #   so = d,
#' #   md_var = "CellType",
#' #   slot1 = "umap"
#' # )
#'
#' @export
sc_umap <- function(
  so,
  md_var,
  slot1 = "umap",
  dims1 = "2D",
  col1 = col_univ(), # nolint
  pos_leg = "none",
  show_lab = TRUE
) {
  # Format input data
  d <- so
  if (
    md_var %in%
      c(rownames(d), names(d@meta.data)) == FALSE
  ) {
    print(
      paste(
        md_var,
        " is not present as a variable in the dataset!",
        " Initializing an empty column for the selected variable...",
        sep = ""
      )
    )
    d <- Seurat::AddMetaData(
      object = d,
      metadata = rep(0, ncol(d)),
      col.name = md_var
    )
  }
  # Extract data from Seurat object
  if(is.null(slot1) == FALSE) { # nolint
    if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3) { # nolint
      d2 <- data.frame(
        Seurat::FetchData(d, vars = c(md_var)),
        `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
        `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
        `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[, 3]
      )
      if(class(d2[[gsub("\\/|\\-|\\ ", ".", md_var)]]) == "character") { # nolint
        d2[[gsub("\\/|\\-|\\ ", ".", md_var)]] <- factor(
          d2[[gsub("\\/|\\-|\\ ", ".", md_var)]],
          levels = sort(unique(d2[[gsub("\\/|\\-|\\ ", ".", md_var)]]))
        )
      }
    }
    if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2) { # nolint
      d2 <- data.frame(
        Seurat::FetchData(d, vars = c(md_var)),
        `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
        `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
      )
      if(class(d2[[gsub("\\/|\\-|\\ ", ".", md_var)]]) == "character") { # nolint
        d2[[gsub("\\/|\\-|\\ ", ".", md_var)]] <- factor(
          d2[[gsub("\\/|\\-|\\ ", ".", md_var)]],
          levels = sort(unique(d2[[gsub("\\/|\\-|\\ ", ".", md_var)]]))
        )
      }
    }
  }
  # Set color scheme
  if (class(d2[[gsub("\\/|\\-|\\ ", ".", md_var)]]) == "factor") {
    scm <- ggplot2::scale_color_manual(
      paste(""),
      values = col1
    )
  }
  if (class(d2[[gsub("\\/|\\-|\\ ", ".", md_var)]]) == "numeric") {
    scm <- ggplot2::scale_color_gradientn(
      name = "Expression/Accessibility",
      colors = col1
    )
  }
  if(dims1 == "3D") { # nolint
    tryCatch(
      {
        # Generate plot
        d2_plot <- plotly::plot_ly(
          d2,
          x = ~`UMAP.1`, # nolint
          y = ~`UMAP.2`, # nolint
          z = ~`UMAP.3`, # nolint
          color = ~.data[[gsub("\\/|\\-|\\ ", ".", md_var)]], # nolint
          colors = col1 # nolint
        ) %>% # nolint
          plotly::add_markers(marker = list(size = 3)) %>%
          plotly::layout(
            autosize = FALSE,
            width = 800,
            height = 600,
            margin = list(
              l = 50,
              r = 50,
              b = 25,
              t = 25,
              pad = 1
            )
          )
        htmlwidgets::saveWidget(
          d2_plot,
          file = paste(
            "analysis/p_umap_3d_", md_var, ".html", sep = ""
          )
        )
        d2_plot <- print(
          paste(
            "3D UMAP plot created as: ",
            "analysis/p_umap_3d_", md_var, ".html", sep = ""
          )
        )
      },
      error = function(e) {
        print("Error: data should contain the chosen overlay variable and three embedding variables named UMAP.1, UMAP.2, and UMAP.3!") # nolint
      }
    )
  }
  if(dims1 == "2D") { # nolint
    tryCatch(
      {
        if (show_lab == TRUE) {
          # Generate plot
          d2_plot <- ggplot2::ggplot(
            d2,
            ggplot2::aes(
              x=`UMAP.1`, # nolint
              y=`UMAP.2`, # nolint
              color = .data[[gsub("\\/|\\-|\\ ", ".", md_var)]], # nolint
              label = .data[[gsub("\\/|\\-|\\ ", ".", md_var)]] # nolint
            )
          ) +
            ggplot2::geom_point(
              shape = 16,
              size = 1,
              alpha = 0.6
            ) +
            ggrepel::geom_text_repel(
              data = setNames(
                aggregate(
                  d2[, c("UMAP.1", "UMAP.2")],
                  list(d2[[gsub("\\/|\\-|\\ ", ".", md_var)]]),
                  FUN = median
                ),
                c(gsub("\\/|\\-|\\ ", ".", md_var), names(
                  d2[, c("UMAP.1", "UMAP.2")]
                )
                )
              ),
              size = 5.5,
              bg.color = "grey0",
              color = "grey55",
              bg.r = 0.075
            ) +
            scm +
            sc_theme1() + # nolint
            ggplot2::ggtitle(md_var) +
            ggplot2::theme(
              panel.grid.major.y = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              plot.margin = ggplot2::unit(
                c(0.1, 0.1, 0.1, 0.1), "cm"
              ),
              legend.position = pos_leg
            ) +
            ggplot2::labs(title = md_var)
        }
        if (show_lab == FALSE) {
          # Generate plot
          d2_plot <- ggplot2::ggplot(
            d2,
            ggplot2::aes(
              x=`UMAP.1`, # nolint
              y=`UMAP.2`, # nolint
              color = .data[[gsub("\\/|\\-|\\ ", ".", md_var)]], # nolint
              label = .data[[gsub("\\/|\\-|\\ ", ".", md_var)]] # nolint
            )
          ) +
            ggplot2::geom_point(
              shape = 16,
              size = 1,
              alpha = 0.6
            ) +
            scm +
            sc_theme1() + # nolint
            ggplot2::ggtitle(md_var) +
            ggplot2::theme(
              panel.grid.major.y = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              plot.margin = ggplot2::unit(
                c(0.1, 0.1, 0.1, 0.1), "cm"
              ),
              legend.position = pos_leg
            ) +
            ggplot2::labs(title = md_var)
        }
      },
      error = function(e) {
        print("Error: data should contain the chosen overlay variable and two embedding variables named UMAP.1 and UMAP.2!") # nolint
      }
    )
  }
  return(d2_plot)
}

#' Violin Plot
#'
#' Generates violin plot(s) with a specific
#' metadata overlay for a Seurat object
#' containing gene expression/chromatin accessibility data.
#'
#' @param so An input Seurat object.
#' @param ct_col Cell Type column or metadata variable to plot.
#' @param gene_col A vector of one or more genes to plot.
#' @param pos_leg Legend position (specify x and y coordinates or type "none").
#' @return A violin plot grouped by a specific metadata column.
#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#' @import ggpubr
#' @examples
#'
#' # p_umap <- sc_violin(
#' #   so = d,
#' #   gene_col = c("DAPI", "EpCAM")
#' # )
#'
#' @export
sc_violin <- function(
  so,
  ct_col = "CellType",
  gene_col,
  pos_leg = "none"
) {
  # Input data
  d <- so
  ab1 <- gene_col[gene_col %in% rownames(d)]
  if (
    length(gene_col[gene_col %in% rownames(d) == FALSE] >= 1)
  ) {
    print("The following genes are not present in the Seurat object:")
    print(gene_col[gene_col %in% rownames(d) == FALSE])
    print("Initializing an empty column for the selected gene...")
    ab1 <- gene_col[gene_col %in% rownames(d) == FALSE]
    d <- Seurat::AddMetaData(
      object = d,
      metadata = rep(0, ncol(d)),
      ab1
    )
  }
  d <- Seurat::FetchData(
    d,
    vars = c(ct_col, ab1)
  )
  if (class(d[[ct_col]]) == "character") {
    d[[ct_col]] <- factor(
      d[[ct_col]],
      levels = sort(unique(d[[ct_col]]))
    )
  }
  # Plot data
  if (length(ab1) == 1) {
    plot_v <- ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = .data[[ct_col]], # nolint
        y = .data[[ab1]],
        fill = .data[[ct_col]]
      )
    ) +
      ggplot2::scale_fill_manual(
        name = ct_col,
        values = col_univ()[1:length( # nolint
          levels(d[[ct_col]])
        )
        ]
      ) +
      # Add violin plot and dotplot
      ggplot2::geom_violin(
        trim = TRUE
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          alpha = 0.2
        ),
        shape = 16,
        size = 0.2,
        position = ggplot2::position_jitter(
          width = 0.4
        ),
        show.legend = FALSE
      ) +
      # Add Theme
      sc_theme1() + # nolint
      ggplot2::labs(
        y = "GEX/ATAC"
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(
            0.1,
            0.1,
            0.1,
            0.1
          ),
          "cm"
        ),
        legend.position = pos_leg
      )
  }
  if (length(ab1) > 1) {
    plot_v <- setNames(lapply(
      seq.int(1, length(ab1), 1),
      function(i) {
        p1 <- ggplot2::ggplot(
          d,
          ggplot2::aes(
            x = .data[[ct_col]], # nolint
            y = .data[[ab1[[i]]]],
            fill = .data[[ct_col]]
          )
        ) +
          ggplot2::scale_fill_manual(
            name = ct_col,
            values = col_univ()[1:length( # nolint
              levels(d[[ct_col]])
            )
            ]
          ) +
          # Add violin plot and dotplot
          ggplot2::geom_violin(
            trim = TRUE
          ) +
          ggplot2::geom_jitter(
            ggplot2::aes(
              alpha = 0.2
            ),
            shape = 16,
            size = 0.2,
            position = ggplot2::position_jitter(
              width = 0.4
            ),
            show.legend = FALSE
          ) +
          # Add Theme
          sc_theme1() + # nolint
          ggplot2::labs(
            y = "GEX/ATAC"
          ) +
          ggplot2::theme(
            plot.margin = ggplot2::unit(
              c(
                0.1,
                0.1,
                0.1,
                0.1
              ),
              "cm"
            ),
            legend.position = pos_leg
          )
        return(p1) # nolint
      }
    ), ab1)
    plot_v <- ggpubr::ggarrange(
      plotlist = plot_v,
      ncol = ifelse(
        length(ab1) > 4,
        4,
        length(ab1)
      ),
      nrow = ifelse(
        length(ab1) > 4,
        ceiling(length(ab1) / 4),
        1
      ),
      labels = names(plot_v),
      common.legend = TRUE
    )
  }
  return(plot_v) # nolint
}

#' Cell Proportion Bar Plot
#'
#' Generates bar plot summarizing the number of cells per each
#' type for all treatment groups.
#'
#' @param dat An input data frame containing counts.
#' @param ct_col Cell Type column or metadata variable to plot.
#' @param xvar An additional metadata variable for stratifying
#' the input data.
#' @param yvar Cell count variable name.
#' @return A bar plot grouped by a specific metadata column.
#' @import ggplot2
#' @import ggpubr
#' @examples
#'
#' # sc_barplot(dat = d)
#'
#' @export
sc_barplot <- function(
  dat,
  ct_col = "CellType",
  xvar = "Group",
  yvar = "Proportion"
) {
  pv <- ggplot2::ggplot(
    data = dat,
    ggplot2::aes(
      x = as.factor(.data[[xvar]]), # nolint
      y = .data[[yvar]],
      fill = .data[[ct_col]]
    )
  ) +
    ggplot2::scale_fill_manual(
      name = ct_col,
      values = col_univ()
    ) +
    # Add barplot
    ggplot2::geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.6,
      color = "grey25"
    ) +
    # Add Theme
    sc_theme1() + # nolint
    ggplot2::labs(
      y = yvar,
      x = xvar
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::unit(
        c(
          0.5,
          0.5,
          0.5,
          0.5
        ),
        "cm"
      ),
      axis.text.y = ggplot2::element_text(
        face = "bold",
        size = 12
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        size = 12
      )
    )
  return(pv) # nolint
}














#' scRNA-Seq UMAP Plot
#'
#' Generates a panel of UMAPs given a
#' Seurat object containing PCA and UMAP results.
#'
#' @param so An object of class Seurat.
#' @param md_list A vector of character strings indicating
#' metadata columns for overlaying on a loadings plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A series of UMAPs with specified metadata overlays.
#' @examples
#'
#' # p_umap <- sc_umap_panel(d_integrated,c("col1","col2","col3"),"wnn.umap")
#'
#' @export
sc_umap_panel <- function(
  so,
  md_list,
  slot1,
  leg_x = 0.9,
  leg_y = 0.25
) {
  d <- so
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3) { #nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
      `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[, 3]
    )
    d2_list <- setNames(
      lapply(
        c(md_list),
        function(x) {
          d2[, c(
            x,
            "UMAP.1",
            "UMAP.2",
            "UMAP.3"
          )
          ]
        }
      ),
      c(md_list)
    )

  }
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2) { #nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
    )
    d2_list <- setNames(
      lapply(
        c(md_list),
        function(x) {
          d2[, c(
            x,
            "UMAP.1",
            "UMAP.2"
          )
          ]
        }
      ),
      c(md_list)
    )
  }
  # Generate plots
  d2_plot <- lapply(
    c(md_list),
    function(x) {
      p <-  ggplot2::ggplot(
        d2_list[[x]],
        ggplot2::aes(
          x=`UMAP.1`, # nolint
          y=`UMAP.2`, # nolint
          color = .data[[x]], # nolint
          label = .data[[x]] # nolint
        )
      ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col_univ() # nolint
        ) +
        # Add points
        ggplot2::geom_point(
          shape = 16,
          size = 1,
          alpha = 0.6
        ) +
        sc_theme1() + # nolint
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(
            c(0.1, 0.1, 0.1, 0.1),
            "cm"
          ),
          legend.position = c(
            leg_x,
            leg_y
          )
        )
      return(p) # nolint
    }
  )
  # Combine output
  if(length(d2_plot) > 1) { # nolint
    d2_out <- ggpubr::ggarrange(
      plotlist = d2_plot,
      labels = c(
        names(
          d2_plot
        )
      ),
      ncol = ifelse(
        length(
          d2_plot
        ) <= 3,
        length(
          d2_plot
        ),
        3
      ),
      nrow = ifelse(
        length(
          d2_plot
        ) > 3,
        ceiling(
          length(
            d2_plot
          ) /
            3
        ),
        1
      )
    )
  }
  if(length(d2_plot) == 1) { # nolint
    d2_out <- d2_plot[[1]]
  }
  return(d2_out)
}

#' Visualize Individual Gene Expression
#'
#' Generates a series of plots to visualize
#' individual gene expression per cluster.
#'
#' @param so An object of class Seurat.
#' @param asy1 Assay to use; either "sct" (default) or "chromvar".
#' @param md_var A character string indicating
#' the clustering column for overlaying on a UMAP plot.
#' @param g_name A character string indicating the gene name.
#' @param col_scheme The color scheme to be used for
#' distinguishing between groups, provided as a vector.
#' @param col_names A vector of the same length as the
#' provided color scheme for assigning colors to each group.
#' @param leg_x A numeric value indicating the placement
#' of the figure legend on the x-axis.
#' @param leg_y A numeric value indicating the placement
#' of the figure legend on the y-axis.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @param col1 Color scheme to use for expression gradient.
#' @param plot_comb Should all panels be displayed? (TRUE/FALSE)
#' @param out1 Panel to return if plot_comb is FALSE
#' (either "umap_gex" or "vio_gex").
#' @return A series of plots stored as a ggplot2 object
#' for visualizing cluster gene expression.
#' @examples
#'
#' # p_umap <- sc_umap_panel_gene(
#' #  so = d,
#' #  md_var = "CellType",
#' #  g_name = "CFTR",
#' #  col_scheme = col_univ(),
#' #  col_names = c("group1","group2"),
#' #  slot1 = "umap"
#' # )
#'
#' @export
sc_umap_panel_gene <- function(
  so,
  asy1 = "sct",
  md_var,
  g_name,
  col_scheme,
  col_names,
  leg_x = 0.9,
  leg_y = 0.2,
  slot1 = "umap",
  col1 = col_grad(), # nolint
  plot_comb = FALSE,
  out1 = "umap_gex"
) {
  # Format input data
  cols <- setNames(col_scheme,
                   col_names)
  d <- so
  SeuratObject::DefaultAssay(d) <- asy1
  if(is.null(slot1)) { # nolint
    d2 <- Seurat::FetchData(
      d,
      vars = c(
        gsub(
          "-",
          "-",
          g_name
        ),
        md_var
      )
    )
  }
  if(is.null(slot1) == FALSE) { # nolint
    d2 <- data.frame(
      Seurat::FetchData(
        d,
        vars = c(
          gsub(
            "-",
            "-",
            g_name
          ),
          md_var
        )
      ),
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
    )
  }
  if(asy1 == "chromvar" | asy1 == "ufy.peaks") { # nolint
    # Add motif names
    colnames(d2) <- c(
      name(TFBSTools::getMatrixByID(JASPAR2020, ID = g_name)), # nolint
      md_var,
      "UMAP.1", # nolint
      "UMAP.2" # nolint
    )
    # Generate plots
    d2_plot <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x=`UMAP.1`, # nolint
        y=`UMAP.2`, # nolint
        color = .data[[ # nolint
          name(TFBSTools::getMatrixByID(JASPAR2020, ID = g_name)) # nolint
        ]]
      )
    ) +
      ggplot2::scale_color_gradientn(
        name = "Accessibility",
        colors = col1
      ) +
      # Add 3D points, axes, and axis-labels
      ggplot2::geom_point(
        shape = 16,
        size = 1,
        alpha = 0.6
      ) +
      ggplot2::ggtitle(name(TFBSTools::getMatrixByID(JASPAR2020, ID = g_name))) + # nolint
      # Add general multivariate plot theme and adjust axis text
      sc_theme1() + # nolint
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(
          c(
            0.1, 0.1, 0.1, 0.1
          ),
          "cm"
        ),
        legend.position = c(
          leg_x,
          leg_y
        )
      )
    ## Metadata overlay
    p_md <-  ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x=`UMAP.1`, # nolint
        y=`UMAP.2`, # nolint
        color = .data[[md_var]], # nolint
        label = .data[[md_var]] # nolint
      )
    ) +
      ggplot2::scale_color_manual(
        paste(""),
        values = cols
      ) +
      ggplot2::geom_point(
        shape = 16,
        size = 1,
        alpha = 0.6
      ) +
      ggrepel::geom_text_repel(data = setNames(
        aggregate(
          d2[, c(
            "UMAP.1",
            "UMAP.2"
          )],
          list(
            d2[[md_var]]
          ),
          FUN = median
        ),
        c(
          md_var,
          names(
            d2[, c(
              "UMAP.1",
              "UMAP.2"
            )]
          )
        )
      ),
      size = 4,
      bg.color = "white") +
      ggplot2::ggtitle(md_var) +
      # Add general multivariate plot theme and adjust axis text
      sc_theme1() + # nolint
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(
          c(
            0.1, 0.1, 0.1, 0.1
          ),
          "cm"
        ),
        legend.position = "none"
      )
    ## Violin Plot
    plot_v <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x = .data[[md_var]], # nolint
        y = .data[[
          name(TFBSTools::getMatrixByID(JASPAR2020, ID = g_name)) # nolint
        ]],
        fill = .data[[md_var]]
      )
    ) +
      ggplot2::scale_fill_manual(
        name = md_var,
        values = col_scheme
      ) +
      # Add violin plot and dotplot
      ggplot2::geom_violin(
        trim = TRUE
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          alpha = 0.2
        ),
        shape = 16,
        size = 0.2,
        position = ggplot2::position_jitter(
          width = 0.4
        ),
        show.legend = FALSE
      ) +
      # Add Theme
      sc_theme1() + # nolint
      ggplot2::labs(
        y = "Accessibility"
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(
            0.1,
            0.1,
            0.1,
            0.1
          ),
          "cm"
        ),
        legend.position = "none"
      )
  }
  if(asy1 == "sct" | asy1 == "PC") { # nolint
    # Generate plots
    d2_plot <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x=`UMAP.1`, # nolint
        y=`UMAP.2`, # nolint
        color = .data[[gsub( # nolint
          "-",
          ".",
          g_name
        )]]
      )
    ) +
      ggplot2::scale_color_gradientn(
        name = "Relative Exp.",
        colors = col1
      ) +
      # Add 3D points, axes, and axis-labels
      ggplot2::geom_point(
        shape = 16,
        size = 1,
        alpha = 0.6
      ) +
      ggplot2::ggtitle(g_name) +
      # Add general multivariate plot theme and adjust axis text
      sc_theme1() + # nolint
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(
          c(
            0.1, 0.1, 0.1, 0.1
          ),
          "cm"
        ),
        legend.position = c(
          leg_x,
          leg_y
        )
      )
    ## Metadata overlay
    p_md <-  ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x=`UMAP.1`, # nolint
        y=`UMAP.2`, # nolint
        color = .data[[md_var]], # nolint
        label = .data[[md_var]] # nolint
      )
    ) +
      ggplot2::scale_color_manual(
        paste(""),
        values = cols
      ) +
      ggplot2::geom_point(
        shape = 16,
        size = 1,
        alpha = 0.6
      ) +
      ggrepel::geom_text_repel(data = setNames(
        aggregate(
          d2[, c(
            "UMAP.1",
            "UMAP.2"
          )],
          list(
            d2[[md_var]]
          ),
          FUN = median
        ),
        c(
          md_var,
          names(
            d2[, c(
              "UMAP.1",
              "UMAP.2"
            )]
          )
        )
      ),
      size = 4,
      bg.color = "white") +
      ggplot2::ggtitle(md_var) +
      # Add general multivariate plot theme and adjust axis text
      sc_theme1() + # nolint
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(
          c(
            0.1, 0.1, 0.1, 0.1
          ),
          "cm"
        ),
        legend.position = "none"
      )
    ## Violin Plot
    plot_v <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x = .data[[md_var]], # nolint
        y = .data[[gsub(
          "-",
          ".",
          g_name
        )]],
        fill = .data[[md_var]]
      )
    ) +
      ggplot2::scale_fill_manual(
        name = md_var,
        values = col_scheme
      ) +
      # Add violin plot and dotplot
      ggplot2::geom_violin(
        trim = TRUE
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          alpha = 0.2
        ),
        shape = 16,
        size = 0.2,
        position = ggplot2::position_jitter(
          width = 0.4
        ),
        show.legend = FALSE
      ) +
      # Add Theme
      sc_theme1() + # nolint
      ggplot2::labs(
        y = "Relative Expression"
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(
            0.1,
            0.1,
            0.1,
            0.1
          ),
          "cm"
        ),
        legend.position = "none"
      )
  }
  if(asy1 == "SCT") { # nolint
    # Generate plots
    d2_plot <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x=`UMAP.1`, # nolint
        y=`UMAP.2`, # nolint
        color = .data[[gsub( # nolint
          "-",
          ".",
          g_name
        )]]
      )
    ) +
      ggplot2::scale_color_gradientn(
        name = "Relative Exp.",
        colors = col1
      ) +
      # Add 3D points, axes, and axis-labels
      ggplot2::geom_point(
        shape = 16,
        size = 1,
        alpha = 0.6
      ) +
      ggplot2::ggtitle(g_name) +
      # Add general multivariate plot theme and adjust axis text
      sc_theme1() + # nolint
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(
          c(
            0.1, 0.1, 0.1, 0.1
          ),
          "cm"
        ),
        legend.position = c(
          leg_x,
          leg_y
        )
      )
    ## Metadata overlay
    p_md <-  ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x=`UMAP.1`, # nolint
        y=`UMAP.2`, # nolint
        color = .data[[md_var]], # nolint
        label = .data[[md_var]] # nolint
      )
    ) +
      ggplot2::scale_color_manual(
        paste(""),
        values = cols
      ) +
      ggplot2::geom_point(
        shape = 16,
        size = 1,
        alpha = 0.6
      ) +
      ggrepel::geom_text_repel(data = setNames(
        aggregate(
          d2[, c(
            "UMAP.1",
            "UMAP.2"
          )],
          list(
            d2[[md_var]]
          ),
          FUN = median
        ),
        c(
          md_var,
          names(
            d2[, c(
              "UMAP.1",
              "UMAP.2"
            )]
          )
        )
      ),
      size = 4,
      bg.color = "white") +
      ggplot2::ggtitle(md_var) +
      # Add general multivariate plot theme and adjust axis text
      sc_theme1() + # nolint
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(
          c(
            0.1, 0.1, 0.1, 0.1
          ),
          "cm"
        ),
        legend.position = "none"
      )
    ## Violin Plot
    plot_v <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x = .data[[md_var]], # nolint
        y = .data[[gsub(
          "-",
          ".",
          g_name
        )]],
        fill = .data[[md_var]]
      )
    ) +
      ggplot2::scale_fill_manual(
        name = md_var,
        values = col_scheme
      ) +
      # Add violin plot and dotplot
      ggplot2::geom_violin(
        trim = TRUE
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          alpha = 0.2
        ),
        shape = 16,
        size = 0.2,
        position = ggplot2::position_jitter(
          width = 0.4
        ),
        show.legend = FALSE
      ) +
      # Add Theme
      sc_theme1() + # nolint
      ggplot2::labs(
        y = "Relative Expression"
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(
            0.1,
            0.1,
            0.1,
            0.1
          ),
          "cm"
        ),
        legend.position = "none"
      )
  }
  # Combine output or choose slot
  if(plot_comb == TRUE) { # nolint
    d2_out <- ggpubr::ggarrange(
      d2_plot,
      p_md,
      plot_v,
      ncol = 3,
      nrow = 1,
      common.legend = FALSE
    )
  }
  if(plot_comb == FALSE && out1 == "umap_gex") { # nolint
    d2_out <- d2_plot
  }
  if(plot_comb == FALSE && out1 == "vio_gex") { # nolint
    d2_out <- plot_v
  }
  return(d2_out)
}

#' Visualize Individual Gene Expression
#'
#' Generates a series of plots to visualize
#' individual gene expression per cluster.
#'
#' @param so An object of class Seurat.
#' @param md_var A character string indicating
#' the clustering column for overlaying on a UMAP plot.
#' @param g_name A character string indicating the gene name.
#' @param col_scheme The color scheme to be used for
#' distinguishing between groups, provided as a vector.
#' @param col_names A vector of the same length as the
#' provided color scheme for assigning colors to each group.
#' @param leg_x A numeric value indicating the placement
#' of the figure legend on the x-axis.
#' @param leg_y A numeric value indicating the placement
#' of the figure legend on the y-axis.
#' @param asy1 Default assay for plotting violin plot ("RNA" by default).
#' for visualizing cluster gene expression.
#' @examples
#'
#' # p_umap <- sc_umap_panel_gene(
#' #  d_integrated,
#' #  c("col1","col2","col3"),
#' #  "CFTR",
#' #  col_univ,
#' #  c("group1","group2"),
#' #  0.95,
#' #  0.95
#' # )
#'
#' @export
sc_vio_panel_gene <- function(
  so,
  md_vars,
  g_name,
  col_scheme = col_univ(), # nolint
  col_names = NULL,
  leg_x = 0.95,
  leg_y = 0.95,
  asy1 = "RNA"
) {
  # Format input data
  if(!is.null(col_names)) { # nolint
    cols <- setNames(
      col_scheme,
      col_names
    )
  }
  d1 <- so
  Seurat::DefaultAssay(d1) <- asy1
  d2 <- data.frame(
    Seurat::FetchData(
      d1,
      vars = c(
        gsub(
          "-",
          "-",
          g_name
        ),
        md_vars
      )
    )
  )
  if(length(md_vars) == 1) { # nolint
    ## Violin Plot
    plot_v <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x = .data[[md_vars]], # nolint
        y = .data[[gsub(
          "-",
          ".",
          g_name
        )]],
        fill = .data[[md_vars]]
      )
    ) +
      ggplot2::scale_fill_manual(
        name = md_vars,
        values = cols
      ) +
      # Add violin plot and dotplot
      ggplot2::geom_violin(
        trim = TRUE
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          alpha = 0.2
        ),
        shape = 16,
        size = 0.2,
        position = ggplot2::position_jitter(
          width = 0.4
        ),
        show.legend = FALSE
      ) +
      # Add Theme
      sc_theme1() + # nolint
      ggplot2::labs(
        y = "Relative Expression",
        x = "Cell Type"
      ) +
      ggplot2::ggtitle(g_name) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(
            0.1,
            0.1,
            0.1,
            0.1
          ),
          "cm"
        ),
        legend.position = "none"
      )
  }
  if(length(md_vars) == 2) { # nolint
    ## Violin Plot
    plot_v <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x = .data[[md_vars[[1]]]], # nolint
        y = .data[[gsub(
          "-",
          ".",
          g_name
        )]],
        fill = .data[[md_vars[[2]]]]
      )
    ) +
      ggplot2::scale_fill_manual(
        name = md_vars[[2]],
        values = col_univ() # nolint
      ) +
      # Add violin plot and dotplot
      ggplot2::geom_violin(
        trim = TRUE
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          alpha = 0.2
        ),
        shape = 16,
        size = 0.8,
        position = ggplot2::position_jitter(
          width = 0.4
        ),
        show.legend = FALSE
      ) +
      # Add Theme
      sc_theme1() + # nolint
      ggplot2::labs(
        y = "Relative Expression",
        x = "Cell Type"
      ) +
      ggplot2::ggtitle(g_name) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(
            0.1,
            0.1,
            0.1,
            0.1
          ),
          "cm"
        )
      )
  }
  return(plot_v)
}

#' Visualize Gene List Expression
#'
#' Generates a series of plots to visualize
#' the expression of a list of genes per cluster.
#'
#' @param list_g A vector containing a list of genes
#' to be plotted for a given Seurat object.
#' @param so An object of class Seurat.
#' @param md_var A character string indicating the clustering
#' column for overlaying on a UMAP plot.
#' @param col_scheme The color scheme to be used for distinguishing
#' between groups, provided as a vector.
#' @param col_names A vector of the same length as the provided
#' color scheme for assigning colors to each group.
#' @param leg_x A numeric value indicating the placement
#' of the figure legend on the x-axis.
#' @param leg_y A numeric value indicating the placement
#' of the figure legend on the y-axis.
#' @param parl Logical indicating whether processing should
#' be run in parallel (Linux and WSL2 only).
#' Set to FALSE if running sequentially.
#' @param core_perc Percentage of available cores to use if
#' running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A list of plots saved as ggplot2 objects for
#' visualizing cluster gene expression.
#' @examples
#'
#' # sc_umap_panel_gene_list(
#' #  list_genes,
#' #  d_seurat,
#' #  "seurat_clusters",
#' #  col_vec,
#' #  col_vec_names,
#' #  0.95,
#' #  0.95,
#' #  TRUE,
#' #  0.5,
#' #  "umap"
#' # )
#'
#' @export
sc_umap_panel_gene_list <- function(
  list_g,
  so,
  md_var,
  col_scheme,
  col_names,
  leg_x,
  leg_y,
  parl,
  core_perc,
  slot1
) {
  lg <- list_g
  d <- so
  lg <- unique(lg[lg %in% SeuratObject::Features(d)])
  lg_abs <- subset(lg, !(lg %in% SeuratObject::Features(d)))
  # Create plots
  if(Sys.info()[["sysname"]] != "Windows" && # nolint
      parl == TRUE
  ) {
    parallel::mclapply(
      mc.cores = ceiling(
        parallel::detectCores() *
          core_perc
      ),
      lg,
      function(x) {
        pg <- sc_umap_panel_gene(
          d,
          md_var,
          x,
          col_scheme,
          col_names,
          leg_x,
          leg_y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste(
            "analysis/gene/plot.umap.exp.",
            x,
            ".png",
            sep = ""
          ),
          pg,
          height = 12,
          width = 36,
          dpi = 700
        )
      }
    )

  }
  if(Sys.info()[["sysname"]] == "Windows") { # nolint
    lapply(
      lg,
      function(x) {
        pg <- sc_umap_panel_gene(
          d,
          md_var,
          x,
          col_scheme,
          col_names,
          leg_x,
          leg_y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste(
            "analysis/gene/plot.umap.exp.",
            x,
            ".png",
            sep = ""
          ),
          pg,
          height = 12,
          width = 36,
          dpi = 700
        )
      }
    )
  }
  if(Sys.info()[["sysname"]] != "Windows" && parl == FALSE) { # nolint
    lapply(
      lg,
      function(x) {
        pg <- sc_umap_panel_gene(
          d,
          md_var,
          x,
          col_scheme,
          col_names,
          leg_x,
          leg_y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste(
            "analysis/gene/plot.umap.exp.",
            x,
            ".png",
            sep = ""
          ),
          pg,
          height = 12,
          width = 36,
          dpi = 700
        )
      }
    )

  }
  if(length(lg_abs) > 0) { # nolint
    print(
      paste(
        lg_abs,
        "was not found; plots for this gene will be 
        excluded from the final list...",
        sep = " "
      )
    )
  }
}

#' Standard UMAP Plot
#'
#' Generates a single UMAP plot with a specific
#' metadata overlay for a Seurat object.
#'
#' @param so An object of class Seurat.
#' @param md_var A character string indicating the clustering
#' column for overlaying on a UMAP plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @param dims1 Should the data be plotted in 2 or 3 dimensions?
#' @param col1 Color scheme to use for plot overlay.
#' @param pos_leg Legend position (specify x and y coordinates or type "none")
#' (either "2D" or "3D").
#' @return A UMAP plot with points grouped by a specific metadata column.
#' @examples
#'
#' # p_umap <- sc_umap_standard(
#' #   so = d_integrated,
#' #   md_var = "CellType",
#' #   slot1 = "wnn.umap"
#' # )
#'
#' @export
sc_umap_standard <- function(
  so,
  md_var,
  slot1 = NULL,
  dims1 = "2D",
  col1 = col_univ(), # nolint
  pos_leg = "none"
) {
  # Format input data
  d <- so
  # If Seurat Object contains embeddings as metadata
  if(is.null(slot1)) { # nolint
    d2 <- d@meta.data
    if(class(d2[[md_var]]) == "character") { # nolint
      d2[[md_var]] <- factor(
        d2[[md_var]],
        levels = sort(unique(d2[[md_var]]))
      )
    }
  }
  # If Seurat Object contains dimension reductions (most often)
  if(is.null(slot1) == FALSE) { # nolint
    if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3) { # nolint
      d2 <- data.frame(
        d@meta.data,
        `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
        `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
        `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[, 3],
        md.var = d@meta.data[[md_var]]
      )
      if(class(d2[[md_var]]) == "character") { # nolint
        d2[[md_var]] <- factor(
          d2[[md_var]],
          levels = sort(unique(d2[[md_var]]))
        )
      }
    }
    if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2) { # nolint
      d2 <- data.frame(
        d@meta.data,
        `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
        `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
        md.var = d@meta.data[[md_var]]
      )
      if(class(d2[[md_var]]) == "character") { # nolint
        d2[[md_var]] <- factor(
          d2[[md_var]],
          levels = sort(unique(d2[[md_var]]))
        )
      }
    }
  }
  if(dims1 == "3D") { # nolint
    tryCatch(
      {
        # Generate plot
        d2_plot <- plotly::plot_ly(
          d2,
          x = ~`UMAP.1`, # nolint
          y = ~`UMAP.2`, # nolint
          z = ~`UMAP.3`, # nolint
          color = ~.data[[md_var]], # nolint
          colors = col_univ() # nolint
        ) %>% # nolint
          plotly::add_markers(marker = list(size = 3)) %>%
          plotly::layout(
            autosize = FALSE,
            width = 800,
            height = 600,
            margin = list(
              l = 50,
              r = 50,
              b = 25,
              t = 25,
              pad = 1
            )
          )
        htmlwidgets::saveWidget(
          d2_plot,
          file = paste(
            "analysis/p_umap_3d_", md_var, ".html", sep = ""
          )
        )
        d2_plot <- print(
          paste(
            "3D UMAP plot created as: ",
            "analysis/p_umap_3d_", md_var, ".html", sep = ""
          )
        )
      },
      error = function(e) {
        print("Error: data should contain the chosen overlay variable and three embedding variables named UMAP.1, UMAP.2, and UMAP.3!") # nolint
      }
    )
  }
  if(dims1 == "2D") { # nolint
    tryCatch(
      {
        # Generate plot
        d2_plot <- ggplot2::ggplot(
          d2,
          ggplot2::aes(
            x=`UMAP.1`, # nolint
            y=`UMAP.2`, # nolint
            color = .data[[md_var]], # nolint
            label = .data[[md_var]] # nolint
          )
        ) +
          ggplot2::geom_point(
            shape = 16,
            size = 1,
            alpha = 0.6
          ) +
          ggrepel::geom_text_repel(
            data = setNames(
              aggregate(
                d2[, c("UMAP.1", "UMAP.2")],
                list(d2[[md_var]]),
                FUN = median
              ),
              c(md_var, names(
                d2[, c("UMAP.1", "UMAP.2")]
              )
              )
            ),
            size = 5.5,
            bg.color = "grey0",
            color = "grey55",
            bg.r = 0.075
          ) +
          ggplot2::scale_color_manual(
            paste(""),
            values = col1
          ) +
          sc_theme1() + # nolint
          ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(
              c(0.1, 0.1, 0.1, 0.1), "cm"
            ),
            legend.position = pos_leg
          )
      },
      error = function(e) {
        print("Error: data should contain the chosen overlay variable and two embedding variables named UMAP.1 and UMAP.2!") # nolint
      }
    )
  }
  return(d2_plot)
}

#' scATAC-Seq Activity UMAP
#'
#' Generates a series of plots given a Seurat scATAC-Seq assay.
#'
#' @param so A Seurat object.
#' @param md_var A character string indicating the clustering
#' column for overlaying on a UMAP plot.
#' @param da_df Differential activity results data frame.
#' @param tf_name Transcription factor name, provided as a character string.
#' @param col_scheme Colors to use for overlaying metadata variables.
#' @param col_names Names for assigning colors from col_scheme.
#' @param leg_x Legend x-axis position.
#' @param leg_y Legend y-axis position.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A UMAP plot with points indicating the activity of a specific
#' transcription factor motif.
#' @examples
#'
#' # p_umap_act <- sc_umap_panel_act(
#' #   d,
#' #   "CellType",
#' #   diff.activity.output,
#' #   "GHRL1",
#' #   col_univ()[1:length(levels(d@meta.data[[md_var]]))],
#' #   c(levels(d@meta.data[[md_var]]))
#' #   0.95,
#' #   0.95,
#' #   "wnn.umap.cor"
#' # )
#'
#' @export
sc_umap_panel_act <- function(
  so,
  md_var,
  da_df,
  tf_name,
  col_scheme,
  col_names,
  leg_x,
  leg_y,
  slot1
) {
  d <- so
  # Format input data
  cols <- setNames(col_scheme, col_names)

  d2 <- data.frame(
    Seurat::FetchData(
      d,
      vars = c(
        da_df[da_df[["gene"]] == tf_name, "ID"][[1]],
        md_var
      )
    ),
    `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
    `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
  )

  # Generate plots
  d2_plot <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[[da_df[da_df[["gene"]] == tf_name, "ID"][[1]]]] # nolint
    )
  ) +
    ggplot2::scale_color_gradientn(
      name = "Motif Activity",
      colors = col_grad() # nolint
    ) +
    # Add 3D points, axes, and axis-labels
    ggplot2::geom_point(
      shape = 16,
      size = 1,
      alpha = 0.6
    ) +
    ggplot2::ggtitle(tf_name) +
    # Add general multivariate plot theme and adjust axis text
    sc_theme1() + # nolint
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1, 0.1, 0.1, 0.1
        ),
        "cm"
      ),
      legend.position = c(
        leg_x,
        leg_y
      )
    )
  ## Metadata overlay
  p_md <-  ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[[md_var]], # nolint
      label = .data[[md_var]] # nolint
    )
  ) +
    ggplot2::scale_color_manual(
      paste(""),
      values = cols
    ) +
    ggplot2::geom_point(
      shape = 16,
      size = 1,
      alpha = 0.6
    ) +

    ggrepel::geom_text_repel(data = setNames(
      aggregate(
        d2[, c(
          "UMAP.1",
          "UMAP.2"
        )],
        list(
          d2[[md_var]]
        ),
        FUN = median
      ),
      c(
        md_var,
        names(
          d2[, c(
            "UMAP.1",
            "UMAP.2"
          )]
        )
      )
    ),
    size = 4,
    bg.color = "white") +
    ggplot2::ggtitle(md_var) +
    # Add general multivariate plot theme and adjust axis text
    sc_theme1() + # nolint
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1, 0.1, 0.1, 0.1
        ),
        "cm"
      ),
      legend.position = "none"
    )
  ## Violin Plot
  plot_v <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x = .data[[md_var]], # nolint
      y = .data[[da_df[da_df[["gene"]] == tf_name, "ID"][[1]]]],
      fill = .data[[md_var]]
    )
  ) +
    ggplot2::scale_fill_manual(
      name = md_var,
      values = col_univ()[1:length( # nolint
        levels(d@meta.data[[md_var]])
      )
      ]
    ) +
    # Add violin plot and dotplot
    ggplot2::geom_violin(
      trim = TRUE
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(
        alpha = 0.2
      ),
      shape = 16,
      size = 0.2,
      position = ggplot2::position_jitter(
        width = 0.4
      ),
      show.legend = FALSE
    ) +
    # Add Theme
    sc_theme1() + # nolint
    ggplot2::labs(
      y = "Motif Activity"
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::unit(
        c(
          0.1,
          0.1,
          0.1,
          0.1
        ),
        "cm"
      ),
      legend.position = "none"
    )
  # Combine output
  d2_out <- ggpubr::ggarrange(
    d2_plot,
    p_md,
    plot_v,
    ncol = 3,
    nrow = 1,
    common.legend = FALSE
  )
  return(d2_out) # nolint
}