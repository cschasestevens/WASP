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
#' @param grp_col Group column if splitting violin plots by
#' an additional variable.
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
  grp_col = NULL,
  gene_col,
  pos_leg = "none"
) {
  # Input data
  d <- so
  ab1 <- gene_col[gene_col %in% rownames(d) | gene_col %in% names(d@meta.data)]
  if (
    length(gene_col[
      gene_col %in% rownames(d) == FALSE &
        gene_col %in% names(d@meta.data) == FALSE
    ] >= 1)
  ) {
    print("The following genes or columns are not present in the Seurat object:") # nolint
    print(gene_col[
      gene_col %in% rownames(d) == FALSE &
        gene_col %in% names(d@meta.data) == FALSE
    ])
    print("Initializing an empty column for the selected gene...")
    ab1 <- gene_col[
      gene_col %in% rownames(d) == FALSE &
        gene_col %in% names(d@meta.data) == FALSE
    ]
    d <- Seurat::AddMetaData(
      object = d,
      metadata = rep(0, ncol(d)),
      ab1
    )
  }
  if (is.null(grp_col)) {
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
  }
  if (!is.null(grp_col)) {
    d <- Seurat::FetchData(
      d,
      vars = c(ct_col, grp_col, ab1)
    )
    if (class(d[[ct_col]]) == "character") {
      d[[ct_col]] <- factor(
        d[[ct_col]],
        levels = sort(unique(d[[ct_col]]))
      )
    }
    if (class(d[[grp_col]]) == "character") {
      d[[grp_col]] <- factor(
        d[[grp_col]],
        levels = sort(unique(d[[grp_col]]))
      )
    }
  }
  # Plot data
  if (is.null(grp_col)) {
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
          y = "GEX/ATAC/Value"
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
              y = "GEX/ATAC/Value"
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
  }
  if (!is.null(grp_col)) {
    if (length(ab1) == 1) {
      plot_v <- ggplot2::ggplot(
        d,
        ggplot2::aes(
          x = .data[[ct_col]], # nolint
          y = .data[[ab1]],
          fill = .data[[grp_col]]
        )
      ) +
        ggplot2::scale_fill_manual(
          name = grp_col,
          values = col_univ()[1:length( # nolint
            levels(d[[grp_col]])
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
          y = "GEX/ATAC/Value"
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
              fill = .data[[grp_col]]
            )
          ) +
            ggplot2::scale_fill_manual(
              name = grp_col,
              values = col_univ()[1:length( # nolint
                levels(d[[grp_col]])
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
              y = "GEX/ATAC/Value"
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

#' UMAP panels
#'
#' Generates a panel of UMAPs from a Seurat object.
#'
#' @param so A Seurat object.
#' @param md_list A vector of metadata columns to plot.
#' @param slot1 Name of the dimension reduction slot to plot.
#' @param pos_leg Legend position provided as a vector (default is "none").
#' @param show_lab Show labels on UMAP?
#' @return A series of UMAPs with specified metadata overlays.
#' @examples
#'
#' # p_umap <- sc_umap_panel(d_integrated,c("col1","col2","col3"),"wnn.umap")
#'
#' @export
sc_umap_panel <- function(
  so,
  md_list,
  slot1 = "umap",
  pos_leg = "none",
  show_lab = TRUE
) {
  d <- so
  lp <- setNames(
    lapply(
      seq.int(1, length(md_list), 1),
      function(i) {
        sc_umap(
          so = d,
          md_var = md_list[[i]],
          slot1 = slot1,
          pos_leg = pos_leg,
          show_lab = show_lab
        )
      }
    ),
    md_list
  )
  # Combine output
  if(length(lp) > 1) { # nolint
    d2_out <- ggpubr::ggarrange(
      plotlist = lp,
      ncol = ifelse(length(lp) <= 3, length(lp), 3),
      nrow = ifelse(length(lp) > 3, ceiling(length(lp) / 3), 1)
    )
  }
  if(length(lp) == 1) { # nolint
    d2_out <- lp[[1]]
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
#' # p_umap <- sc_vio_panel_gene(
#' #   so = d, md_vars = c("CellType"), g_name = "SFTPC"
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