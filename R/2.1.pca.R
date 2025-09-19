#' Seurat PCA
#'
#' Performs PCA on a Seurat object.
#'
#' @param so An object on class Seurat.
#' @param slot1 Character string indicating which assay to use for PCA.
#' @return A Seurat object containing PCA results.
#' @examples
#'
#' # d_integrated <- sc_pca(d_integrated, "RNA")
#'
#' @export
sc_pca <- function(
  so,
  slot1
) {
  d <- so
  # Change assay, scale, and run PCA
  if(length(SeuratObject::Layers(d)) < 2) { # nolint
    d <- Seurat::NormalizeData(d)
    d <- Seurat::FindVariableFeatures(d)
  }
  Seurat::DefaultAssay(
    object = d
  ) <- slot1
  d <- Seurat::RunPCA(
    object = Seurat::ScaleData(
      object = d,
      verbose = TRUE
    ),
    verbose = TRUE
  )
  # Visualize dim loadings and elbow plot
  Seurat::VizDimLoadings(
    object = d,
    dims = 1:2,
    reduction = "pca"
  )
  Seurat::ElbowPlot(
    d,
    reduction = "pca",
    ndims = 50
  )
  return(d) # nolint
}

#' scRNA-Seq PCA Loadings Plot
#'
#' Generates a panel of loadings plots given
#' a Seurat object containing PCA results.
#'
#' @param so An object of class Seurat.
#' @param md_list A vector of character strings indicating
#' metadata columns for overlaying on a loadings plot.
#' @param red1 Reduction to use for plotting PCA components.
#' @return A series of loadings plots with specified metadata overlays.
#' @examples
#'
#' # p_pca <- sc_pca_plot(d_integrated,c("col1","col2","col3"), "pca_cor")
#'
#' @export
sc_pca_plot <- function(so, md_list, red1) { # nolint
  # Format input data
  d <- so
  d2 <- data.frame(
    d@meta.data,
    PC1 = d@reductions[[red1]]@cell.embeddings[, 1],
    PC2 = d@reductions[[red1]]@cell.embeddings[, 2]
  )
  # Generate df for each grouping
  d2_list <- setNames(
    lapply(
      c(
        md_list
      ),
      function(x) {
        d2[, c(
          x,
          "PC1",
          "PC2"
        )
        ]
      }
    ),
    c(
      md_list
    )
  )
  # Generate plots
  d2_plot <- lapply(
    c(
      md_list
    ),
    function(x) {
      ggplot2::ggplot(
        d2_list[[x]],
        ggplot2::aes(
          x = PC1, # nolint
          y = PC2, # nolint
          color = .data[[x]] # nolint
        )
      ) +
        ggplot2::geom_point(
          shape = 16,
          size = 2,
          alpha = 0.5
        ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col_univ() # nolint
        ) +
        sc_theme1() # nolint
    }
  )
  # Combine output
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

  return(d2_out)
}
