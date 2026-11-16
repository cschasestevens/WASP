#' General Heatmap
#'
#' Generates a heatmap from a Seurat Object,
#' which can be filtered based on differentially
#' expressed genes.
#'
#' @param so A Seurat object.
#' @param deg_plot Generate heatmap based on DGEA/DA results instead
#' of expression/activity scores from a Seurat object. Supports DGEA/DA
#' result tables generated from sc_diff.
#' @param deg_sig Should only genes with statistically significant
#' differential expression/activity be plotted?
#' @param filt_deg If TRUE, filters genes from a Seurat object based
#' on differential expression in at least one group/cluster.
#' @param top_deg Top DEGs to plot if filtering data by DEGs.
#' @param l_cstm (optional) Filter data based on custom gene list.
#' @param asy Assay to use (ex. "RNA").
#' @param cl_var Character string containing the name
#' of the cluster variable for cell type predictions.
#' @param var_cl Custom gene list variable for heatmap annotation.
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param cl_c Cluster columns?
#' @param cl_r Cluster rows?
#' @param rc Show column names?
#' @param rn Show row names?
#' @param rot_c Rotation of column names.
#' @param col1 Gradient color scheme to use.
#' @param mark_dir Directory for saving the calculated marker gene table.
#' @param mark_tab Path to an existing marker table for skipping the
#' calculation of marker genes if l_cstm is not provided.
#' @param transp Transpose heatmap input matrix (FALSE by default).
#' @param var_as_factor Should cluster variable be converted to factor?
#' @param ret_input Returns input heatmap matrix instead of the heatmap plot.
#' @param hmlinewidth Numeric value for heatmap cell border width.
#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' # p_hmap <- sc_top10_marker_heatmap(
#' #   so = d_recluster[["data"]],
#' #   cl_var = "recluster"
#' # )
#'
#' @export
sc_heatmap <- function(
  so = NULL,
  deg_plot = NULL,
  deg_sig = FALSE,
  filt_deg = TRUE,
  top_deg = 10,
  l_cstm = NULL,
  asy = "SCT",
  cl_var = "CellType",
  var_cl = NULL,
  h_w = 40,
  h_h = 20,
  fs_c = 6,
  fs_r = 8,
  cl_c = FALSE,
  cl_r = FALSE,
  rot_c = 45,
  rc = TRUE,
  rn = TRUE,
  col1 = col_grad(scm = 5), # nolint
  mark_dir = "analyze/",
  mark_tab = NULL,
  transp = FALSE,
  var_as_factor = FALSE,
  ret_input = FALSE,
  hmlinewidth = 0.0
) {
  if (is.null(deg_plot)) {
    d <- so
    if (is.null(l_cstm)) {
      if (is.null(mark_tab)) {
        if (asy == "RNA" || asy == "sct" || asy == "SCT") {
          print("Calculating marker genes for each cluster...")
          Seurat::Idents(d) <- cl_var
          Seurat::DefaultAssay(d) <- asy
          if(length(SeuratObject::Layers(d)) < 2) { # nolint
            d <- Seurat::NormalizeData(d)
          }
          cl_mark <- Seurat::PrepSCTFindMarkers(d)
          cl_mark <- Seurat::FindAllMarkers(cl_mark, verbose = TRUE)
        }
        if(asy == "ATAC") { # nolint
          print(
            "Calculating marker motifs for each cluster..."
          )
          Seurat::Idents(d) <- cl_var
          Seurat::DefaultAssay(d) <- asy
          cl_mark <- Seurat::FindAllMarkers(
            d,
            min.pct = 0.05,
            verbose = TRUE
          )
          names_motif <- data.frame(
            "gene" = seq.int(1, nrow(d@assays$ATAC@meta.features), 1),
            "near.gene" = paste(
              d@assays$ATAC@meta.features[["nearestGene"]],
              seq.int(1, nrow(d@assays$ATAC@meta.features), 1),
              sep = "."
            ),
            "motif" = paste(
              d@assays$ATAC@meta.features[["seqnames"]],
              paste(
                d@assays$ATAC@meta.features[["start"]],
                d@assays$ATAC@meta.features[["end"]],
                sep = "-"
              ),
              sep = ":"
            )
          )
          cl_mark <- dplyr::left_join(
            cl_mark,
            names_motif,
            by = "gene"
          )
        }
        ## Marker gene input matrix (top10 per cell type)
        if(class(cl_mark[["cluster"]]) == "character") { # nolint
          cl_mark <- cl_mark[gtools::mixedorder(cl_mark[["cluster"]]), ]
        }
        cl_mark[["CellType.no"]] <- cl_mark[["cluster"]]
        cl_mark <- dplyr::group_by(
          cl_mark,
          .data[["CellType.no"]] # nolint
        )
        write.table(
          cl_mark,
          paste(mark_dir, "markers.txt", sep = ""),
          sep = "\t",
          row.names = FALSE
        )
      }
      if (!is.null(mark_tab)) {
        cl_mark <- read.table(
          mark_tab,
          header = TRUE,
          sep = "\t"
        )
        cl_mark[["cluster"]] <- factor(
          as.character(cl_mark[["cluster"]]),
          levels = gtools::mixedsort(
            unique(as.character(cl_mark[["cluster"]]))
          )
        )
      }
      ### DEGs per cluster
      if (filt_deg == TRUE) {
        cl_mark <- cl_mark[cl_mark[["p_val_adj"]] < 0.05, ]
        cl_mark <- dplyr::group_by(
          cl_mark,
          .data[["cluster"]] # nolint
        )
        cl_mark <- dplyr::slice_max(
          cl_mark,
          order_by = .data[["avg_log2FC"]], # nolint
          n = top_deg
        )
      }
    }
    Seurat::DefaultAssay(d) <- asy
    ### Subset seurat and scale
    if (is.null(l_cstm)) {
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          unique(cl_mark[["gene"]])
        )
      )
    }
    if (!is.null(l_cstm)) {
      if (length(l_cstm[["Gene"]]) > 1) {
        h <- SeuratObject::FetchData(
          d,
          vars = c(
            cl_var,
            unique(l_cstm[["Gene"]][l_cstm[["Gene"]] %in% rownames(d)])
          )
        )
      }
      if (length(l_cstm[["Gene"]]) == 1) {
        h <- SeuratObject::FetchData(
          d,
          vars = c(
            cl_var,
            l_cstm
          )
        )
      }
    }
  }
  if (!is.null(deg_plot)) {
    h <- deg_plot
  }
  if (var_as_factor == TRUE) {
    h[[1]] <- factor(
      as.character(h[[1]]),
      levels = gtools::mixedsort(unique(as.character(h[[1]])))
    )
  }
  ### Scale and plot average expression/accessibility per cell type
  if (is.null(deg_plot)) {
    h_in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[, 2:ncol(
                  h
                )],
                function(x) {
                  dplyr::select(
                    aggregate(
                      x,
                      list(
                        h[, 1]
                      ),
                      FUN = mean
                    ),
                    c(2)
                  )
                }
              )
            ),
            names(h[, 2:ncol(h)])
          ),
          levels(h[, 1])
        )
      ),
      center = TRUE
    )
  }
  if (!is.null(deg_plot)) {
    # Format DGEA/DA results as matrix
    # and include only significant genes
    if (deg_sig == TRUE) {
      h_in <- h[h[["H.qval"]] < 0.05, ]
    }
    if (deg_sig == FALSE) {
      h_in <- h
    }
    # reassign factor levels
    if (var_as_factor == TRUE) {
      h_in[["CellType"]] <- factor(
        as.character(h_in[["CellType"]]),
        levels = gtools::mixedsort(
          unique(as.character(h_in[["CellType"]]))
        )
      )
    }
    # Create matrix
    h_in <- reshape2::dcast(
      h_in[, c("CellType", "GENE", "log2FC")],
      CellType ~ GENE,
      value.var = "log2FC"
    )
    h_in <- as.matrix(
      magrittr::set_rownames(
        setNames(
          h_in[, 2:ncol(h_in)],
          names(h_in[, 2:ncol(h_in)])
        ),
        levels(h_in[[1]])
      )
    )
  }
  if (!is.null(deg_plot)) {
    h_in[is.na(h_in)] <- 0
  }
  qs <- quantile(
    h_in,
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )
  if (is.null(deg_plot)) {
    h_in <- as.matrix(
      as.data.frame(h_in)[, unlist(
        lapply(
          seq.int(1, ncol(as.data.frame(h_in)), 1),
          function(x) {
            !anyNA(as.data.frame(h_in)[x])
          }
        )
      )
      ]
    )
  }
  if (transp == TRUE) {
    h_in <- t(h_in)
  }
  if (qs[[1]] < 0) {
    fun_hm_col <- circlize::colorRamp2(
      c(
        qs[[1]],
        0,
        qs[[2]]
      ),
      colors = col1
    )
  }
  if (qs[[1]] > 0) {
    fun_hm_col <- circlize::colorRamp2(
      c(
        qs[[1]],
        round(mean(c(qs[[1]], qs[[2]])), digits = 0),
        qs[[2]]
      ),
      colors = col_grad(scm = 4) # nolint
    )
  }
  #---- Include heatmap annotation ----
  if (!is.null(var_cl)) {
    if (is.null(l_cstm)) {
      print("Gene annotation only compatible with custom gene lists. Please provide a formatted list using the l_cstm parameter...") # nolint
    }
    if (!is.null(l_cstm)) {
      if (transp == FALSE) {
        cl2 <- as.factor(
          unlist(
            lapply(
              seq.int(1, length(colnames(h_in)), 1),
              function(i) {
                unique(
                  l_cstm[l_cstm[["Gene"]] == colnames(h_in)[[i]], ][[var_cl]]
                )[[1]]
              }
            )
          )
        )
        # Heatmap annotation colors
        fun_hm_bar <- list(
          Set = setNames(
            col_univ()[1:length(unique(cl2))], # nolint
            gtools::mixedsort(unique(cl2))
          )
        )
        # Annotations
        hm_anno_col <- ComplexHeatmap::HeatmapAnnotation( # nolint
          `Set` = cl2, # nolint
          col = fun_hm_bar,
          show_annotation_name = FALSE,
          show_legend = TRUE
        )
      }
      if (transp == TRUE) {
        cl2 <- as.factor(
          unlist(
            lapply(
              seq.int(1, length(rownames(h_in)), 1),
              function(i) {
                unique(
                  l_cstm[l_cstm[["Gene"]] == rownames(h_in)[[i]], ][[var_cl]]
                )[[1]]
              }
            )
          )
        )
        fun_hm_bar_row <- list(
          Set = setNames(
            col_univ()[1:length(unique(cl2))], # nolint
            gtools::mixedsort(unique(cl2))
          )
        )
        ### Annotations
        hm_anno_row <- ComplexHeatmap::rowAnnotation( # nolint
          `Set` = cl2, # nolint
          col = fun_hm_bar_row,
          show_annotation_name = FALSE,
          show_legend = TRUE
        )
      }
    }
  }
  # Create Plot
  h_out <- ComplexHeatmap::Heatmap(
    h_in,
    top_annotation = if (!is.null(var_cl) && transp == FALSE) hm_anno_col else NULL, # nolint
    left_annotation = if (!is.null(var_cl) && transp == TRUE) hm_anno_row else NULL, # nolint
    col = fun_hm_col,
    name = "Scaled Expression",
    show_column_names = rc,
    show_row_names = rn,
    heatmap_width = ggplot2::unit(h_w, "cm"),
    heatmap_height = ggplot2::unit(h_h, "cm"),
    column_names_rot = rot_c,
    column_names_gp = grid::gpar(fontsize = fs_c),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = fs_r),
    cluster_columns = cl_c,
    cluster_rows = cl_r,
    rect_gp = grid::gpar(col = "black", lwd = hmlinewidth)
  )
  if (ret_input == FALSE) {
    return(h_out) # nolint
  }
  if (ret_input == TRUE) {
    return(h_in) # nolint
  }
}

#' Dot Plot
#'
#' Generates a dot plot from a Seurat Object and a DGEA/DA list.
#'
#' @param l_gene A vector of genes/TFs for plotting. If a data frame
#' is provided, the first column is used for selecting genes to plot.
#' @param so A Seurat object.
#' @param asy1 Assay to use for the selected Seurat object. Note
#' that the assay must match the names in the seleceted
#' list type (ex. DGEA/DA).
#' @param ct_filt Filter based on cell type?
#' @param ct A pattern provided as a character string for
#' matching a specific cell type or types if ct_filt is TRUE.
#' @param ct_var Character string indicating the name of the
#' cluster/cell type variable.
#' @param list_var A vector of character strings containing
#' up to two group variables for stratifying cell types.
#' @param col1 Gradient color scheme to use.
#' @param include_all_genes Logical indicating whether genes/TFs
#' missing from the Seurat object should be included as blank columns
#' in the plot. Useful for comparing a common gene list across data sets.
#' @param groups_xaxis Grouping variable derived from the input gene list
#' that is used to highlight functionally related genes/TFs in the plot.
#' @param groups_yaxis A vector containing one or more integers that
#' are used for plotting dotted lines along the y-axis.
#' Used to separate cell types or related cell groups.
#' @return A dotplot and formatted data frame for the selected assay.
#' @examples
#'
#' # sc_dotplot(
#' #   so = d1,
#' #   l_gene = lgene
#' # )
#'
#' @export
sc_dotplot <- function( # nolint
  l_gene,
  so,
  asy1 = "RNA",
  ct_filt = FALSE,
  ct = NULL,
  ct_var = "CellType",
  list_var = NULL,
  col1 = col_grad(scm = 4), # nolint
  include_all_genes = TRUE,
  groups_xaxis = NULL,
  groups_yaxis = NULL
) {
  # Load data and define plot variables
  d <- so # data
  SeuratObject::DefaultAssay(d) <- asy1 # assay
  c <- ct # cell type column
  l1 <- l_gene # gene list
  varct <- ct_var
  vargrp <- list_var
  # GEX assays
  if (grepl("RNA|sct|SCT", asy1)) {
    gvar <- "GENE"
    exvar <- "avg.exp"
    exvar2 <- "Average Expression"
  }
  # ATAC assays
  if (grepl("ufy.peaks|chromvar", asy1)) {
    gvar <- "ID"
    exvar <- "avg.activity"
    exvar2 <- "Average Activity"
  }
  #---- format data based on custom gene list ----
  # Plot data based on a custom list
  # check for presence/absence of genes in list
  if (class(l1) == "data.frame") {
    top10 <- as.vector(l1[[1]])
  } else {
    top10 <- as.vector(l1)
  }
  top10_pres <- unique(top10[top10 %in% SeuratObject::Features(d)])
  top10_abs <- subset(top10, !(top10 %in% SeuratObject::Features(d)))
  cat("The following are not present in the Seurat object:", "\n")
  cat(top10_abs, "\n")
  # pull genes from Seurat object
  if (!is.null(vargrp)) {
    d1 <- cbind(
      SeuratObject::FetchData(
        d,
        vars = c(
          varct,
          vargrp,
          top10_pres
        )
      )
    )
  } else {
    d1 <- cbind(
      SeuratObject::FetchData(
        d,
        vars = c(
          varct,
          top10_pres
        )
      )
    )
  }
  # Initialize empty columns for missing genes if desired
  if (include_all_genes == TRUE && length(top10_abs) >= 1) {
    cat("Initializing columns for missing genes/TFs", "\n")
    d1 <- dplyr::bind_cols(
      d1,
      setNames(as.data.frame(
        lapply(
          seq.int(1, length(top10_abs), 1),
          function(i) {
            rep(0, nrow(d1))
          }
        )
      ), top10_abs)
    )
    d1 <- d1[, c(varct, vargrp, top10)]
  } else {
    d1 <- d1
  }
  #---- Cell Type Filter ----
  if(ct_filt == TRUE) { # nolint
    d1 <- d1[grepl(c, d1[[varct]]), ]
  }
  #---- Calculate percent expressing cells and ratios ----
  ## Count/ratio table for creating dot plots
  d1_prc <- dplyr::bind_rows(
    setNames(
      lapply(
        top10,
        function(x) {
          # Determine average expression of each gene
          d_avg <- setNames(
            aggregate(
              d1[[x]],
              by = lapply(
                seq.int(1, length(c(varct, vargrp)), 1),
                function(j) {
                  d1[, c(varct, vargrp)[[j]]]
                }
              ),
              function(y) mean(y)
            ),
            c(
              c(varct, vargrp),
              exvar
            )
          )
          d_avg[[exvar]] <- round(
            d_avg[[exvar]],
            digits = 2
          )
          ### percent expressed
          if (is.null(vargrp)) {
            d_prc <- dplyr::count(
              d1, .data[[varct]], # nolint
              .data[[x]] > 0
            )
          } else {
            if (length(vargrp) == 2) {
              d_prc <- dplyr::count(
                d1, .data[[varct]], # nolint
                .data[[vargrp[[1]]]],
                .data[[vargrp[[2]]]],
                .data[[x]] > 0
              )
            } else {
              d_prc <- dplyr::count(
                d1, .data[[varct]], # nolint
                .data[[vargrp[[1]]]],
                .data[[x]] > 0
              )
            }
          }
          d_prc <- setNames(
            dplyr::filter(
              d_prc, d_prc[[length(c(varct, vargrp)) + 1]] == TRUE
            ),
            c(
              c(varct, vargrp),
              "pres", "n"
            )
          )
          if (is.null(vargrp)) {
            d_cnt <- setNames(
              dplyr::count(
                d1, .data[[varct]], # nolint
              ),
              c(
                c(varct, vargrp),
                "n"
              )
            )
          } else {
            if (length(vargrp) == 2) {
              d_cnt <- setNames(
                dplyr::count(
                  d1, .data[[varct]], # nolint
                  .data[[vargrp[[1]]]],
                  .data[[vargrp[[2]]]]
                ),
                c(
                  c(varct, vargrp),
                  "n"
                )
              )
            } else {
              d_cnt <- setNames(
                dplyr::count(
                  d1, .data[[varct]], # nolint
                  .data[[vargrp[[1]]]]
                ),
                c(
                  c(varct, vargrp),
                  "n"
                )
              )
            }
          }
          # combine percents and counts to determine ratio expressed
          d_comb <- dplyr::left_join(
            d_cnt,
            d_prc,
            by = c(varct, vargrp)
          )
          d_comb[is.na(d_comb)] <- 0
          d_comb <- dplyr::mutate(
            d_comb,
            "perc.exp" = round(
              d_comb[["n.y"]] /
                d_comb[["n.x"]],
              digits = 2
            )
          )
          ### combined
          d_comb_out <- dplyr::left_join(
            d_avg,
            d_comb[,
                   c(
                     c(varct, vargrp),
                     "perc.exp"
                   )],
            by = c(varct, vargrp)
          )
          return(d_comb_out) # nolint
        }
      ),
      c(top10)
    ),
    .id = gvar
  )

  #---- Format labels and variable column for plotting ----
  d1_prc <- data.frame(
    d1_prc,
    "labs" = factor(
      do.call(paste, c(d1_prc[c(varct, vargrp)], sep = " ")),
      levels = gtools::mixedsort(
        unique(
          do.call(paste, c(d1_prc[c(varct, vargrp)], sep = " "))
        )
      )
    )
  )
  d1_prc[[gvar]] <- factor(
    d1_prc[[gvar]],
    levels = unique(
      d1_prc[[gvar]]
    )
  )
  #---- Plot annotations ----
  if (!is.null(groups_xaxis)) {
    if (class(l1[[groups_xaxis]]) == "character") {
      l1[[groups_xaxis]] <- factor(
        l1[[groups_xaxis]], levels = unique(l1[[groups_xaxis]])
      )
    }
    annox <- dplyr::count(
      l1,
      .data[[groups_xaxis]]
    )
    annox <- lapply(
      seq.int(1, length(unique(annox[[groups_xaxis]])), 1),
      function(k) {
        n1 <- k
        prop1 <- data.frame(
          "tot_genes" = length(unique(d1_prc[[gvar]])),
          "tot_sections" = nrow(annox),
          "section" = n1,
          "name" = annox[n1, ][[groups_xaxis]],
          "span" = annox[n1, ][["n"]]
        )
        prop1[["min"]] <- sum(annox[1:n1, ][["n"]]) + 0.5 - prop1[["span"]]
        prop1[["max"]] <- sum(annox[1:n1, ][["n"]]) + 0.5
        list(
          # Background shading for groups
          ggplot2::annotate(
            geom = "rect",
            xmin = prop1[["min"]],
            xmax = prop1[["max"]],
            ymin = -Inf, ymax = Inf,
            alpha = 0.25,
            fill = col_univ()[[n1]]
          ),
          # Labels
          ggplot2::annotate(
            geom = "label", x = (prop1[["max"]] + prop1[["min"]]) / 2,
            y = 0.6,
            label = annox[n1, ][[groups_xaxis]],
            color = col_univ()[[n1]],
            fontface = "bold"
          )
        )
      }
    )
  }
  if (!is.null(groups_yaxis)) {
    # add lines to at each position specified by groups_yaxis
    annoy <- lapply(
      seq.int(1, length(groups_yaxis), 1),
      function(k) {
        # add hline to distinguish related cell types
        ggplot2::geom_hline(
          yintercept = groups_yaxis[[k]], linetype = "dashed"
        )
      }
    )
  }
  #---- Final plot output ----
  p_dot <- ggplot2::ggplot(
    data = d1_prc,
    ggplot2::aes(
      x = .data[[gvar]],
      y = .data[["labs"]]
    )
  )
  p_elements <- list(
    "theme_main" = sc_theme1(),
    "theme_extra" =  ggplot2::theme(
      plot.margin = ggplot2::unit(
        c(.2, .2,
          .2, .2),
        "cm"
      )
    ),
    "coords" =  ggplot2::geom_point(
      shape = 21,
      size = 0
    ),
    "group_shading" =  if (exists("annox")) annox else NULL,
    "celltype_lines" = if (exists("annoy")) annoy else NULL,
    "points" =  ggplot2::geom_point(
      shape = 21,
      mapping = ggplot2::aes(
        x = .data[[gvar]],
        y = .data[["labs"]],
        fill = .data[[exvar]],
        size = .data[["perc.exp"]]
      )
    ),
    "scale_gradient" = ggplot2::scale_fill_gradientn(
      colors = col1 # nolint
    ),
    "scale_size" = ggplot2::scale_size_area(max_size = 12),
    "labels" = ggplot2::labs(
      fill = exvar2,
      size = "Percent Expressed",
      y = ""
    )
  )
  pout <- p_dot + p_elements
  return( # nolint
    list(
      "Input" = d1_prc,
      "Plot" = pout
    )
  )
}
