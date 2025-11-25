#' General Heatmap
#'
#' Generates a heatmap from a Seurat Object,
#' which can be filtered based on differentially
#' expressed genes.
#'
#' @param so A Seurat object.
#' @param title1 Heatmap title.
#' @param filt_deg If TRUE, filters genes from a Seurat object based
#' on differential expression in at least one group/cluster.
#' @param l_cstm (optional) Filter data based on custom gene list.
#' @param asy Assay to use (ex. "RNA").
#' @param cl_var Character string containing the name
#' of the cluster variable for cell type predictions.
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param cl_c Cluster columns?
#' @param cl_r Cluster rows?
#' @param rc Show column names?
#' @param rn Show row names?
#' @param rot_c Rotation of column names.
#' @param col1 Gradient color scheme to use
#' (must be exactly 4 colors in length).
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
  so,
  filt_deg = TRUE,
  top_deg = 10,
  l_cstm = NULL,
  asy = "SCT",
  cl_var = "CellType",
  h_w = 40,
  h_h = 20,
  fs_c = 6,
  fs_r = 8,
  cl_c = FALSE,
  cl_r = FALSE,
  rot_c = 45,
  rc = TRUE,
  rn = TRUE,
  col1 = col_grad(scm = 3), # nolint
  mark_dir = "analyze/"
) {
  d <- so
  if (is.null(l_cstm)) {
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
    ### DEGs per cluster
    if (filt_deg == TRUE) {
      cl_mark <- cl_mark[cl_mark[["p_val_adj"]] < 0.05, ]
      cl_mark <- dplyr::group_by(
        cl_mark,
        .data[["cluster"]] # nolint
      )
      cl_mark <- dplyr::slice_max(
        cl_mark,
        order_by = .data[["avg_log2FC"]],
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
    h <- SeuratObject::FetchData(
      d,
      vars = c(
        cl_var,
        unique(l_cstm[l_cstm %in% rownames(d)])
      )
    )
  }
  h[[1]] <- factor(
    as.character(h[[1]]),
    levels = gtools::mixedsort(unique(as.character(h[[1]])))
  )
  ### Scale and plot average expression/accessibility per cell type
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
  qs <- quantile(
    h_in,
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )
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
  fun_hm_col <- circlize::colorRamp2(
    c(
      qs[[1]],
      (qs[[1]]) / 2,
      (qs[[2]]) / 2,
      qs[[2]]
    ),
    colors = col1
  )
  # Create Plot
  h_out <- ComplexHeatmap::Heatmap(
    h_in,
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
    cluster_rows = cl_r
  )
  return(h_out)
}

#' Dot Plot
#'
#' Generates a dot plot from a Seurat Object based on a DGEA/DA list based on
#' a threshold or custom gene list for a specific cell type and comparison.
#'
#' @param p_type Should a custom gene list or DGEA results object
#' be used for plotting?
#' Type "cstm" for custom gene lists and "threshold"
#' for using DGEA/DA results.
#' @param topn If p_type is "threshold," select the number of
#' top DEGs/TFs to plot.
#' @param l_gene A list of DGEA/DA results
#' or a vector of selected genes for plotting.
#' @param filt_var If plotting based on a threshold,
#' which variable should be used for filtering?
#' @param so An object of class Seurat.
#' @param asy1 Assay to use for the selected Seurat object. Note
#' that the assay must match the names in the seleceted
#' list type (ex. DGEA/DA).
#' @param ct_filt Filter based on cell type?
#' @param ct A pattern provided as a character string for
#' matching a specific cell type or types if ct_filt is TRUE.
#' @param ct_var Character string indicating the name of the
#' cluster/cell type variable.
#' @param split_var Logical indicating whether the cluster variable
#' should be stratified by additional group variables.
#' @param list_var If split_var is TRUE, provide a vector of character strings
#' containing two group variables for stratifying plot points.
#' @param col1 Gradient color scheme to use.
#' @param vline1 Add a vertical line to plot?
#' @param xint If vline1 is TRUE, specify the x-intercept of the vertical line
#' to place on the dot plot.
#' @return A dotplot and formatted data frame for the selected assay.
#' @examples
#'
#' # dot1 <- sc_dotplot( # nolint
#' #   p_type = "cstm",
#' #   l_gene = c("IRF1", "IRF2", "STAT1"),
#' #   so = d,
#' #   asy1 = "SCT",
#' #   ct_var = "CellType",
#' #   split_var = TRUE,
#' #   list_var = "Group",
#' #   vline1 = FALSE
#' # )
#'
#' @export
sc_dotplot <- function( # nolint
  p_type = "threshold",
  topn = 10,
  l_gene,
  filt_var = NULL,
  so,
  asy1 = "SCT",
  ct_filt = FALSE,
  ct = NULL,
  ct_var,
  split_var = FALSE,
  list_var = NULL,
  col1 = col_grad(), # nolint
  vline1 = FALSE,
  xint = NULL
) {
  # Load data and define plot variables
  d <- so
  SeuratObject::DefaultAssay(d) <- asy1
  c <- ct
  l1 <- l_gene
  if(grepl("RNA|sct|SCT", asy1)) { # nolint
    gvar <- "GENE"
    exvar <- "avg.exp"
    exvar2 <- "Average Expression"
  }
  if(grepl("ufy.peaks|chromvar", asy1)) { # nolint
    gvar <- "ID"
    exvar <- "avg.activity"
    exvar2 <- "Average Activity"
    if(p_type == "cstm") { # nolint
      list_tf <- unlist(
        lapply(
          row.names(d),
          function(x) {
            name(TFBSTools::getMatrixByID(JASPAR2020, ID = x)) # nolint
          }
        )
      )
      list_tf <- data.frame(
        "ID" = rownames(d),
        "Name" = list_tf
      )
      list_tf <- list_tf[order(list_tf[["Name"]]), ]
      l1 <- list_tf[list_tf[["Name"]] %in% l1, "ID"]
      l2 <- unlist(
        lapply(
          l1,
          function(x) {
            name(TFBSTools::getMatrixByID(JASPAR2020, ID = x)) # nolint
          }
        )
      )
    }
  }
  # Select data based on threshold
  if(p_type == "threshold") { # nolint
    ## Return specified genes/motifs for selected comparison and cell types
    top10 <- dplyr::select(
      dplyr::slice_max(
        dplyr::group_by(
          l1,
          dplyr::across(
            dplyr::all_of(
              c(
                "Comparison",
                ct_var
              )
            )
          )
        ),
        order_by = -.data[[filt_var]], # nolint
        n = topn
      ),
      c(
        ct_var,
        "Comparison",
        gvar,
        filt_var
      )
    )
    ## filter list for chosen cell types
    top10 <- top10[grepl(c, top10[[ct_var]]), ]
    ## return list of unique genes
    top10_pres <- unique(
      top10[[gvar]][top10[[gvar]] %in% SeuratObject::Features(d)]
    )
    top10_abs <- subset(
      top10[[gvar]],
      !(top10[[gvar]] %in% SeuratObject::Features(d))
    )
  }
  # Plot data based on a custom list
  if(p_type == "cstm") { # nolint
    top10 <- as.vector(l1)
    top10_pres <- unique(top10[top10 %in% SeuratObject::Features(d)])
    top10_abs <- subset(top10, !(top10 %in% SeuratObject::Features(d)))
  }
  ## select genes from Seurat object
  if(split_var == TRUE) { # nolint
    d1 <- cbind(
      SeuratObject::FetchData(
        d,
        vars = c(
          ct_var,
          list_var,
          top10_pres
        )
      )
    )
    ## Set column names if plotting activity scores
    if(grepl("ufy.peaks|chromvar", asy1)) { # nolint
      if(p_type == "cstm") { # nolint
        d1 <- setNames(
          d1,
          c(ct_var, list_var, l2)
        )
        top10_pres <- l2
      }
      if(p_type == "threshold") { # nolint
        l2 <- unlist(
          lapply(
            top10_pres,
            function(x) {
              name(TFBSTools::getMatrixByID(JASPAR2020, ID = x)) # nolint
            }
          )
        )
        d1 <- setNames(
          d1,
          c(ct_var, list_var, l2)
        )
        top10_pres <- l2
      }
    }
  }
  if(split_var == FALSE) { # nolint
    d1 <- cbind(
      SeuratObject::FetchData(
        d,
        vars = c(
          ct_var,
          top10_pres
        )
      )
    )
    ## Set column names if plotting activity scores
    if(grepl("ufy.peaks|chromvar", asy1)) { # nolint
      d1 <- setNames(
        d1,
        c(ct_var, l2)
      )
      top10_pres <- l2
    }
  }
  ## Subset based on cell type if ct_filt is TRUE
  if(ct_filt == TRUE) { # nolint
    d1 <- d1[grepl(c, d1[[ct_var]]), ]
  }
  ## Count/ratio table for creating dot plots
  d1_prc <- dplyr::bind_rows(
    setNames(
      lapply(
        top10_pres,
        function(x) {
          # Determine average expression of each gene
          ## for 2 variables:
          if(split_var == TRUE && length(list_var) == 2) { # nolint
            ### average expression
            d_avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[,  c(ct_var)],
                  d1[,  c(list_var[[1]])],
                  d1[,  c(list_var[[2]])]
                ),
                function(y) mean(y)
              ),
              c(
                ct_var, c(list_var),
                exvar
              )
            )
            d_avg[[exvar]] <- round(
              d_avg[[exvar]],
              digits = 2
            )
            ### percent expressed
            d_prc <- dplyr::count(
              d1, .data[[ct_var]], # nolint
              .data[[list_var[[1]]]],
              .data[[list_var[[2]]]],
              .data[[x]] > 0
            )
            d_prc <- setNames(
              dplyr::filter(
                d_prc, d_prc[4] == TRUE
              ),
              c(
                ct_var, c(list_var),
                "pres", "n"
              )
            )
            d_cnt <- setNames(
              dplyr::count(
                d1, .data[[ct_var]], # nolint
                .data[[list_var[[1]]]],
                .data[[list_var[[2]]]]
              ),
              c(
                ct_var, c(list_var),
                "n"
              )
            )
            d_comb <- dplyr::left_join(
              d_cnt,
              d_prc,
              by = c(
                ct_var, c(list_var)
              )
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
                       ct_var, c(list_var),
                       "perc.exp"
                     )],
              by = c(
                ct_var,
                c(list_var)
              )
            )
          }
          if(split_var == TRUE && length(list_var) < 2) { # nolint
            ### average expression
            d_avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[, c(ct_var)],
                  d1[, c(list_var[[1]])]
                ),
                function(y) mean(y)
              ),
              c(
                ct_var, c(list_var),
                exvar
              )
            )
            d_avg[[exvar]] <- round(
              d_avg[[exvar]],
              digits = 2
            )
            ### percent expressed
            d_prc <- dplyr::count(
              d1,.data[[ct_var]], # nolint
              .data[[list_var[[1]]]],
              .data[[x]] > 0
            )
            d_prc <- setNames(
              dplyr::filter(
                d_prc, d_prc[3] == TRUE
              ),
              c(
                ct_var, c(list_var),
                "pres", "n"
              )
            )
            d_cnt <- setNames(
              dplyr::count(
                d1, .data[[ct_var]], # nolint
                .data[[list_var[[1]]]]
              ),
              c(
                ct_var, c(list_var),
                "n"
              )
            )
            d_comb <- dplyr::left_join(
              d_cnt,
              d_prc,
              by = c(
                ct_var, c(list_var)
              )
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
                       ct_var, c(list_var),
                       "perc.exp"
                     )],
              by = c(
                ct_var,
                c(list_var)
              )
            )
          }
          if(split_var == FALSE) { # nolint
            ### average expression/activity
            d_avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[, c(ct_var)]
                ),
                function(y) mean(y)
              ),
              c(
                ct_var,
                exvar
              )
            )
            d_avg[[exvar]] <- round(
              d_avg[[exvar]],
              digits = 2
            )
            ### percent expressed
            d_prc <- dplyr::count(
              d1,.data[[ct_var]], # nolint
              .data[[x]] > 0
            )
            d_prc <- setNames(
              dplyr::filter(
                d_prc, d_prc[2] == TRUE
              ),
              c(
                ct_var,
                "pres", "n"
              )
            )
            d_cnt <- setNames(
              dplyr::count(
                d1, .data[[ct_var]] # nolint
              ),
              c(
                ct_var,
                "n"
              )
            )
            d_comb <- dplyr::left_join(
              d_cnt,
              d_prc,
              by = c(
                ct_var
              )
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
                       ct_var,
                       "perc.exp"
                     )],
              by = c(
                ct_var
              )
            )
          }
          return(d_comb_out) # nolint
        }
      ),
      c(top10_pres)
    ),
    .id = gvar
  )
  ## Add row name labels and convert GENE/Motif column to factor
  if(split_var == TRUE && length(list_var) == 2) { # nolint
    d1_prc <- data.frame(
      d1_prc,
      "labs" = factor(
        paste(
          d1_prc[[ct_var]],
          d1_prc[[list_var[[1]]]],
          d1_prc[[list_var[[2]]]],
          sep = " "
        ),
        levels = gtools::mixedsort(
          unique(
            paste(
              d1_prc[[ct_var]],
              d1_prc[[list_var[[1]]]],
              d1_prc[[list_var[[2]]]],
              sep = " "
            )
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
  }
  ## Add row name labels and convert GENE column to factor
  if(split_var == TRUE && length(list_var) < 2) { # nolint
    d1_prc <- data.frame(
      d1_prc,
      "labs" = factor(
        paste(
          d1_prc[[ct_var]],
          d1_prc[[list_var[[1]]]],
          sep = " "
        ),
        levels = gtools::mixedsort(
          unique(
            paste(
              d1_prc[[ct_var]],
              d1_prc[[list_var[[1]]]],
              sep = " "
            )
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
  }
  if(split_var == FALSE) { # nolint
    d1_prc <- data.frame(
      d1_prc,
      "labs" = d1_prc[[ct_var]]
    )
    d1_prc[[gvar]] <- factor(
      d1_prc[[gvar]],
      levels = unique(
        d1_prc[[gvar]]
      )
    )
  }
  ## Plot
  if(vline1 == FALSE) { # nolint
    p_dot <- ggplot2::ggplot(
      d1_prc,
      ggplot2::aes(
        x = .data[[gvar]], # nolint
        y = .data[["labs"]],
        fill = .data[[exvar]],
        size = .data[["perc.exp"]]
      )
    ) +
      ggplot2::geom_point(
        shape = 21
      ) +
      ggplot2::scale_fill_gradientn(
        colors = col1 # nolint
      ) +
      ggplot2::scale_size_area(max_size = 12) +
      sc_theme1() + # nolint
      ggplot2::labs(
        fill = exvar2,
        size = "Percent Expressed",
        y = ""
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(.2, .2,
            .2, .2),
          "cm"
        )
      )
  }
  if(vline1 == TRUE) { # nolint
    p_dot <- ggplot2::ggplot(
      d1_prc,
      ggplot2::aes(
        x = .data[[gvar]], # nolint
        y = .data[["labs"]],
        fill = .data[[exvar]],
        size = .data[["perc.exp"]]
      )
    ) +
      ggplot2::geom_point(
        shape = 21
      ) +
      ggplot2::geom_vline(xintercept = xint, linetype = "dashed") +
      ggplot2::scale_fill_gradientn(
        colors = col1 # nolint
      ) +
      ggplot2::scale_size_area(max_size = 12) +
      sc_theme1() + # nolint
      ggplot2::labs(
        fill = exvar2,
        size = "Percent Expressed",
        y = ""
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(.2, .2,
            .2, .2),
          "cm"
        )
      )
  }
  if(exists("top10.abs") && length(top10_abs) == 1) { # nolint
    print(
      paste(
        top10_abs,
        "was not found in the provided Seurat object;
        points for this gene/motif were excluded from the dot plot...",
        sep = " "
      )
    )
  }
  if(exists("top10.abs") && length(top10_abs) > 1) { # nolint
    print(
      paste(
        top10_abs,
        "were not found in the provided Seurat object;
        these genes/motifs were excluded from the dot plot...",
        sep = " "
      )
    )
  }
  return(
    list(
      "Input" = d1_prc,
      "Plot" = p_dot
    )
  )
}

#' Top Differentially Active Motif Heatmap
#'
#' Generates a heatmap from a Seurat Object and
#' differential activity analysis result.
#'
#' @param so An object of class Seurat.
#' @param mot_m Differential activity results object.
#' @param mot_list (optional) A vector of motif names to plot on the
#' heatmap. Note that supplying a list will override filtering based
#' on the top n motifs for a specified activity results object.
#' @param filt_ct Logical indicating whether the data should
#' be filtered based on cell type
#' @param ct_name Character string containing the cell type name to filter.
#' @param cl_var Character string containing the name
#' of the clustering variable.
#' @param split_var Logical indicating whether the cluster variable
#' should be stratified by additional group variables.
#' @param inc_var Logical indicating whether two separate clustering variables
#' should be plotted in a single figure (useful for comparing summarized and
#' split expression together). Uses list_var if TRUE.
#' @param list_var (Optional) A vector of character strings
#' indicating the name(s) of up to two group variables
#' for stratifying plot points.
#' @param top_n Number of motifs to use (filters top-100 motifs per assay or
#' specified cell type and subsequently filters by scaled activity score).
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param cl_c Cluster columns?
#' @param cl_r Cluster rows?
#' @param col1 Gradient color scheme to use.
#' @return A ComplexHeatmap object containing a top motif heatmap.
#' @examples
#'
#' # tf_heatmap <- sc_top_motif_heatmap(
#' #   # Seurat object
#' #   so = d,
#' #   # Differential activity results
#' #   mot_m = diff_output_activity,
#' #   # Filter based on CellType?
#' #   filt_ct = TRUE,
#' #   # Cell type name
#' #   ct_name = "11.Secretory",
#' #   # Clustering column
#' #   cl_var = "CellType",
#' #   # Split by additional variable(s)?
#' #   split_var = TRUE,
#' #   # Additional variable(s) to split plot (if split_var = TRUE)
#' #   list_var = c("Airway"),
#' #   # Number of motifs to use
#' #   top_n = 50,
#' #   # Heatmap width
#' #   h_w = 36,
#' #   # Heatmap height
#' #   h_h = 12,
#' #   # Column font size
#' #   fs_c = 6,
#' #   # Row font size
#' #   fs_r = 8,
#' #   # Gradient color scheme
#' #   col_grad()[c(3, 6, 9, 12)]
#' # )
#'
#' @export
sc_top_motif_heatmap <- function( # nolint
  so,
  mot_m,
  mot_list = NULL,
  filt_ct = FALSE,
  ct_name = NULL,
  cl_var = NULL,
  split_var = FALSE,
  inc_var = FALSE,
  list_var = NULL,
  top_n = 10,
  ord1 = FALSE,
  ord_c = NULL,
  h_w = 44,
  h_h = 14,
  fs_c = 10,
  fs_r = 10,
  cl_c = TRUE,
  cl_r = TRUE,
  col1 = col_grad # nolint
) {
  d <- so
  Seurat::DefaultAssay(d) <- "chromvar"

  ## Motif input matrix (top motifs per cell type)
  d_mark <- mot_m
  if(class(d_mark[["cluster"]]) == "character") { # nolint
    d_mark <- d_mark[gtools::mixedorder(d_mark[["cluster"]]), ]
  }

  if(!is.null(mot_list)) { # nolint
    if(filt_ct == TRUE) { # nolint
      d_mark <- dplyr::bind_rows(setNames(
        lapply(
          seq.int(1, length(unique(d_mark[["gene"]])), 1),
          function(x) {
            d1 <- d_mark[
              d_mark[["gene"]] == unique(d_mark[["gene"]])[[x]] &
                d_mark[["avg_diff"]] > 0,
            ]
            d1 <- d1[order(d1[["avg_diff"]], decreasing = TRUE), ]
            d1 <- d1[1, ]
            return(d1) # nolint
          }
        ),
        c(unique(d_mark[["gene"]]))
      ))
      d_mark <- d_mark[
        d_mark[["cluster"]] == ct_name &
          !is.na(d_mark[["cluster"]]),
      ]
      d_mark[["ID"]] <- gsub(
        "^([^.]*\\.)|\\..*",
        "\\1",
        d_mark[["ID"]]
      )
      #### Save table
      write.table(
        d_mark,
        paste(
          "analysis/table.diff.active.motifs.",
          ct_name,
          ".txt",
          sep = ""
        ),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t"
      )
    }
    if(filt_ct == FALSE) { # nolint
      d_mark[["ID"]] <- gsub(
        "^([^.]*\\.)|\\..*",
        "\\1",
        d_mark[["ID"]]
      )
      mot_l <- gsub(
        "^([^.]*\\.)|\\..*",
        "\\1",
        mot_list[["ID"]]
      )
    }
    ### Subset seurat and scale
    ## select motifs from Seurat object
    if(split_var == TRUE && ord1 == FALSE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          list_var,
          mot_l
        )
      )
    }
    if(split_var == FALSE && ord1 == FALSE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          mot_l
        )
      )
    }
    if(split_var == TRUE && ord1 == TRUE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          list_var,
          ord_c
        )
      )
    }
    if(split_var == FALSE && ord1 == TRUE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          ord_c
        )
      )
    }
    if(split_var == FALSE && inc_var == TRUE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          list_var,
          mot_l
        )
      )
    }
    if(split_var == TRUE && length(list_var) == 2) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              as.data.frame(
                lapply(
                  h[, 4:ncol(
                    h
                  )
                  ],
                  function(x) {
                    dplyr::select(
                      aggregate(
                        x,
                        list(
                          h[, 1],
                          h[, 2],
                          h[, 3]
                        ),
                        FUN = mean
                      ),
                      c(4)
                    )
                  }
                )
              ),
              names(h[, 4:ncol(h)])
            ),
            paste(
              aggregate(
                h[, 4],
                list(
                  h[, 1],
                  h[, 2],
                  h[, 3]
                ),
                FUN = mean
              )[[1]],
              aggregate(
                h[, 4],
                list(
                  h[, 1],
                  h[, 2],
                  h[, 3]
                ),
                FUN = mean
              )[[2]],
              aggregate(
                h[, 4],
                list(
                  h[, 1],
                  h[, 2],
                  h[, 3]
                ),
                FUN = mean
              )[[3]],
              sep = "."
            )
          )
        ),
        center = TRUE
      )
    }
    if(split_var == TRUE && length(list_var) < 2) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              as.data.frame(
                lapply(
                  h[, 3:ncol(
                    h
                  )
                  ],
                  function(x) {
                    dplyr::select(
                      aggregate(
                        x,
                        list(
                          h[, 1],
                          h[, 2]
                        ),
                        FUN = mean
                      ),
                      c(3)
                    )
                  }
                )
              ),
              names(h[, 3:ncol(h)])
            ),
            paste(
              aggregate(
                h[, 3],
                list(
                  h[, 1],
                  h[, 2]
                ),
                FUN = mean
              )[[1]],
              aggregate(
                h[, 3],
                list(
                  h[, 1],
                  h[, 2]
                ),
                FUN = mean
              )[[2]],
              sep = "."
            )
          )
        ),
        center = TRUE
      )
    }
    if(split_var == FALSE && inc_var == FALSE) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              as.data.frame(
                lapply(
                  h[, 2:ncol(
                    h
                  )
                  ],
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
    if(split_var == FALSE && inc_var == TRUE) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              rbind(
                as.data.frame(
                  lapply(
                    h[, 3:ncol(
                      h
                    )
                    ],
                    function(x) {
                      dplyr::select(
                        aggregate(
                          x,
                          list(
                            h[, 2]
                          ),
                          FUN = mean
                        ),
                        c(2)
                      )
                    }
                  )
                ),
                as.data.frame(
                  lapply(
                    h[, 3:ncol(
                      h
                    )
                    ],
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
                )
              ),
              names(h[, 3:ncol(h)])
            ),
            c(levels(as.factor(h[, 2])), levels(h[, 1]))
          )
        ),
        center = TRUE
      )
    }
    qs <- quantile(
      h_in,
      probs = c(
        0.05,
        0.95
      ),
      na.rm = TRUE
    )
    h_in <- as.matrix(
      as.data.frame(h_in)[, unlist(
        lapply(
          seq.int(1, ncol(as.data.frame(h_in)), 1),
          function(x) !anyNA(as.data.frame(h_in)[x])
        )
      )
      ]
    )
    if(split_var == TRUE && length(list_var) == 2) { # nolint
      h_in <- h_in[
        gtools::mixedsort(
          unique(
            paste(
              h[["CellType"]],
              h[[list_var[[1]]]], # nolint
              h[[list_var[[2]]]],
              sep = "."
            )
          )
        ),
      ]
    }
    if(split_var == TRUE && length(list_var) < 2) { # nolint
      h_in <- h_in[
        gtools::mixedsort(
          unique(
            paste(
              h[["CellType"]],
              h[[list_var[[1]]]], # nolint
              sep = "."
            )
          )
        ),
      ]
    }
    tf_names <- data.frame("ID" = colnames(h_in))
    tf_names <- dplyr::left_join(
      tf_names,
      unique(d_mark[, c("ID", "gene")]),
      by = "ID"
    )
    colnames(h_in) <- tf_names[["gene"]]
    fun_hm_col <- circlize::colorRamp2(
      c(
        qs[[1]],
        (qs[[1]]) / 2,
        (qs[[2]]) / 2,
        qs[[2]]
      ),
      colors = col1
    )
    # Create Plot
    if(filt_ct == TRUE) { # nolint
      # Plot
      h_out <- ComplexHeatmap::Heatmap(
        h_in,
        col = fun_hm_col,
        name = "Activity Score",
        show_column_names = TRUE,
        show_row_names = TRUE,
        heatmap_width = ggplot2::unit(h_w, "cm"),
        heatmap_height = ggplot2::unit(h_h, "cm"),
        column_title = paste("Top Motifs:", ct_name),
        column_names_rot = 45,
        column_names_gp = grid::gpar(fontsize = fs_c),
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = fs_r),
        cluster_columns = cl_c,
        cluster_rows = cl_r
      )
    }
    if(filt_ct == FALSE) { # nolint
      h_out <- ComplexHeatmap::Heatmap(
        h_in,
        col = fun_hm_col,
        name = "Activity Score",
        show_column_names = TRUE,
        show_row_names = TRUE,
        heatmap_width = ggplot2::unit(h_w, "cm"),
        heatmap_height = ggplot2::unit(h_h, "cm"),
        column_title = paste("Selected TF Motifs"),
        column_names_rot = 45,
        column_names_gp = grid::gpar(fontsize = fs_c),
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = fs_r),
        cluster_columns = cl_c,
        cluster_rows = cl_r
      )
    }
  }
  if(is.null(mot_list)) { # nolint
    if(filt_ct == TRUE) { # nolint
      d_mark <- dplyr::bind_rows(setNames(
        lapply(
          seq.int(1, length(unique(d_mark[["gene"]])), 1),
          function(x) {
            d1 <- d_mark[
              d_mark[["gene"]] == unique(d_mark[["gene"]])[[x]] &
                d_mark[["avg_diff"]] > 0,
            ]
            d1 <- d1[order(d1[["avg_diff"]], decreasing = TRUE), ]
            d1 <- d1[1, ]
            return(d1) # nolint
          }
        ),
        c(unique(d_mark[["gene"]]))
      ))
      d_mark <- d_mark[
        d_mark[["cluster"]] == ct_name &
          !is.na(d_mark[["cluster"]]),
      ]
      d_mark[["ID"]] <- gsub(
        "^([^.]*\\.)|\\..*",
        "\\1",
        d_mark[["ID"]]
      )
      #### Save table
      write.table(
        d_mark,
        paste(
          "analysis/table.diff.active.motifs.",
          ct_name,
          ".txt",
          sep = ""
        ),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t"
      )
    }
    if(filt_ct == FALSE) { # nolint
      d_mark <- dplyr::group_by(
        d_mark,
        .data[["cluster"]] # nolint
      )
      ## Top motifs per cluster (by p value then by fold change)
      d_mark <- dplyr::slice_max(
        d_mark[d_mark[["avg_diff"]] > 0 & d_mark[["p_val_adj"]] < 0.05, ],
        order_by = -.data[["p_val_adj"]], # nolint
        n = 100
      )[, c(
        "gene",
        "ID",
        "cluster",
        "avg_diff",
        "p_val_adj"
      )
      ]
      d_mark <- dplyr::group_by(
        d_mark,
        .data[["cluster"]] # nolint
      )
      d_mark <- dplyr::slice_max(
        d_mark,
        order_by = .data[["avg_diff"]], # nolint
        n = top_n
      )[, c(
        "gene",
        "ID",
        "cluster"
      )
      ]
      d_mark[["ID"]] <- gsub(
        "^([^.]*\\.)|\\..*",
        "\\1",
        d_mark[["ID"]]
      )
      #### Save table
      write.table(
        d_mark,
        paste(
          "analysis/table.marker.motifs.Top",
          paste(as.character(top_n)),
          ".txt",
          sep = ""
        ),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t"
      )
    }
    ### Subset seurat and scale
    ## select genes from Seurat object
    if(split_var == TRUE && ord1 == FALSE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          list_var,
          unique(
            d_mark[["ID"]]
          )
        )
      )
    }
    if(split_var == FALSE && ord1 == FALSE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          unique(
            d_mark[["ID"]]
          )
        )
      )
    }
    if(split_var == TRUE && ord1 == TRUE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          list_var,
          ord_c
        )
      )
    }
    if(split_var == FALSE && ord1 == TRUE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          ord_c
        )
      )
    }
    if(split_var == FALSE && inc_var == TRUE) { # nolint
      h <- SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          list_var,
          unique(
            d_mark[["ID"]]
          )
        )
      )
    }
    if(split_var == TRUE && length(list_var) == 2) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              as.data.frame(
                lapply(
                  h[, 4:ncol(
                    h
                  )
                  ],
                  function(x) {
                    dplyr::select(
                      aggregate(
                        x,
                        list(
                          h[, 1],
                          h[, 2],
                          h[, 3]
                        ),
                        FUN = mean
                      ),
                      c(4)
                    )
                  }
                )
              ),
              names(h[, 4:ncol(h)])
            ),
            paste(
              aggregate(
                h[, 4],
                list(
                  h[, 1],
                  h[, 2],
                  h[, 3]
                ),
                FUN = mean
              )[[1]],
              aggregate(
                h[, 4],
                list(
                  h[, 1],
                  h[, 2],
                  h[, 3]
                ),
                FUN = mean
              )[[2]],
              aggregate(
                h[, 4],
                list(
                  h[, 1],
                  h[, 2],
                  h[, 3]
                ),
                FUN = mean
              )[[3]],
              sep = "."
            )
          )
        ),
        center = TRUE
      )
    }
    if(split_var == TRUE && length(list_var) < 2) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              as.data.frame(
                lapply(
                  h[, 3:ncol(
                    h
                  )
                  ],
                  function(x) {
                    dplyr::select(
                      aggregate(
                        x,
                        list(
                          h[, 1],
                          h[, 2]
                        ),
                        FUN = mean
                      ),
                      c(3)
                    )
                  }
                )
              ),
              names(h[, 3:ncol(h)])
            ),
            paste(
              aggregate(
                h[, 3],
                list(
                  h[, 1],
                  h[, 2]
                ),
                FUN = mean
              )[[1]],
              aggregate(
                h[, 3],
                list(
                  h[, 1],
                  h[, 2]
                ),
                FUN = mean
              )[[2]],
              sep = "."
            )
          )
        ),
        center = TRUE
      )
    }
    if(split_var == FALSE && inc_var == FALSE) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              as.data.frame(
                lapply(
                  h[, 2:ncol(
                    h
                  )
                  ],
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
    if(split_var == FALSE && inc_var == TRUE) { # nolint
      h_in <- scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              rbind(
                as.data.frame(
                  lapply(
                    h[, 3:ncol(
                      h
                    )
                    ],
                    function(x) {
                      dplyr::select(
                        aggregate(
                          x,
                          list(
                            h[, 2]
                          ),
                          FUN = mean
                        ),
                        c(2)
                      )
                    }
                  )
                ),
                as.data.frame(
                  lapply(
                    h[, 3:ncol(
                      h
                    )
                    ],
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
                )
              ),
              names(h[, 3:ncol(h)])
            ),
            c(levels(as.factor(h[, 2])), levels(h[, 1]))
          )
        ),
        center = TRUE
      )
    }
    qs <- quantile(
      h_in,
      probs = c(
        0.05,
        0.95
      ),
      na.rm = TRUE
    )
    h_in <- as.matrix(
      as.data.frame(h_in)[, unlist(
        lapply(
          seq.int(1, ncol(as.data.frame(h_in)), 1),
          function(x) !anyNA(as.data.frame(h_in)[x])
        )
      )
      ]
    )
    if(split_var == TRUE && length(list_var) == 2) { # nolint
      h_in <- h_in[
        gtools::mixedsort(
          unique(
            paste(
              h[["CellType"]],
              h[[list_var[[1]]]], # nolint
              h[[list_var[[2]]]],
              sep = "."
            )
          )
        ),
      ]
    }
    if(split_var == TRUE && length(list_var) < 2) { # nolint
      h_in <- h_in[
        gtools::mixedsort(
          unique(
            paste(
              h[["CellType"]],
              h[[list_var[[1]]]], # nolint
              sep = "."
            )
          )
        ),
      ]
    }
    tf_names <- data.frame("ID" = colnames(h_in))
    tf_names <- dplyr::left_join(
      tf_names,
      unique(d_mark[, c("ID", "gene")]),
      by = "ID"
    )
    colnames(h_in) <- tf_names[["gene"]]
    fun_hm_col <- circlize::colorRamp2(
      c(
        qs[[1]],
        (qs[[1]]) / 2,
        (qs[[2]]) / 2,
        qs[[2]]
      ),
      colors = col1
    )
    # Create Plot
    if(filt_ct == TRUE) { # nolint
      # Plot
      h_out <- ComplexHeatmap::Heatmap(
        h_in,
        col = fun_hm_col,
        name = "Activity Score",
        show_column_names = TRUE,
        show_row_names = TRUE,
        heatmap_width = ggplot2::unit(h_w, "cm"),
        heatmap_height = ggplot2::unit(h_h, "cm"),
        column_title = paste("Top Motifs:", ct_name),
        column_names_rot = 45,
        column_names_gp = grid::gpar(fontsize = fs_c),
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = fs_r),
        cluster_columns = cl_c,
        cluster_rows = cl_r
      )
    }
    if(filt_ct == FALSE) { # nolint
      h_out <- ComplexHeatmap::Heatmap(
        h_in,
        col = fun_hm_col,
        name = "Activity Score",
        show_column_names = TRUE,
        show_row_names = TRUE,
        heatmap_width = ggplot2::unit(h_w, "cm"),
        heatmap_height = ggplot2::unit(h_h, "cm"),
        column_title = paste("Top", top_n, "Motifs"),
        column_names_rot = 45,
        column_names_gp = grid::gpar(fontsize = fs_c),
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = fs_r),
        cluster_columns = cl_c,
        cluster_rows = cl_r
      )
    }
  }
  return(h_out)
}