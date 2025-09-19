#' CellChat Based Cell-cell Communication
#'
#' Conducts CellChat for a selected Seurat object, either provided as
#' a whole data set or split by a treatment variable for comparison.
#'
#' @param so A Seurat object. Must contain a gene expression assay.
#' @param asy1 RNA assay to use.
#' @param g_name Group name for comparing cell-cell communication between
#' treatments, provided as a character string. If only one group is present,
#' include an additional variable in the selected Seurat data set indicating
#' that all cells belong to the same group.
#' @param s_name Variable name containing individual sample names.
#' @param ct_col Cell type column name to use in the selected Seurat object.
#' @return A data frame including either pathway or individual
#' ligand-receptor interactions.
#' @examples
#'
#' # ccmerge <- sc_cc_run(
#' #   so = d,
#' #   asy1 = "RNA",
#' #   g_name = "Group",
#' #   s_name = "Code",
#' #   ct_col = "CellType"
#' # )
#'
#' @export
sc_cc_run <- function(
  so,
  asy1 = "SCT",
  g_name,
  s_name,
  ct_col = "CellType"
) {
  # Load and extract data (subset before creating CellChat if desired)
  dcc <- so
  Seurat::DefaultAssay(dcc) <- asy1
  ## Subset by condition
  dcc_l <- setNames(
    parallel::mclapply(
      mc.cores = 2,
      seq.int(1, length(unique(dcc@meta.data[[g_name]])), 1),
      function(x) {
        dcc <- dcc[
          ,
          dcc[[g_name]] == as.character(unique(dcc@meta.data[[g_name]]))[[x]]
        ]
        dcc@meta.data[["CellType"]] <- factor(
          as.character(dcc@meta.data[["CellType"]]),
          levels = gtools::mixedsort(
            unique(as.character(dcc@meta.data[["CellType"]]))
          )
        )
        dcc_d <- dcc[[asy1]]$data
        dcc_m <- dcc@meta.data
        dcc_m[["samples"]] <- as.factor(dcc_m[[s_name]])
        ct_col1 <- ct_col
        gc(reset = TRUE)
        # Create CellChat object
        cc1 <- CellChat::createCellChat(
          object = dcc_d,
          meta = dcc_m,
          group.by = ct_col1
        )
        # Set database
        CellChatDB <- CellChatDB.human # nolint
        CellChatDB.use <- CellChatDB # nolint
        cc1@DB <- CellChatDB.use
        # Pre-process to determine over-expressed ligands and receptors
        cc1 <- CellChat::subsetData(cc1)
        future::plan("multisession", workers = parallel::detectCores() * 0.25)
        cc1 <- CellChat::identifyOverExpressedGenes(cc1)
        cc1 <- CellChat::identifyOverExpressedInteractions(cc1)
        cc1 <- CellChat::smoothData(cc1, adj = PPI.human) # nolint
        # Compute communication probability
        options(future.globals.maxSize = 3000 * 1024^2)
        cc1 <- CellChat::computeCommunProb(cc1, type = "triMean")
        cc1 <- CellChat::computeCommunProbPathway(cc1)
        # Calculate aggregated cc communication network
        cc1 <- CellChat::aggregateNet(cc1)
        options(future.globals.maxSize = 500 * 1024^2)
        return(cc1)
      }
    ),
    unique(dcc@meta.data[[g_name]])
  )
  saveRDS(dcc_l, "analysis/data.cellchat.list.rds")
  return(dcc_l)
}

#' CellChat Chord Diagram
#'
#' Creates a chord diagram from a provided CellChat results data frame.
#'
#' @param ccdf A CellChat data frame generat.
#' @param title1 Plot title, provided as a character string.
#' @param cc_type Type of plot to use (either "total" of "comp"). Selecting
#' "total" plots the total number of interactions between cell types, whereas
#' "comp" compares the fold difference between two treatment groups
#' @param pw_sel Character string of a pathway to use for filtering the
#' provided CellChat object.
#' @param g_name Group name to split data by if cc_type is "comp."
#' @param pw Plot width (in cm).
#' @param ph Plot height (in cm).
#' @param fs1 Font size (between 0.1 and 1).
#' @param fd1 Gap between cell type labels and plot (generally between 1 and 3).
#' @return A chord diagram displaying either pathway or individual
#' ligand-receptor interactions.
#' @examples
#'
#' # sc_cc_chrd(
#' #   ccdf = ccm[[1]],
#' #   title1 = "plot.cc.LR.sec",
#' #   title2 = "L-R Interactions: Secretory",
#' #   cc_type = "total",
#' #   sel_tar = c("5.Secretory", "16.Secretory"),
#' #   pw = 12,
#' #   ph = 12,
#' #   fs1 = 0.6,
#' #   fd1 = 2.8
#' # )
#'
#' @export
sc_cc_chrd <- function( # nolint
  ccdf,
  title1,
  title2,
  cc_type,
  col_g = NULL,
  sel_pw = NULL,
  sel_src = NULL,
  sel_tar = NULL,
  spl_pw = FALSE,
  pw,
  ph,
  fs1,
  fd1,
  lgx = 0.1,
  lgy = 0.25
) {
  # Input CC data frame
  chrd1 <- ccdf
  # Set global circlize parameters
  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    gap.degree = 5,
    track.margin = c(-0.1, 0.1),
    points.overflow.warning = FALSE
  )
  par(mar = rep(0.5, 4))

  # Plots
  ## Total interactions
  if(cc_type == "total") { # nolint
    ## Subset by pathway, source, or target
    if(!is.null(sel_pw)) { # nolint
      chrd1 <- chrd1[chrd1[["pathway_name"]] == sel_pw, ]
    }
    if(!is.null(sel_src) && is.null(sel_tar)) { # nolint
      chrd1 <- chrd1[chrd1[["source"]] %in% sel_src, ]
    }
    if(!is.null(sel_tar) && is.null(sel_src)) { # nolint
      chrd1 <- chrd1[chrd1[["target"]] %in% sel_tar, ]
    }
    if(!is.null(sel_tar) && !is.null(sel_src)) { # nolint
      chrd1 <- chrd1[
        chrd1[["target"]] %in% sel_tar &
          chrd1[["source"]] %in% sel_src,
      ]
    }
    ## Set colors
    col1 <- setNames(
      col_univ()[1:length(unique( # nolint
        c(
          as.character(chrd1[["target"]]),
          as.character(chrd1[["source"]])
        )
      ))],
      gtools::mixedsort(unique(
        c(
          as.character(chrd1[["target"]]),
          as.character(chrd1[["source"]])
        )
      ))
    )
    # Plot
    png(
      paste(
        "analysis/",
        title1,
        ".png",
        sep = ""
      ),
      width = pw,
      height = ph,
      res = 1200,
      units = "cm"
    )
    circlize::chordDiagram(
      x = chrd1,
      grid.col = col1, # nolint
      transparency = 0.25,
      order = gtools::mixedsort(
        unique(c(
          as.character(chrd1[["target"]]),
          as.character(chrd1[["source"]])
        ))
      ),
      directional = 1,
      direction.type = c("arrows", "diffHeight"),
      diffHeight  = -0.04,
      annotationTrack = "grid",
      annotationTrackHeight = c(0.05, 0.1),
      link.arr.type = "big.arrow",
      link.sort = TRUE,
      link.largest.ontop = TRUE
    )
    circlize::circos.trackPlotRegion(
      track.index = 1,
      bg.border = NA,
      panel.fun = function(x, y) {
        xlim <- circlize::get.cell.meta.data("xlim")
        sector_index <- circlize::get.cell.meta.data("sector.index")
        circlize::circos.text(
          x = mean(xlim),
          y = fd1,
          labels = sector_index,
          facing = "bending",
          cex = fs1
        )
      }
    )
    circlize::circos.clear()
    text(-0, 1.02, title2, cex = 1)
    dev.off()

  }
  ## Compare between treatments
  if(cc_type == "comp") { # nolint
    ## Subset by pathway, source, or target
    if(!is.null(sel_pw)) { # nolint
      chrd1 <- chrd1[chrd1[["pathway_name"]] == sel_pw, ]
    }
    if(!is.null(sel_src) && is.null(sel_tar)) { # nolint
      chrd1 <- chrd1[chrd1[["source"]] %in% sel_src, ]
    }
    if(!is.null(sel_tar) && is.null(sel_src)) { # nolint
      chrd1 <- chrd1[chrd1[["target"]] %in% sel_tar, ]
    }
    if(!is.null(sel_tar) && !is.null(sel_src)) { # nolint
      chrd1 <- chrd1[
        chrd1[["target"]] %in% sel_tar &
          chrd1[["source"]] %in% sel_src,
      ]
    }

    # calculate ratio column
    ## add fold change column
    if(spl_pw == TRUE) { # nolint
      cc_comp2 <- dplyr::count(
        chrd1,
        .data[[col_g]], # nolint
        source, # nolint
        target, # nolint
        pathway_name # nolint
      )
      cc_comp3 <- reshape2::dcast(
        cc_comp2,
        source + target + pathway_name ~ cc_comp2[["Group"]], value.var = "n"
      )
      head(chrd1)
      cc_comp3[is.na(cc_comp3)] <- 0
      cc_comp3[["ratio"]] <- cc_comp3[[4]] / (cc_comp3[[4]] + cc_comp3[[5]])
      chrd1 <- cc_comp3[, c("source", "target", "ratio", "pathway_name")]
      # color schemes and grouping
      nm1 <- sort(unique(chrd1[["pathway_name"]]))
      g1 <- structure(chrd1[["pathway_name"]], names = nm1)
      col2 <- circlize::colorRamp2(
        c(min(chrd1[["ratio"]]), 0.5, max(chrd1[["ratio"]])),
        c("dodgerblue3", "white", "firebrick2"),
        transparency = 0.25
      )
      col1 <- setNames(
        col_univ()[ # nolint
          1:( # nolint
            length(unique(c(
              as.character(chrd1[["target"]]),
              as.character(chrd1[["source"]])
            )))
          )
        ],
        gtools::mixedsort(
          unique(c(
            as.character(chrd1[["target"]]),
            as.character(chrd1[["source"]])
          ))
        )
      )

      # Plot
      png(
        paste(
          "analysis/",
          title1,
          ".png",
          sep = ""
        ),
        width = pw,
        height = ph,
        res = 1200,
        units = "cm"
      )
      circlize::chordDiagram(
        x = chrd1,
        grid.col = col1, # nolint
        col = col2,
        group = g1,
        order = gtools::mixedsort(
          unique(c(
            as.character(chrd1[["target"]]),
            as.character(chrd1[["source"]])
          ))
        ),
        transparency = 0.25,
        directional = 1,
        direction.type = c("arrows", "diffHeight"),
        diffHeight  = -0.04,
        annotationTrack = "grid",
        annotationTrackHeight = c(0.05, 0.1),
        link.arr.type = "big.arrow",
        link.sort = TRUE,
        link.largest.ontop = TRUE
      )
      circlize::circos.trackPlotRegion(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
          xlim <- circlize::get.cell.meta.data("xlim")
          sector_index <- circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(
            x = mean(xlim),
            y = fd1,
            labels = sector_index,
            facing = "bending",
            cex = fs1
          )
        }
      )
      lgd <- ComplexHeatmap::Legend(
        col_fun = col2,
        title = "L-R Interaction Ratio"
      )
      ComplexHeatmap::draw(
        lgd,
        x = grid::unit(1, "npc") - grid::unit(lgx, "mm"), # nolint
        y = grid::unit(lgy, "mm"),
        just = c("right", "bottom")
      )
      circlize::circos.clear()
      text(-0, 1.02, title2, cex = 1)
      dev.off()
    }

    if(spl_pw == FALSE) { # nolint
      cc_comp2 <- dplyr::count(
        chrd1,
        .data[[col_g]], # nolint
        source, # nolint
        target # nolint
      )
      cc_comp3 <- reshape2::dcast(
        cc_comp2,
        source + target ~ cc_comp2[[col_g]], value.var = "n"
      )
      cc_comp3[is.na(cc_comp3)] <- 0
      cc_comp3[["ratio"]] <- cc_comp3[[3]] / (cc_comp3[[3]] + cc_comp3[[4]])
      chrd1 <- cc_comp3[, c("source", "target", "ratio")]
      # color schemes
      col2 <- circlize::colorRamp2(
        c(min(chrd1[["ratio"]]), 0.5, max(chrd1[["ratio"]])),
        c("dodgerblue3", "white", "firebrick2"),
        transparency = 0.25
      )
      col1 <- setNames(
        col_univ()[ # nolint
          1:( # nolint
            length(unique(c(
              as.character(chrd1[["target"]]),
              as.character(chrd1[["source"]])
            )))
          )
        ],
        gtools::mixedsort(
          unique(c(
            as.character(chrd1[["target"]]),
            as.character(chrd1[["source"]])
          ))
        )
      )

      # Plot
      png(
        paste(
          "analysis/",
          title1,
          ".png",
          sep = ""
        ),
        width = pw,
        height = ph,
        res = 1200,
        units = "cm"
      )
      circlize::chordDiagram(
        x = chrd1,
        grid.col = col1, # nolint
        col = col2,
        order = gtools::mixedsort(
          unique(c(
            as.character(chrd1[["target"]]),
            as.character(chrd1[["source"]])
          ))
        ),
        transparency = 0.25,
        directional = 1,
        direction.type = c("arrows", "diffHeight"),
        diffHeight  = -0.04,
        annotationTrack = "grid",
        annotationTrackHeight = c(0.05, 0.1),
        link.arr.type = "big.arrow",
        link.sort = TRUE,
        link.largest.ontop = TRUE
      )
      circlize::circos.trackPlotRegion(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
          xlim <- circlize::get.cell.meta.data("xlim")
          sector_index <- circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(
            x = mean(xlim),
            y = fd1,
            labels = sector_index,
            facing = "bending",
            cex = fs1
          )
        }
      )
      lgd <- ComplexHeatmap::Legend(
        col_fun = col2,
        title = "L-R Interaction Ratio"
      )
      ComplexHeatmap::draw(
        lgd,
        x = grid::unit(1, "npc") - grid::unit(lgx, "mm"), # nolint
        y = grid::unit(lgy, "mm"),
        just = c("right", "bottom")
      )
      circlize::circos.clear()
      text(-0, 1.02, title2, cex = 1)
      dev.off()
    }
  }
}
