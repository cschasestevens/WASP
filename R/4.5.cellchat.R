#' CellChat Based Cell-cell Communication
#'
#' Conducts CellChat for a selected Seurat object, either provided as
#' a whole dataset or split by a treatment variable for comparison.
#'
#' @param so A Seurat object. Must contain a gene expression assay.
#' @param db_type Database type (either "native" or "custom"). Formats
#' custom lists into CellChat database format. If "custom," appends custom
#' list to the selected native CellChat database.
#' @param db_native Native CellChat database to use.
#' Can provide a custom list of ligand-receptor pairs (see example data).
#' @param db_custom Path to data frame containing custom ligand receptor pairs.
#' @param asy1 Normalized RNA assay to use.
#' @param g_name Group name for comparing cell-cell communication between
#' treatments, provided as a character string. If only one group is present,
#' include an additional variable in the selected Seurat data set indicating
#' that all cells belong to the same group.
#' @param s_name Variable name containing individual sample names.
#' @param ct_col Cell type column name to use in the selected Seurat object.
#' @param cores If not using Windows, how many cores should be used
#' for processing?
#' @param fsize Future globals size (in GB). Only applicable if cores > 1.
#' @return A list containing individual L-R and pathway level communications.
#' @examples
#'
#' # ccmerge <- sc_cc_run(
#' #   so = d,
#' #   db_custom = "ref/20261014_lipid_pairs_for_cellchat_v01.csv",
#' #   g_name = "Group",
#' #   s_name = "Code",
#' #   ct_col = "CellType"
#' # )
#'
#' @export
sc_cc_run <- function(
  so,
  db_type = "native",
  db_native = CellChatDB.human, # nolint
  db_custom = NULL,
  asy1 = "SCT",
  g_name = "Group",
  s_name = "Code",
  ct_col = "CellType",
  cores = 1,
  fsize = 3000
) {
  # Load data
  cat("Setting up environment...", "\n")
  dcc <- so
  Seurat::DefaultAssay(dcc) <- asy1
  # Create/set database
  if (db_type == "native") {
    base_db <- CellChatDB.human # nolint
    CellChatDB <- base_db # nolint
    CellChatDB.use <- CellChatDB # nolint
  }
  if (db_type == "custom") {
    base_db <- CellChatDB.human # nolint
    ccdb <- read.csv(db_custom)
    ccdb <- list(
      interaction = ccdb,
      complex = data.frame(),
      cofactor = data.frame(),
      geneInfo = data.frame(Symbol = unique(c(ccdb$ligand, ccdb$receptor)))
    )
    missing_cols <- setdiff(colnames(base_db$interaction), colnames(ccdb[["interaction"]])) # nolint
    for (col in missing_cols) {
      ccdb[["interaction"]][[col]] <- NA
    }
    ccdb[["interaction"]] <- ccdb[["interaction"]][, colnames(base_db$interaction)] # nolint
    missing_cols <- setdiff(colnames(base_db$geneInfo), colnames(ccdb[["geneInfo"]])) # nolint
    for (col in missing_cols) {
      ccdb[["geneInfo"]][[col]] <- NA
    }
    ccdb[["geneInfo"]] <- ccdb[["geneInfo"]][, colnames(base_db$geneInfo)]
    # Merge databases
    ccdb2 <- CellChat::updateCellChatDB(
      db = ccdb[["interaction"]],
      gene_info = ccdb[["geneInfo"]],
      merged = TRUE,
      species_target = "human"
    )
    CellChatDB <- ccdb2 # nolint
    CellChatDB.use <- CellChatDB # nolint
  }
  ## Subset by condition and run Cell Chat
  if (Sys.info()[["sysname"]] != "Windows" && cores > 1) {
    dcc_l <- setNames(
      lapply(
        seq.int(1, length(unique(dcc@meta.data[[g_name]])), 1),
        function(x) {
          cat("Subsetting data by group column...", "\n")
          dcc <- dcc[
            ,
            dcc[[g_name]] == as.character(unique(dcc@meta.data[[g_name]]))[[x]]
          ]
          dcc@meta.data[[ct_col]] <- factor(
            as.character(dcc@meta.data[[ct_col]]),
            levels = gtools::mixedsort(
              unique(as.character(dcc@meta.data[[ct_col]]))
            )
          )
          dcc_d <- dcc[[asy1]]$data
          dcc_m <- dcc@meta.data
          dcc_m[["samples"]] <- as.factor(dcc_m[[s_name]])
          ct_col1 <- ct_col
          gc(reset = TRUE)
          # Create CellChat object
          print("Creating CellChat object...")
          cc1 <- CellChat::createCellChat(
            object = dcc_d,
            meta = dcc_m,
            group.by = ct_col1
          )
          # Set database
          print("Set CellChat database...")
          cc1@DB <- CellChatDB.use
          # Pre-process to determine over-expressed ligands and receptors
          print("Preprocess CellChat data...")
          cc1 <- CellChat::subsetData(cc1)
          future::plan(
            "multicore",
            workers = cores
          )
          options(future.globals.maxSize = fsize * 1024^2)
          cc1 <- CellChat::identifyOverExpressedGenes(cc1)
          cc1 <- CellChat::identifyOverExpressedInteractions(cc1)
          cc1 <- CellChat::smoothData(cc1, adj = PPI.human) # nolint
          # Compute communication probability
          print("Calculate communication probabilities...")
          cc1 <- CellChat::computeCommunProb(cc1, type = "triMean")
          cc1 <- CellChat::computeCommunProbPathway(cc1)
          # Calculate aggregated cc communication network
          print("Generate communication network...")
          cc1 <- CellChat::aggregateNet(cc1)
          options(future.globals.maxSize = 500 * 1024^2)
          future::plan("sequential")
          print("CellChat completed!")
          return(cc1) # nolint
        }
      ),
      unique(dcc@meta.data[[g_name]])
    )
  } else {
    dcc_l <- setNames(
      lapply(
        seq.int(1, length(unique(dcc@meta.data[[g_name]])), 1),
        function(x) {
          print("Subsetting data by group column...")
          dcc <- dcc[
            ,
            dcc[[g_name]] == as.character(unique(dcc@meta.data[[g_name]]))[[x]]
          ]
          dcc@meta.data[[ct_col]] <- factor(
            as.character(dcc@meta.data[[ct_col]]),
            levels = gtools::mixedsort(
              unique(as.character(dcc@meta.data[[ct_col]]))
            )
          )
          dcc_d <- dcc[[asy1]]$data
          dcc_m <- dcc@meta.data
          dcc_m[["samples"]] <- as.factor(dcc_m[[s_name]])
          ct_col1 <- ct_col
          gc(reset = TRUE)
          # Create CellChat object
          print("Creating CellChat object...")
          cc1 <- CellChat::createCellChat(
            object = dcc_d,
            meta = dcc_m,
            group.by = ct_col1
          )
          # Set database
          print("Set CellChat database...")
          cc1@DB <- CellChatDB.use
          # Pre-process to determine over-expressed ligands and receptors
          print("Preprocess CellChat data...")
          cc1 <- CellChat::subsetData(cc1)
          options(future.globals.maxSize = fsize * 1024^2)
          cc1 <- CellChat::identifyOverExpressedGenes(cc1)
          cc1 <- CellChat::identifyOverExpressedInteractions(cc1)
          cc1 <- CellChat::smoothData(cc1, adj = PPI.human) # nolint
          # Compute communication probability
          print("Calculate communication probabilities...")
          cc1 <- CellChat::computeCommunProb(cc1, type = "triMean")
          cc1 <- CellChat::computeCommunProbPathway(cc1)
          # Calculate aggregated cc communication network
          print("Generate communication network...")
          cc1 <- CellChat::aggregateNet(cc1)
          options(future.globals.maxSize = 500 * 1024^2)
          future::plan("sequential")
          print("CellChat completed!")
          return(cc1) # nolint
        }
      ),
      unique(dcc@meta.data[[g_name]])
    )
  }
  return(dcc_l) # nolint
}

#' CellChat Chord Diagram
#'
#' Creates a chord diagram from a provided CellChat results data frame.
#'
#' @param ccdf A CellChat data frame generat.
#' @param title1 Plot file name, provided as a character string.
#' @param title2 Plot title name.
#' @param cc_type Type of plot to use (either "total" of "comp"). Selecting
#' "total" plots the total number of interactions between cell types, whereas
#' "comp" compares the fold difference between two treatment groups.
#' @param sel_src vector for filtering source communication cell types.
#' @param sel_tar vector for filtering target communication cell types.
#' @param spl_pw Should the data be filtered by signaling pathway? If TRUE,
#' pw_sel must be specified.
#' @param sel_pw Character string of a pathway to use for filtering the
#' provided CellChat object.
#' @param col_g Group name to split data by if cc_type is "comp."
#' @param pw Plot width (in cm).
#' @param ph Plot height (in cm).
#' @param fs1 Font size (between 0.1 and 1).
#' @param fd1 Gap between cell type labels and plot (generally between 1 and 3).
#' @param lgx x coordinate for legend position.
#' @param lgy y coordinate for legend position.
#' @param cc_dir Directory for saving output plot.
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
  pw = 18,
  ph = 18,
  fs1 = 0.6,
  fd1 = 1.8,
  lgx = 0.1,
  lgy = 0.25,
  cc_dir = "analyze/"
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
        cc_dir,
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
      if (min(chrd1[["ratio"]]) > 0.5) {
        col2 <- circlize::colorRamp2(
          c(0.5, max(chrd1[["ratio"]])),
          c("white", "firebrick2"),
          transparency = 0.25
        )
      } else {
        col2 <- circlize::colorRamp2(
          c(min(chrd1[["ratio"]]), 0.5, max(chrd1[["ratio"]])),
          c("dodgerblue3", "white", "firebrick2"),
          transparency = 0.25
        )
      }
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
          cc_dir,
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
      if (min(chrd1[["ratio"]]) > 0.5) {
        col2 <- circlize::colorRamp2(
          c(0.5, max(chrd1[["ratio"]])),
          c("white", "firebrick2"),
          transparency = 0.25
        )
      } else {
        col2 <- circlize::colorRamp2(
          c(min(chrd1[["ratio"]]), 0.5, max(chrd1[["ratio"]])),
          c("dodgerblue3", "white", "firebrick2"),
          transparency = 0.25
        )
      }
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
          cc_dir,
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

#' CellChat Communication Heatmap
#'
#' Generates a CellChat heatmap from a CellChat results data frame,
#' which visualize either ligand-receptor or pathway level results.
#'
#' @param ccdf A CellChat results data frame.
#' @param gname CellChat group for filtering results data frame. The
#' filtered data frame should only contain one group.
#' @param colgrp CellChat data frame group column.
#' @param hmthresh Upper threshold for scaling communcation scores on heatmap.
#' @param hmcol Heatmap gradient color scheme, containing 3 colors.
#' @param transp Transpose the resulting heatmap.
#' @param hmwidth Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param hmheight Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param hmfscol Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param hmfsrow Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param hmclusrow Cluster rows?
#' @param hmrotcol Rotation of column names.
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
sc_cc_heatmap <- function(
  ccdf,
  gname,
  cctype = "PW",
  colgrp = "Group",
  hmthresh = 0.9,
  hmcol = col_grad(scm = 4),
  transp = FALSE,
  hmwidth = 24,
  hmheight = 16,
  hmrotcol = 45,
  hmfscol = 8,
  hmfsrow = 8,
  hmclusrow = TRUE,
  hmlinewidth = 0.1
) {
  # Load CellChat data frame
  h1 <- ccdf
  # Determine CellChat data frame type
  if ("ligand" %in% names(h1)) {
    cat("Detected L-R CellChat results; plotting individual interaction heatmap.", "\n") # nolint
    cctype <- "LR"
  } else {
    cat("Detected pathway CellChat results; plotting aggregated pathway communication heatmap.", "\n") # nolint
    cctype <- "PW"
  }
  if (cctype == "PW") {
    pwlist <- gtools::mixedsort(unique(h1[["pathway_name"]]))
    colpw <- "pathway_name"
  }
  if (cctype == "LR") {
    pwlist <- gtools::mixedsort(unique(h1[["interaction_name"]]))
    colpw <- "interaction_name"
  }
  ctlist <- gtools::mixedsort(unique(c(h1[["source"]], h1[["target"]])))
  # Subset by group
  h2 <- h1[
    h1[[colgrp]] == gname,
    c(colpw, "source", "target", "prob")
  ]
  # Generate matrices for each pathway and combine
  cat("Generating heatmap input matrices", "\n")
  h3 <- setNames(lapply(
    seq.int(1, length(pwlist), 1),
    function(i) {
      tryCatch(
        {
          hin1 <- h2[
            h2[[colpw]] == pwlist[[i]],
            c("source", "target", "prob")
          ]
          hin2 <- reshape2::dcast(
            data = hin1,
            formula = source ~ target,
            value.var = "prob"
          )
          hin3 <- as.matrix(
            magrittr::set_rownames(
              setNames(
                hin2[, 2:ncol(hin2)],
                names(hin2[, 2:ncol(hin2)])
              ),
              hin2[[1]]
            )
          )
          # check for missing rows, append, and sort
          hin4 <- magrittr::set_rownames(
            matrix(
              NA,
              nrow = length(ctlist[ctlist %in% rownames(hin3) == FALSE]),
              ncol = ncol(hin3)
            ),
            ctlist[ctlist %in% rownames(hin3) == FALSE]
          )
          colnames(hin4) <- colnames(hin3)
          hin4 <- rbind(hin4, hin3)
          hin4 <- hin4[ctlist, ]
          # append pathway name to output targets
          colnames(hin4) <- paste(colnames(hin4), pwlist[[i]], sep = ":")
          return(hin4) # nolint
        },
        error = function(e) {
          cat(pwlist[[i]], "matrix generation failed to produce a result.", "\n") # nolint
        }
      )
    }
  ), pwlist)
  ## define heatmap annotation legend based on lengths of individual
  ## heatmap matrices
  ### calculate number of interactions per list element
  cat("Creating heatmap annotation labels", "\n")
  anno1 <- unlist(lapply(
    h3, function(j) if (!is.null(j)) ncol(j) else 0
  ))
  ### create annotation label
  anno1 <- as.factor(unlist(lapply(
    seq.int(1, length(anno1), 1),
    function(k) {
      if (cctype == "PW") {
        plab <- rep(names(anno1)[[k]], anno1[[k]])
      }
      if (cctype == "LR") {
        plab <- rep(
          unique(h1[h1[["interaction_name"]] == names(anno1)[[k]], ][["pathway_name"]]), # nolint
          anno1[[k]]
        )
      }
      return(plab) # nolint
    }
  )))
  cat("Combining inputs", "\n")
  h3 <- h3[lengths(h3) > 0]
  hdata <- do.call(cbind, h3)
  # set color scale
  cat("Defining color scale", "\n")
  qs <- quantile(
    hdata,
    probs = c(
      0.05,
      hmthresh
    ),
    na.rm = TRUE
  )
  fun_hm_col <- circlize::colorRamp2(
    c(
      0,
      round(mean(c(qs[[1]], qs[[2]])), digits = 5),
      qs[[2]]
    ),
    colors = hmcol # nolint
  )
  # add pathway annotation
  cat("Constructing heatmap annotation", "\n")
  if (transp == FALSE) {
    cl2 <- anno1
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
  } else {
    cl2 <- anno1
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
  # output heatmap
  cat("Plotting heatmap", "\n")
  h_out <- ComplexHeatmap::Heatmap(
    matrix = if (transp == TRUE) t(hdata) else hdata,
    bottom_annotation = if (transp == FALSE) hm_anno_col else NULL, # nolint
    left_annotation = if (transp == TRUE) hm_anno_row else NULL, # nolint
    col = fun_hm_col,
    name = "Communication Probability",
    heatmap_width = ggplot2::unit(hmwidth, "cm"),
    heatmap_height = ggplot2::unit(hmheight, "cm"),
    column_names_rot = hmrotcol,
    column_names_gp = grid::gpar(fontsize = hmfscol, fontface = "bold"),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = hmfsrow, fontface = "bold"),
    cluster_columns = FALSE,
    cluster_rows = hmclusrow,
    rect_gp = grid::gpar(col = "black", lwd = hmlinewidth),
    na_col = "white"
  )
  h_out <- list("plot" = h_out, "data" = hdata)
  cat("Successful", "\n")
  return(h_out) # nolint
}
