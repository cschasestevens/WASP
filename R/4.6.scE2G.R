#' scE2G Pre-processing
#'
#' Creates individual fragment and UMI files for processing
#' with scE2G pipeline, which is run separately using the CLI.
#'
#' @param d_path Character string indicating the path
#' to a set of CellRanger files for data processing.
#' @param study_md A data frame containing sample
#' metadata variables specific to an individual study.
#' @return Data frame containing a list of parameters
#' to use for scRNA-Seq data processing by Seurat.
#' @examples
#'
#' # sc_e2g_input("path/to/data/folder", df1)
#'
#' @export
sc_e2g_input <- function(
  d_path,
  study_md
) {
  list_params <- data.frame(
    # Universal columns
    data.frame(
      # Sample Number
      Sample.No = seq(1:length( # nolint
        basename(
          list.files(d_path)
        )
      )
      ),
      # Individual file names (uses CellRanger folder name by default)
      File.ID = basename(list.files(d_path)),
      # Data file paths (Location of CellRanger files: 'Data/' by default)
      Path = paste(d_path, list.files(d_path), sep = ""),
      # Path to feature files
      Path.feat = paste(
        d_path,
        list.files(d_path),
        "/filtered_feature_bc_matrix/features.tsv.gz",
        sep = ""
      )
    ),
    # Dataset-specific columns
    study_md
  )
  return(list_params) # nolint
}

#' scE2G Locus Plot
#'
#' Creates locus plot for selected gene region from an
#' scE2G results directory.
#'
#' @param so Seurat object containing paths to fragment files.
#' @param asy Assay within Seurat object containing ATAC peaks.
#' @param path_proj Project directory location.
#' @param path_e2g Path to prediction results tsv gz file.
#' @param path_gtf Path to reference genome annotation file.
#' @param gene_name Gene name (character string).
#' @param genome Genome name.
#' @param col_ct_e2g Cell type column in prediction results file.
#' @param col_ct Cell type column in input Seurat object.
#' @param grps Optionally select specific groups for subsetting
#' output data tracks.
#' @param track_window Gene track plotting window, relative to
#' transcription start site for selected gene region.
#' @param track_promoter Promoter site window (in bp).
#' @param thm_disc Discrete color scheme to use.
#' @param thm_cont Continuous color scheme to use.
#' @param font_col Plot font color.
#' @param font_size Plot font size.
#' @param plot_int_anchors Logical indicating whether to plot
#' anchors on the interaction track plot.
#' @param atac_interval bp interval for binning and smoothing
#' ATAC peaks.
#' @param atac_smooth Smoothing factor for ATAC peaks. Higher
#' numbers result in smoother peaks but may lose some information
#' for smaller neighboring peaks.
#' @param font_panel_width Width scaling factor for title panel.
#' @param show_scale_axis Show axis scales for ATAC peak signal
#' and heatmap interaction scores.
#' @return Combined track plot visualizing enhancer-gene expression
#' interactions for a specified gene region.
#' @examples
#'
#' # sc_e2g_locusplot(
#' #   path_proj = "/directory/to/results/",
#' #   path_e2g = "path_proj/path/to/prediction.tsv.gz",
#' #   path_gtf = "path_proj/path/to/gtffile.gtf.gz",
#' #   gene_name = "gene name"
#' # )
#'
#' @import Gviz
#' @export
sc_e2g_locusplot <- function(
  so,
  asy = "ufy.peaks",
  path_proj,
  path_e2g = NULL,
  path_gtf,
  gene_name,
  genome = "hg38",
  col_ct_e2g = "CellType",
  col_ct = "CellType",
  grps = NULL,
  track_window = 60000,
  track_promoter = 2000,
  thm_disc = col_univ(),
  thm_cont = col_grad(scm = 4),
  font_col = "grey25",
  font_size = 0.8,
  plot_int_anchors = FALSE,
  atac_interval = 100,
  atac_smooth = 1,
  font_panel_width = 1,
  show_scale_axis = FALSE
) {
  #---- Set parameters ----
  listp <- list(
    #-- Seurat Object --
    "obj_seurat" = so,
    "asy" = asy,
    #-- Directories --
    "path_proj" = path_proj,
    ## E2G results can be a vector of paths (for multiple comparisons)
    "path_e2g" = path_e2g,
    "path_gtf" = path_gtf, # nolint
    #-- Global params --
    "genome" = genome,
    "gene_name" = gene_name,
    "col_ct_e2g" = col_ct_e2g,
    "col_ct" = col_ct,
    "grps" = grps,
    "window_track" = track_window,
    "window_promoter" = track_promoter,
    "theme_color_discrete" = thm_disc,
    "theme_color_continuous" = thm_cont,
    "font_color" = font_col,
    "font_size" = font_size,
    #-- Interaction plot params --
    "plot_int_anchors" = plot_int_anchors,
    "plot_int_type" = "arcs",
    #-- ATAC plot params --
    "atac_interval" = atac_interval,
    "atac_smoothing" = atac_smooth,
    #-- Other plotting parameters --
    "font_panel_width" = font_panel_width,
    "show_scale_axis" = show_scale_axis
  )

  #---- 1. Gene Axis Track ----
  # Axis track doesn't change with increasing numbers of samples
  fun.track.axis <- function(lparam) { # nolint
    Gviz::GenomeAxisTrack(
      col = lparam[["font_color"]],
      fontcolor = lparam[["font_color"]],
      cex = lparam[["font_size"]],
      lwd = 1,
      labelPos = "below",
      distFromAxis = 25,
      background.title = "white",
      col.title = lparam[["font_color"]]
    )
  }
  #---- 2. Gene Region Track ----
  fun.track.gene <- function(lparam, fillp) { # nolint
    if (class(lparam[["path_gtf"]]) == "character") {
      cat("Importing gtf file...", "\n")
      gtf1 <- rtracklayer::import(paste0(lparam[["path_proj"]], lparam[["path_gtf"]])) # nolint
    } else {
      cat("Loading existing gtf GRanges object...", "\n")
      gtf1 <- lparam[["path_gtf"]]
    }
    # Select exon regions for selected gene
    gtf1sel <- gtf1[
      gtf1$type == "exon" &
        gtf1$gene_name == lparam[["gene_name"]]
    ]
    gtf1sel <- range(gtf1sel) + lparam[["window_track"]]
    ## Filter by ranges including specified gene
    gtf1ex <- IRanges::subsetByOverlaps(gtf1, gtf1sel)
    ## Filter by exon and protein coding regions
    gtf1ex <- gtf1ex[
      gtf1ex$type == "exon" &
        gtf1ex$gene_type == "protein_coding"
    ]
    # format GRanges and select primary entry for each gene in range
    gtf_select <- unlist(lapply(
      seq.int(1, length(unique(gtf1ex$gene_name)), 1),
      function(i) {
        gtf_name <- gtf1ex[
          gtf1ex$gene_name == unique(gtf1ex$gene_name)[[i]]
        ]
        gtf_ranges <- range(split(gtf_name, gtf_name$transcript_id))
        gtf_lengths <- unlist(width(gtf_ranges))
        gtf_select <- names(which.max(gtf_lengths))
        return(gtf_select) # nolint
      }
    ))
    gtf1ex <- gtf1ex[
      grepl(paste(gtf_select, collapse = "|"), gtf1ex$transcript_id)
    ]
    # set plotting window for track (centered around TSS for specific gene)
    options(ucscChromosomeNames = FALSE)
    chr1 <- as.character(seqnames(gtf1ex)[1])
    gene_start <- min(start(gtf1ex[gtf1ex$gene_name == lparam[["gene_name"]]]))
    gene_end   <- max(end(gtf1ex[gtf1ex$gene_name == lparam[["gene_name"]]]))
    gene_strand <- as.character(strand(gtf1ex[gtf1ex$gene_name == lparam[["gene_name"]]])[1]) # nolint
    gene_tss <- if (gene_strand == "+") {
      gene_start
    } else {
      gene_end
    }
    # create track
    track_gene <- Gviz::GeneRegionTrack(
      gtf1ex,
      genome = lparam[["genome"]],
      chromosome = chr1,
      transcript = gtf1ex$transcript_id,
      gene = gtf1ex$gene_id,
      feature = gtf1ex$type,
      transcriptAnnotation = "symbol",
      collapseTranscripts = FALSE,
      symbol = gtf1ex$gene_name,
      stacking = "squish",
      fill = fillp,
      col = lparam[["font_color"]],
      col.line = lparam[["font_color"]],
      cex.group = lparam[["font_size"]],
      fontface.group = 3,
      background.title = "white",
      col.title = lparam[["font_color"]],
      name = chr1,
      rotation.title = if (length(lparam[["path_e2g"]]) > 1) 0 else 90, # nolint
      just.title = "center",
      vjust.title = 0.5,
      cex.title = lparam[["font_size"]]
    )
    return(list( # nolint
      "ranges" = gtf1ex,
      "track_gene" = track_gene,
      "chromosome" = chr1,
      "strand" = gene_strand,
      "TSS" = gene_tss
    ))
  }

  #---- 3. Enhancer Interaction Arcs ----
  fun.track.interact <- function(lparam) { # nolint
    list_track <- setNames(
      lapply(
        seq.int(1, length(lparam[["path_e2g"]]), 1),
        function(i) {
          # Read E2G results file
          intp <- read.delim(
            gzfile(paste0(lparam[["path_proj"]], lparam[["path_e2g"]][[i]])),
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          # Filter by selected gene
          intp2 <- intp[
            grepl(paste(unique(track_gene[["ranges"]]$gene_name), collapse = "|"), intp[["TargetGene"]]), # nolint
          ]
          # Construct GenomicInteractions object (two anchor sets and metadata)
          ## Construct GRanges from filtered results
          genint_enhancer <- GenomicRanges::GRanges(
            seqnames = unique(intp2[["chr"]]),
            ranges = IRanges::IRanges(
              start = intp2[["start"]],
              end = intp2[["end"]]
            )
          )
          ## Construct TSS GRanges
          genint_tss <- rep(
            GenomicRanges::GRanges(
              seqnames = track_gene[["chromosome"]],
              ranges = IRanges::IRanges(
                start = max(1, track_gene[["TSS"]] - lparam[["window_promoter"]]), # nolint
                end = track_gene[["TSS"]] + lparam[["window_promoter"]]
              ),
              strand = track_gene[["strand"]]
            ),
            if (length(genint_enhancer) == 0) 1 else length(genint_enhancer)
          )
          ## Construct metadata DataFrame
          ### NOTE: DataFrame is different from data.frame
          genint_meta <- S4Vectors::DataFrame(
            scE2G_score = if (length(genint_enhancer) == 0) 0 else intp2[["E2G.Score.qnorm"]], # nolint
            ATAC.counts = if (length(genint_enhancer) == 0) 0 else intp2[["normalizedATAC_enh"]] # nolint
          )
          ## Combine elements to create GenomicInteractions object
          genint <- GenomicInteractions::GenomicInteractions(
            anchor1 = if (length(genint_enhancer) == 0) genint_tss else genint_enhancer, # nolint
            anchor2 = genint_tss
          )
          S4Vectors::mcols(genint) <- cbind(
            S4Vectors::mcols(genint), genint_meta
          )
          ## Create gradient color scale
          interaction_scores <- S4Vectors::mcols(genint)[["scE2G_score"]]
          interaction_scores <- pmax(0, pmin(1, interaction_scores))
          color_bins <- cut(
            interaction_scores, breaks = seq(0, 1, length.out = 101),
            labels = FALSE, include.lowest = TRUE
          )
          custom_colors <- sapply(interaction_scores, function(s) {
            alpha_floor <- 0.1 + (0.85 * s)
            rgb(red = 0.2 * s,
                green = 0.4 * s,
                blue = 0.9 * s,
                alpha = alpha_floor)
          })
          S4Vectors::mcols(genint)$value <- paste0("bin_", color_bins)
          named_vector <- custom_colors
          names(named_vector) <- S4Vectors::mcols(genint)$value
          named_vector <- named_vector[!duplicated(names(named_vector))]
          ## Plot interaction track
          track_int <- GenomicInteractions::InteractionTrack(
            genint,
            chromosome = track_gene[["chromosome"]],
            name = if (length(lparam[["path_e2g"]]) > 1) paste("") else "Enhancer-Gene\nInteractions" # nolint
          )
          ## Adjust plotting parameters
          Gviz::displayPars(track_int) <- list(
            # Interaction Track
            col.interactions = named_vector,
            plot.anchors = lparam[["plot_int_anchors"]],
            interaction.type = lparam[["plot_int_type"]],
            background.title = "white",
            col.title = lparam[["font_color"]],
            cex.title = lparam[["font_size"]]
          )
          return(list("track_int" = track_int, "ranges" = if (length(genint_enhancer) == 0) genint_tss else genint_enhancer, "metadata" = genint_meta)) # nolint
        }
      ),
      unlist(lapply(
        seq.int(1, length(lparam[["path_e2g"]]), 1),
        function(i) {
          unique(
            read.delim(
              gzfile(paste0(lparam[["path_proj"]], lparam[["path_e2g"]][[i]])),
              header = TRUE,
              sep = "\t",
              stringsAsFactors = FALSE,
              check.names = FALSE
            )[[lparam[["col_ct_e2g"]]]]
          )
        }
      ))
    )
    return(list_track) # nolint
  }
  #---- 4. Plot Normalized ATAC Signal ----
  ## For all peaks included within the track window
  fun.track.atac <- function(lparam, fillp) { # nolint
    # Calculate scaling factors for normalizing ATAC signal per
    # cell groups included in plot
    cell_depths <- lparam[["obj_seurat"]][[
      paste("nCount", lparam[["asy"]], sep = "_")
    ]]
    groups <- as.character(unique(
      lparam[["obj_seurat"]]@meta.data[[lparam[["col_ct"]]]]
    ))
    # match order of group names if also plotting scE2G interaction data
    if (!is.null(lparam[["path_e2g"]])) {
      groups <- groups[match(names(track_int), groups)]
    }
    scale_factors <- setNames(
      unlist(lapply(
        seq.int(1, length(groups), 1),
        function(g) {
          cells_in_group <- rownames(lparam[["obj_seurat"]]@meta.data)[
            lparam[["obj_seurat"]]@meta.data[[lparam[["col_ct"]]]] == groups[[g]] # nolint
          ]
          n_cells <- length(cells_in_group)
          mean_depth <- mean(cell_depths[cells_in_group, ])
          return(n_cells * mean_depth) # nolint
        }
      )),
      groups
    )
    median_scale <- median(scale_factors)
    ## Define plotting window
    atac_window <- GenomicRanges::GRanges(
      seqnames = track_gene[["chromosome"]],
      ranges = IRanges::IRanges(
        start = max(1, track_gene[["TSS"]] - lparam[["window_track"]]),
        end = track_gene[["TSS"]] + lparam[["window_track"]]
      )
    )
    ## Extract fragments for each sample
    # Extract fragments from fragment file
    frag1 <- dplyr::bind_rows(lapply(
      seq.int(
        1,
        length(Signac::Fragments(
          lparam[["obj_seurat"]][[lparam[["asy"]]]]
        )),
        1
      ),
      function(k) {
        # Return fragments and cells
        frag1 <- Signac::Fragments(
          lparam[["obj_seurat"]][[lparam[["asy"]]]]
        )[[k]]
        frag_path <- Signac::GetFragmentData(frag1, slot = "path")
        cell_map <- Signac::GetFragmentData(frag1, slot = "cells")
        # Query tabix-indexed fragment file
        tabix_res <- Rsamtools::scanTabix(
          frag_path, param = atac_window
        )[[1]]
        frag_df <- read.table(
          text = tabix_res, sep = "\t",
          col.names = c("chr", "start", "end", "barcode", "count")
        )
        # append donor barcode prior to combining list
        # harmonize cell names by appending donor code
        frag_df[["barcode"]] <- paste(unique(gsub(
          "_[^_]*$", "", names(cell_map)
        )), frag_df[["barcode"]], sep = "_")
        # Filter for valid cells in your Seurat object and map barcodes
        frag_df <- frag_df[frag_df$barcode %in% names(cell_map), ]
        frag_df$cell <- cell_map[frag_df$barcode]
        frag_df$group <- lparam[["obj_seurat"]]@meta.data[
          frag_df[["barcode"]],
        ][[lparam[["col_ct"]]]]
        return(frag_df) # nolint
      }
    ))
    # Convert fragment ends into individual Tn5 insertion events
    # Each fragment contributes two 1-bp insertion sites: start and end
    insertions <- rbind(
      data.frame(pos = frag1$start, group = frag1$group),
      data.frame(pos = frag1$end, group = frag1$group)
    )
    # Restrict to exact region boundaries
    reg_start <- start(atac_window)
    reg_end <- end(atac_window)
    insertions <- insertions[
      insertions$pos >= reg_start & insertions$pos <= reg_end,
    ]
    region_length <- reg_end - reg_start + 1
    ## Create tracks for all groups present in Seurat object
    list_atac <- setNames(
      lapply(
        seq.int(1, length(groups), 1),
        function(j) {
          # Subset insertions for selected group
          grp_pos <- insertions$pos[insertions$group == groups[[j]]]
          # Count insertions at every base pair (1-indexed relative to_start)
          bp_counts <- tabulate(grp_pos - reg_start + 1, nbins = region_length)
          # Normalize by group scaling factor and rescale to dataset median
          norm_signal <- (bp_counts / scale_factors[j]) * median_scale
          # Create bins for track window
          atac_window2 <- unlist(IRanges::tile(atac_window, width = lparam[["atac_interval"]])) # nolint
          # perform cubic spline interpolation to smooth peaks
          coord_x <- start(atac_window2) + (width(atac_window2) / 2)
          # bin signal per atac window
          bins <- ceiling(seq_along(norm_signal) / lparam[["atac_interval"]])
          # Calculate sum per bin
          binned_signal <- as.numeric(tapply(
            norm_signal, bins, sum, na.rm = TRUE
          ))
          S4Vectors::mcols(atac_window2)[["counts.atac"]] <- binned_signal
          coord_y <- atac_window2$counts.atac
          # Smoothing factor (higher number increases smoothing)
          coord_int <- length(coord_x) / lparam[["atac_smoothing"]]
          coord_spline <- spline(
            x = coord_x,
            y = coord_y,
            n = coord_int
          )
          # Replace negative values
          coord_spline[["y"]][coord_spline[["y"]] < 0] <- 0
          # Reconstruct smoothed granges
          atac_grange <- GenomicRanges::GRanges(
            seqnames = track_gene[["chromosome"]],
            ranges = IRanges::IRanges(
              start = round(coord_spline[["x"]]),
              end = round(coord_spline[["x"]])
            ),
            counts.atac = coord_spline[["y"]]
          )
          # Plot track
          track_atac <- Gviz::DataTrack(
            range = atac_grange,
            genome = lparam[["genome"]],
            chromosome = track_gene[["chromosome"]],
            type = "polygon",
            data = "counts.atac",
            background.title = "white",
            col.title = lparam[["font_color"]],
            cex.title = lparam[["font_size"]],
            name = if (length(groups) > 1) paste0(groups[[j]], "\n", "Peaks") else "Normalized\nATAC Signal", # nolint
            col.axis = lparam[["font_color"]],
            cex.axis = lparam[["font_size"]],
            fill.mountain = rep(fillp[[j]], 2),
            col.mountain = lparam[["font_color"]]
          )
          # Highlight enhancer peaks for separate track if significant
          if (!is.null(lparam[["path_e2g"]])) { # nolint
            if (length(track_int) > 1) {
              Gviz::displayPars(track_atac) <- list(
                rotation.title = 0,
                just.title = "center",
                vjust.title = 0
              )
            }
            if (max(track_int[[groups[[j]]]][["metadata"]][["scE2G_score"]]) > 0) { # nolint
              atac_grange_enh <- track_int[[groups[[j]]]][["ranges"]]
              track_atac <- Gviz::HighlightTrack(
                trackList = list(track_atac),
                start = start(atac_grange_enh),
                end = end(atac_grange_enh),
                chromosome = track_gene[["chromosome"]],
                inBackground = FALSE,
                fill = adjustcolor(c("grey25"), alpha.f = 0.3),
                col = "grey10"
              )
            }
          }
          return(track_atac) # nolint
        }
      ),
      groups
    )
    return(list_atac) # nolint
  }

  #---- 5. Generate E2G Scores Heatmap ----
  fun.track.heat <- function(lparam) { # nolint
    list_heat <- setNames(
      lapply(
        seq.int(1, length(track_int), 1),
        function(k) {
          heat_grange <- track_int[[k]][["ranges"]]
          S4Vectors::mcols(heat_grange)[["score"]] <- track_int[[k]][["metadata"]][["scE2G_score"]] # nolint
          ## Plot data track
          track_heat <- Gviz::DataTrack(
            range = heat_grange,
            genome = lparam[["genome"]],
            chromosome = track_gene[["chromosome"]],
            type = "heatmap",
            data = "score",
            gradient = if (length(S4Vectors::mcols(heat_grange)[["score"]]) == 1) c("grey75") else lparam[["theme_color_continuous"]], # nolint
            background.title = "white",
            col.title = lparam[["font_color"]],
            cex.title = lparam[["font_size"]],
            name = if (length(track_int) > 1) paste0(names(track_int)[[k]]) else "Score", # nolint
            col.axis = lparam[["font_color"]],
            cex.axis = lparam[["font_size"]]
          )
          if (length(track_int) > 1) {
            Gviz::displayPars(track_heat) <- list(
              rotation.title = 0,
              just.title = "center",
              vjust.title = 0
            )
          }
          return(track_heat) # nolint
        }
      ),
      names(track_int)
    )
    return(list_heat) # nolint
  }

  #---- 6. Plot Output ----
  # Create axis
  cat("Creating genome track axis\n")
  track_axis <- fun.track.axis(listp)
  # Create gene region
  cat("Plotting genome region on axis\n")
  track_gene <- fun.track.gene(listp, listp[["theme_color_discrete"]][[2]])
  if (!is.null(listp[["path_e2g"]])) {
    cat("Enhancer region predictions present; overlaying predicted enhancer regions on plot\n") # nolint
    # Create Interaction arcs
    track_int <- fun.track.interact(listp)
    # Create heatmaps
    track_heat <- fun.track.heat(listp)
    if (!is.null(listp[["grps"]])) {
      track_int <- track_int[listp[["grps"]]]
    }
    if (!is.null(listp[["grps"]])) {
      track_heat <- track_heat[listp[["grps"]]]
    }
  }
  # Create ATAC peak plots
  cat("Plotting ATAC peaks within selected region\n")
  track_atac <- fun.track.atac(listp, listp[["theme_color_discrete"]])
  # filter to include specific groups
  if (!is.null(listp[["grps"]])) {
    track_atac <- track_atac[listp[["grps"]]]
  }
  cat("Generating final plot with all tracks\n")
  ## Final output plot
  if (length(track_atac) == 1) {
    plot_output <- Gviz::plotTracks(
      list(
        # Interaction arcs
        if (!is.null(listp[["path_e2g"]])) track_int[[1]][["track_int"]] else NULL, # nolint
        # ATAC signal
        track_atac[[1]],
        # scE2G heatmap
        if (!is.null(listp[["path_e2g"]])) track_heat[[1]] else NULL,
        # Gene region with exons displayed
        track_gene[["track_gene"]],
        # Main genome axis track
        track_axis
      ),
      chromosome = track_gene[["chromosome"]],
      from = max(1, track_gene[["TSS"]] - listp[["window_track"]]),
      to = track_gene[["TSS"]] + listp[["window_track"]],
      sizes = c(
        if (!is.null(listp[["path_e2g"]])) 2 else NULL,
        3,
        if (!is.null(listp[["path_e2g"]])) 1 else NULL,
        1,
        2
      ),
      panel.margin = 25,
      title.width = listp[["font_panel_width"]],
      showAxis = listp[["show_scale_axis"]]
    )
  } else {
    ## Format output list and provide
    ## to plotTracks function
    track_list <- do.call(
      c,
      list(
        # Interaction/ATAC pairs
        unlist(lapply(
          seq.int(1, length(track_atac), 1),
          function(x) {
            if (!is.null(listp[["path_e2g"]])) {
              c(
                track_int[[x]][["track_int"]],
                track_atac[[x]]
              )
            } else {
              track_atac[[x]]
            }
          }
        )),
        # All e2g heatmaps
        if (!is.null(listp[["path_e2g"]])) track_heat else NULL,
        # Gene track and axis
        track_gene[["track_gene"]],
        track_axis
      )
    )
    plot_output <- Gviz::plotTracks(
      track_list,
      chromosome = track_gene[["chromosome"]],
      from = max(1, track_gene[["TSS"]] - listp[["window_track"]]),
      to = track_gene[["TSS"]] + listp[["window_track"]],
      sizes = c(
        # Interaction/ATAC pairs
        if (!is.null(listp[["path_e2g"]])) rep(c(2, 3), length(c(track_atac))) else rep(3, length(c(track_atac))), # nolint
        # E2G heatmaps
        if (!is.null(listp[["path_e2g"]])) rep(2, length(track_heat)) else NULL, # nolint
        # Gene Region
        1.5,
        # Gene Axis
        2
      ),
      panel.margin = 25,
      title.width = listp[["font_panel_width"]],
      showAxis = listp[["show_scale_axis"]]
    )
  }
  cat("Complete\n")
  return(plot_output) # nolint
}
