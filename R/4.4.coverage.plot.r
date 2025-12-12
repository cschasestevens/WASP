#' Add JASPAR TF Motifs
#'
#' Creates a ChromatinAssay from a scATAC-Seq counts assay
#' and annotates transcription factor motifs based on
#' the JASPAR database.
#'
#' @param so An object of class Seurat. Must contain an ATAC assay.
#' @param asy1 Assay to use for annotating transcription factors.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @param dir1 Output directory for exporting motif names.
#' @return An annotated ChromatinAssay object.
#' @examples
#'
#' # d_chrmasy <- sc_atac_motifs(
#' #   d,
#' #   "ufy.peaks",
#' #   rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' #    )
#' # )
#'
#' @export
sc_atac_motifs <- function(
  so,
  asy1 = "ufy.peaks",
  dref,
  dir1 = "processing/"
) {
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Hsapiens.UCSC.hg38)
  ## Extract counts from Seurat ATAC assay
  d <- so
  Seurat::DefaultAssay(d) <- asy1
  dc <- SeuratObject::GetAssayData(d, slot = "counts")
  ## Format gene locations and create ChromatinAssay
  ref1 <- dref
  gene_coords <- ref1[ref1$type == "gene"]
  gene_coords$gene_biotype <- gene_coords$gene_type
  GenomeInfoDb::seqlevelsStyle(gene_coords) <- "UCSC"
  gene_coords <- GenomeInfoDb::keepStandardChromosomes(
    gene_coords,
    pruning.mode = "coarse"
  )
  d1 <- Signac::CreateChromatinAssay(
    counts = dc,
    sep = c("-", "-"),
    verbose = TRUE,
    meta.data = d@meta.data
  )
  ## Use CORE, CNE, PBM, PBM_HLH, PHYLOFACTS, POLII, SPLICE, and PBM_HOMEO
  list_clt1 <- c(
    "CORE", "CNE", "PBM", "PBM_HLH",
    "PBM_HOMEO", "PHYLOFACTS", "POLII",
    "SPLICE"
  )
  list_clt <- setNames(
    lapply(
      list_clt1,
      function(x) {
        opt1 <- list(
          collection = x,
          tax_group = "vertebrates",
          all_versions = FALSE
        )
        return(opt1) # nolint
      }
    ),
    list_clt1
  )
  ## Query JASPAR and combine position frequency matrices
  list_pfm <- setNames(
    lapply(
      list_clt,
      function(y) {
        pfm1 <- TFBSTools::getMatrixSet(
          x = JASPAR2020, # nolint
          opts = y
        )
        return(pfm1) # nolint
      }
    ),
    list_clt1
  )
  list_pfm <- list_pfm[lengths(list_pfm) > 0]
  ## 4 JASPAR collections have PFMs
  list_pfm <- c(
    list_pfm[[1]],
    list_pfm[[2]],
    list_pfm[[3]],
    list_pfm[[4]]
  )
  ## Add motifs to chromatinassay
  d1 <- Signac::AddMotifs(
    object = d1,
    genome = BSgenome.Hsapiens.UCSC.hg38, # nolint
    pfm = list_pfm
  )
  ## Save motif list
  tf_list <- Signac::GetMotifData(d1)
  tf_list <- data.frame(
    "ID" = colnames(tf_list),
    "name" = unlist(
      lapply(
        colnames(tf_list),
        function(x) {
          name(TFBSTools::getMatrixByID(JASPAR2020, ID = x)) # nolint
        }
      )
    )
  )
  write.table(
    tf_list,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    file = paste(
      dir1, "motifs_jaspar2020.txt", sep = ""
    )
  )
  return(d1)
}

#' scATAC-Seq Coverage Plot
#'
#' Generates a coverage plot from a Signac ChromatinAssay.
#' Requires scATAC-Seq peak information and a reference GRanges object
#' for plotting gene positions. The Signac function CoveragePlot() is
#' used to generate the final plot.
#'
#' @param so Input Seurat Object.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @param asy1 ATAC peaks assay to use.
#' @param asy2 GEX assay to use.
#' @param g_name Gene or region to plot, provided as a character string.
#' @param g_region (optional) Region within the coverage plot to highlight.
#' @param bp_window Vector containing the upstream and downstream limits
#' for extending the plotting window (in base pairs).
#' @param ct_col Cell type column name to use in the selected Seurat object.
#' Must be a factor variable.
#' @param ct_spl Logical indicating whether cell types should be stratified
#' by an additional variable.
#' @param ct_sel Pattern for selecting individual cell types if ct_spl == TRUE.
#' @param md_list1 Name of variable to split by group.
#' @param p_type Final coverage plot type
#' (either "ct_only", "gn_only", or "gn_ct");
#' Note that "gn_only" and "gn_ct" use the variable
#' from md_list1 for plotting. This can be useful for summarizing differences
#' in accessbility between treatment groups.
#' @return A coverage plot including the specified gene track and all other
#' genes present within the specified window.
#' @examples
#'
#' # p_cov <- sc_coverage_plot2(
#' #   so = d,
#' #   dref = rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' #   ),
#' #   asy1 = "ufy.peaks",
#' #   g_name = "SCGB3A2",
#' #   bp_window = 5000,
#' #   ct_col = "CellType",
#' #   ct_spl = TRUE,
#' #   ct_sel = "Secretory",
#' #   md_list1 = c("Group"),
#' #   p_type = "gn_ct"
#' # )
#'
#' @export
sc_coverage_plot <- function( # nolint
  so,
  dref,
  asy1 = "ufy.peaks",
  asy2 = "RNA",
  asy_dat = "scale.data",
  g_name,
  g_region = NULL,
  bp_window = c(5000, 5000),
  ct_col = "CellType",
  ct_spl = FALSE,
  ct_filt = FALSE,
  ct_sel = NULL,
  md_list1 = NULL,
  p_type = "ct_only"
) {
  # Set parameters
  d <- so
  a1 <- asy1
  a2 <- asy2
  gn <- g_name
  adat <- asy_dat
  rg <- dref
  bpw <- bp_window
  ctc <- ct_col
  cts <- ct_spl
  ctf <- ct_filt
  ctl <- ct_sel
  mdl <- md_list1
  plt <- p_type
  ## Highlight region (optional)
  if(!is.null(g_region)) { # nolint
    g_reg <- Signac::StringToGRanges(g_region)
    g_reg$color <- "orange"
  }
  if(is.null(g_region)) { # nolint
    g_reg <- NULL
  }
  # Format data for plot input
  ## Change assay
  Seurat::DefaultAssay(d) <- a1
  ## Set gene annotation
  Signac::Annotation(d[[a1]]) <- rg
  ## Calculate linkages for selected genes
  d <- Signac::RegionStats(
    d,
    genome = BSgenome.Hsapiens.UCSC.hg38 # nolint
  )
  dlnk <- tryCatch(
    {
      d <- Signac::LinkPeaks(
        object = d,
        peak.assay = a1,
        expression.assay = a2,
        genes.use = c(
          gn
        )
      )
    },
    error = function(e) {
      print("No linked peaks; insufficient peak counts for selected gene...") # nolint
    }
  )
  ### Set links to FALSE if unsuccessful
  if (class(dlnk) != "character") {
    d <- dlnk
    lnk <- TRUE
  }
  if (class(dlnk) == "character") {
    lnk <- FALSE
  }
  remove(dlnk)
  ## Filter data by cell type (optional)
  if (ctf == TRUE) {
    d2 <- d[, grepl(ctl, d@meta.data[[ctc]])]
  }
  if (ctf == FALSE) {
    d2 <- d
  }
  # If creating coverage plots split by group
  if (cts == TRUE) {
    # Add split column to existing Seurat object
    d_spl <- data.frame(
      "CellType.split" = gtools::mixedsort(
        paste(
          unique(d2@meta.data[, c(mdl, ctc)])[[ctc]],
          unique(d2@meta.data[, c(mdl, ctc)])[[mdl]],
          sep = "."
        )
      )
    )
    d_spl[[ctc]] <- factor(
      gsub(
        paste(
          ".", unique(d2@meta.data[[mdl]])[[1]],
          "|.", unique(d2@meta.data[[mdl]])[[2]],
          sep = ""
        ),
        "",
        d_spl[["CellType.split"]]
      ),
      levels = levels(d2@meta.data[[ctc]])
    )
    d_spl[[mdl]] <- gsub(
      ".*\\.",
      "",
      d_spl[["CellType.split"]]
    )
    d_spl[["CellType.split"]] <- factor(
      d_spl[["CellType.split"]],
      levels = gtools::mixedsort(d_spl[["CellType.split"]])
    )
    d_ct2 <- dplyr::select(
      dplyr::left_join(
        d2@meta.data[, c(ctc, mdl)],
        d_spl,
        by = c(ctc, mdl)
      ),
      .data[["CellType.split"]] # nolint
    )
    d2 <- SeuratObject::AddMetaData(
      d2,
      d_ct2[["CellType.split"]],
      col.name = "CellType.split"
    )
    ctc2 <- "CellType.split"
    Seurat::Idents(d2) <- ctc2
  }
  # If creating standard coverage plots
  if (cts == FALSE) {
    ctc2 <- ctc
    Seurat::Idents(d2) <- ctc2
  }
  # Create final plots
  ## Combine group variable and cell type panels
  if (plt == "gn_ct") {
    p1 <- Signac::CoveragePlot(
      object = d2,
      region = gn,
      features = gn,
      expression.assay = a2,
      expression.slot = adat,
      annotation = TRUE,
      peaks = TRUE,
      links = lnk,
      idents = unique(d2@meta.data[[ctc2]]),
      extend.upstream = bpw[[1]],
      extend.downstream = bpw[[2]],
      region.highlight = g_reg
    ) &
      ggplot2::scale_fill_manual(values = col_univ()) & # nolint
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 8),
        axis.title.y = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 6),
        plot.margin = ggplot2::unit(
          c(0.25, 0.05, 0.25, 0.25),
          "cm"
        )
      )
    Seurat::Idents(d2) <- mdl
    p2 <- Signac::CoveragePlot(
      object = d2,
      region = gn,
      features = gn,
      expression.assay = a2,
      expression.slot = adat,
      annotation = TRUE,
      peaks = TRUE,
      links = lnk,
      idents = sort(unique(d2@meta.data[[mdl]])),
      extend.upstream = bpw[[1]],
      extend.downstream = bpw[[2]],
      region.highlight = g_reg
    ) &
      ggplot2::scale_fill_manual(values = col_univ()) & # nolint
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 8),
        axis.title.y = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(size = 6),
        axis.text.x = ggplot2::element_text(size = 6),
        plot.margin = ggplot2::unit(
          c(0.25, 0, 0.25, 0),
          "cm"
        )
      )
    p_comb <- p2 +
      p1 +
      patchwork::plot_layout(ncol = 1, heights = c(0.75, 4)) + # nolint
      ggplot2::theme(
        plot.margin = ggplot2::unit(
          c(0.25, 0.0, 0.25, 0),
          "cm"
        )
      )
  }
  if (plt == "gn_only") {
    Seurat::Idents(d2) <- mdl
    p2 <- Signac::CoveragePlot(
      object = d2,
      region = gn,
      features = gn,
      expression.assay = a2,
      expression.slot = adat,
      annotation = TRUE,
      peaks = TRUE,
      links = lnk,
      idents = sort(unique(d2@meta.data[[mdl]])),
      extend.upstream = bpw[[1]],
      extend.downstream = bpw[[2]],
      region.highlight = g_reg
    ) &
      ggplot2::scale_fill_manual(values = col_univ()) & # nolint
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 8),
        axis.title.y = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(size = 6),
        axis.text.x = ggplot2::element_text(size = 6),
        plot.margin = ggplot2::unit(
          c(0.25, 0, 0.25, 0),
          "cm"
        )
      )
    p_comb <- p2
  }
  if (plt == "ct_only") {
    p1 <- Signac::CoveragePlot(
      object = d2,
      region = gn,
      features = gn,
      expression.assay = a2,
      expression.slot = adat,
      annotation = TRUE,
      peaks = TRUE,
      links = lnk,
      idents = sort(unique(d2@meta.data[[ctc2]])),
      extend.upstream = bpw[[1]],
      extend.downstream = bpw[[2]],
      region.highlight = g_reg
    ) &
      ggplot2::scale_fill_manual(values = col_univ()) & # nolint
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 8),
        axis.title.y = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(size = 6),
        axis.text.x = ggplot2::element_text(size = 6),
        plot.margin = ggplot2::unit(
          c(0.25, 0, 0.25, 0),
          "cm"
        )
      )
    p_comb <- p1
  }
  return(p_comb)
}
