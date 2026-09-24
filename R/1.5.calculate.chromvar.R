#' ChromVAR
#'
#' Calculates and adds a chromVAR assay to the selected Seurat object.
#'
#' @param out_file Output RDS file name containing Seurat object.
#' @param so Input Seurat object. Must have a unified peaks assay to
#' run properly.
#' @param file_jaspar Path to a JASPAR PFM .txt file downloaded from
#' https://jaspar.elixir.no/downloads/.
#' @param asy Name of unified ATAC peak assay in Seurat object.
#' @param genome_use Name of genome to use for creating motif matix.
#' @param bsgenome_sel BSgenome library name (provided as object).
#' @return A Seurat object containing a chromVAR assay.
#' @import Seurat
#' @import Signac
#' @import BSgenome
#' @import TFBSTools
#' @import chromVAR
#' @import SeuratObject
#' @examples
#'
#' # d1 <- sc_chromvar(
#' #   out_file = "data_replication_cohort.rds",
#' #   so = d1,
#' #   file_jaspar = "ref/JASPAR2026_CORE_vertebrates.txt"
#' # )
#'
#' @export
sc_chromvar <- function( # nolint
  out_file,
  so,
  file_jaspar,
  asy = "ufy.peaks",
  genome_use = "hg38",
  bsgenome_sel = BSgenome.Hsapiens.UCSC.hg38 # nolint
) {
  #---- setup ----
  cat("Loading JASPAR file...", "\n")
  # Format JASPAR file
  fjas <- file_jaspar
  jaspar <- TFBSTools::readJASPARMatrix(
    fjas, matrixClass = "PFM"
  )
  cat("Loading Seurat object...", "\n")
  d1 <- so
  Seurat::DefaultAssay(d1) <- asy
  #---- create motif matrix ----
  cat("Creating motif matrix from PFM file...", "\n")
  motif.matrix <- CreateMotifMatrix( # nolint
    features = granges(d1),
    pwm = jaspar,
    genome = genome_use,
    use.counts = FALSE
  )
  cat("Adding motif data to unified peaks assay...", "\n")
  motif.object <- CreateMotifObject( # nolint
    data = motif.matrix,
    pwm = jaspar
  )
  d1 <- SetAssayData(
    d1,
    assay = asy,
    layer = "motifs",
    new.data = motif.object
  )
  #---- calculate chromVAR ----
  cat("Formatting data for chromVAR...", "\n")
  peak_counts <- GetAssayData(d1, assay = asy, layer = "counts")
  chromvar_input <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak_counts),
    rowRanges = granges(d1)
  )
  chromvar_input <- chromVAR::addGCBias(
    chromvar_input,
    genome = bsgenome_sel
  )
  motif_ix <- GetMotifData(d1, assay = asy, slot = "data")
  bg_peaks <- chromVAR::getBackgroundPeaks(chromvar_input)

  cat("Computing chromVAR...", "\n")
  dev <- chromVAR::computeDeviations(
    object = chromvar_input,
    annotations = motif_ix,
    background_peaks = bg_peaks
  )
  chromvar_scores <- chromVAR::deviationScores(dev)
  #---- output ----
  cat("Adding chromVAR assay to Seurat object...", "\n")
  d1[["chromvar"]] <- CreateAssayObject(counts = chromvar_scores)
  d1 <- SetAssayData(
    d1,
    assay = "chromvar",
    layer = "data",
    new.data = chromvar_scores
  )
  cat("Saving Seurat object...", "\n")
  saveRDS(d1, out_file)
  cat("Complete.", "\n")
  return(d1) # nolint
}
