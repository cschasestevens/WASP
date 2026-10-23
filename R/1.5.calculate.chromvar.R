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
#' @param slt Slot from the chosen assay to use for filtering out
#' empty peaks.
#' @param bsgenome_sel BSgenome library name (provided as object).
#' 
#' @return A Seurat object containing a chromVAR assay.
#' @import Seurat
#' @import Signac
#' @import BSgenome
#' @import TFBSTools
#' @import chromVAR
#' @import SeuratObject
#' @import motifmatchr
#' @import BiocParallel
#' @import GenomicRanges
#' @import Matrix
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
  slt = "counts",
  bsgenome_sel = BSgenome.Hsapiens.UCSC.hg38, # nolint
  cores = 12
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
  cat("Filtering out empty peaks from assay...", "\n")
  filt1 <- Seurat::GetAssayData(d1, assay = asy, layer = slt)
  filt1 <- rowSums(filt1) > 0
  d1 <- d1[filt1, ]
  if (Sys.info()[["sysname"]] != "Windows" && cores > 1) {
    # set parallel parameters
    BiocParallel::register(BiocParallel::MulticoreParam(workers = cores))
    peaks <- Signac::granges(d1[["ufy.peaks"]])
    chunkn <- cores
    # split peaks into even chunks
    chunks <- split(peaks, cut(seq_along(peaks), chunkn, labels = FALSE))
    cat("Scanning motifs in parallel using specified number of cores...\n") # nolint
    # match motifs for each chunk
    chunked_results <- BiocParallel::bplapply(chunks, function(chunk_peaks) {
      # Inside bplapply, matchMotifs runs cleanly on its specific sub-chunk
      motifmatchr::matchMotifs(
        pwms = jaspar,
        subject = chunk_peaks,
        genome = bsgenome_sel, # nolint
        out = "matches"
      )
    })
    # combine and create motif object to add to existing peaks assay
    combined_matches <- do.call(
      rbind,
      lapply(chunked_results, motifmatchr::motifMatches)
    )
    rownames(combined_matches) <- rownames(d1[["ufy.peaks"]])
    motif_obj <- Signac::CreateMotifObject(
      data = combined_matches, pwm = jaspar
    )
    d1[[asy]] <- SeuratObject::SetAssayData(
      object = d1[[asy]],
      slot = "motifs",
      new.data = motif_obj
    )
  } else {
    # set parameters
    peaks <- Signac::granges(d1[["ufy.peaks"]])
    chunkn <- cores
    # split peaks into even chunks
    chunks <- split(peaks, cut(seq_along(peaks), chunkn, labels = FALSE))
    cat("Scanning motifs sequentially...\n") # nolint
    # match motifs for each chunk
    chunked_results <- lapply(chunks, function(chunk_peaks) {
      # Inside bplapply, matchMotifs runs cleanly on its specific sub-chunk
      motifmatchr::matchMotifs(
        pwms = jaspar,
        subject = chunk_peaks,
        genome = bsgenome_sel, # nolint
        out = "matches"
      )
    })
    # combine and create motif object to add to existing peaks assay
    combined_matches <- do.call(
      rbind,
      lapply(chunked_results, motifmatchr::motifMatches)
    )
    rownames(combined_matches) <- rownames(d1[["ufy.peaks"]])
    motif_obj <- Signac::CreateMotifObject(
      data = combined_matches, pwm = jaspar
    )
    d1[[asy]] <- SeuratObject::SetAssayData(
      object = d1[[asy]],
      slot = "motifs",
      new.data = motif_obj
    )
  }
  #---- calculate chromVAR ----
  cat("Formatting data for chromVAR...", "\n")
  peak_counts <- SeuratObject::GetAssayData(
    d1, assay = asy, layer = "counts"
  )
  chromvar_input <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak_counts),
    rowRanges = granges(d1)
  )
  # chromVAR background peaks
  chromvar_input <- chromVAR::addGCBias(
    chromvar_input,
    genome = bsgenome_sel
  )
  motif_ix <- Signac::GetMotifData(d1, assay = asy, slot = "data")
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
  d1[["chromvar"]] <- SeuratObject::CreateAssayObject(counts = chromvar_scores)
  Seurat::DefaultAssay(d1) <- "chromvar"
  cat("Saving Seurat object...", "\n")
  saveRDS(d1, out_file)
  cat("Complete.", "\n")
  return(d1) # nolint
}
