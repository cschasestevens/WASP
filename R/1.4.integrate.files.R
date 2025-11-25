#' Batch Integration of scRNA-Seq and Multiome Datasets
#'
#' Integrates a list of samples processed as Seurat objects.
#'
#' @param l_so List of processed Seurat objects to be integrated
#' from any sc_process function.
#' @param l_par List of processing parameters obtained from a
#' sc_params function.
#' @param proc_mode Processing mode to be used
#' (either "scRNA-Seq" or "multiome").
#' @param parl Logical indicating whether processing should be run
#' in parallel. Set to FALSE if running sequentially.
#' @param core_num Number of available cores
#' to use if running in parallel.
#' @param res_clus Seurat clustering resolution.
#' @return A Seurat object containing integrated data for all samples
#' present in a scRNA-Seq or multiome experiment.
#' @import Seurat
#' @import Signac
#' @import BSgenome
#' @import parallel
#' @import harmony
#' @import SeuratObject
#' @examples
#'
#' # d_integrated <- sc_integrate_data(
#' #   l_so = scdata,
#' #   l_par = scparam,
#' #   proc_mode = "multiome"
#' # )
#'
#' @export
sc_integrate_data <- function( # nolint
  l_so,
  l_par,
  proc_mode,
  parl = FALSE,
  core_num = NULL,
  res_clus = 0.5
) {
  # Standard multiome data integration
  if(proc_mode == "multiome") { # nolint
    print("Integrating samples using multiome method...")
    # Integration for Seurat versions 5*
    if (unlist(packageVersion("Seurat"))[1] < 5) {
      print("Integrating data using Seurat v4 and older was deprecated in WASP 4.04; Upgrade Seurat version and re-run data integration.") # nolint
    }
    # For Seurat versions 5*
    if(unlist(packageVersion("Seurat"))[1] == 5) { # nolint
      if (file.exists("processing/int_unified.rds") && !file.exists("processing/int_merged.rds")) { # nolint
        print("Unified peaks Seurat object already exists; loading existing object...") # nolint
        ld_int <- readRDS("processing/int_unified.rds")
      }
      if(!file.exists("processing/int_unified.rds")) { # nolint
        print("---- Step 1: Pull processed objects from list ----")
        ld_int <- l_so
        # Extract Seurat objects for integration
        ld_int <- setNames(
          lapply(
            seq.int(1, length(ld_int), 1),
            function(x) {
              d1 <- ld_int[[x]][[1]]
              Seurat::DefaultAssay(d1) <- "SCT"
              d1 <- Seurat::RunPCA(d1)
              Seurat::DefaultAssay(d1) <- "peaks"
              return(d1) # nolint
            }
          ),
          names(ld_int)
        )
        # Unify ATAC peaks across data sets
        print("---- Step 2: Unify ATAC peaks ----")
        ufy_peaks <- Signac::UnifyPeaks(
          object.list = ld_int,
          mode = "reduce"
        )
        ## Filter low quality peaks
        ufy_width <- BSgenome::width(ufy_peaks)
        print(range(ufy_width))
        print(quantile(ufy_width))
        print(median(ufy_width))
        # Create new peak assays using unified list
        if (Sys.info()[["sysname"]] != "Windows" && parl == TRUE) {
          ld_int <- setNames(
            parallel::mclapply(
              mc.cores = core_num,
              seq.int(1, length(ld_int), 1),
              function(x) {
                d1 <- ld_int[[x]]
                Seurat::DefaultAssay(d1) <- "peaks"
                d1_cnts <- Signac::FeatureMatrix(
                  fragments = Signac::Fragments(d1),
                  features = ufy_peaks,
                  cells = colnames(d1)
                )
                d1[["ufy.peaks"]] <- Signac::CreateChromatinAssay(
                  counts = d1_cnts,
                  fragments = Signac::Fragments(d1),
                  annotation = l_par[["ref.gtf"]]
                )
                Seurat::DefaultAssay(d1) <- "ufy.peaks"
                return(d1) # nolint
              }
            ),
            names(ld_int)
          )
        }
        if (Sys.info()[["sysname"]] == "Windows" || parl == FALSE) {
          ld_int <- setNames(
            lapply(
              seq.int(1, length(ld_int), 1),
              function(x) {
                d1 <- ld_int[[x]]
                Seurat::DefaultAssay(d1) <- "peaks"
                d1_cnts <- Signac::FeatureMatrix(
                  fragments = Signac::Fragments(d1),
                  features = ufy_peaks,
                  cells = colnames(d1)
                )
                d1[["ufy.peaks"]] <- Signac::CreateChromatinAssay(
                  counts = d1_cnts,
                  fragments = Signac::Fragments(d1),
                  annotation = l_par[["ref.gtf"]]
                )
                Seurat::DefaultAssay(d1) <- "ufy.peaks"
                return(d1) # nolint
              }
            ),
            names(ld_int)
          )
        }
        saveRDS(ld_int, "processing/int_unified.rds")
        # Compute WNN for each sample, then integrate layers
        ## Run Harmony to correct for batch effects introduced by Code
        ### Remove previous objects prior to integration steps
        gc(reset = TRUE)
      }
      # Weighted nearest neighbors (Integrate RNA and unified ATAC peak assays)
      if (file.exists("processing/int_merged.rds")) {
        print("Merged Seurat object already exists; loading existing object...") # nolint
        d_int <- readRDS("processing/int_merged.rds")
      }
      if(!file.exists("processing/int_merged.rds")) { # nolint
        print("---- Step 3: WNN integration of GEX/ATAC assays ----")
        if (Sys.info()[["sysname"]] != "Windows" && parl == TRUE) {
          ld_int <- setNames(
            parallel::mclapply(
              mc.cores = 2,
              seq.int(1, length(ld_int), 1),
              function(x) {
                d1 <- ld_int[[x]]
                d1 <- Seurat::FindMultiModalNeighbors(
                  d1,
                  reduction.list = list("pca", "lsi"),
                  dims.list = list(1:50, 2:50)
                )
                Seurat::DefaultAssay(d1) <- "ufy.peaks"
                return(d1) # nolint
              }
            ),
            names(ld_int)
          )
        }
        if (Sys.info()[["sysname"]] == "Windows" || parl == FALSE) {
          ld_int <- setNames(
            lapply(
              seq.int(1, length(ld_int), 1),
              function(x) {
                d1 <- ld_int[[x]]
                d1 <- Seurat::FindMultiModalNeighbors(
                  d1,
                  reduction.list = list("pca", "lsi"),
                  dims.list = list(1:50, 2:50)
                )
                Seurat::DefaultAssay(d1) <- "ufy.peaks"
                return(d1) # nolint
              }
            ),
            names(ld_int)
          )
        }
        gc(reset = TRUE)
        print("---- Step 4: Merge datasets ----")
        if (length(ld_int) == 1) {
          print("Only a single object is present in the dataset; skipping merge step...") # nolint
          d_int <- ld_int[[1]]
          saveRDS(d_int, "processing/int_merged.rds")
        }
        if (length(ld_int) > 1) {
          d_int <- merge(
            x = ld_int[[1]],
            y = ld_int[2:length(ld_int)],
            add.cell.ids = names(ld_int)
          )
          saveRDS(d_int, "processing/int_merged.rds")
        }
        remove(ld_int)
        gc(reset = TRUE)
      }
      print("---- Step 5: Normalize data ----")
      # Normalization
      ## GEX
      Seurat::DefaultAssay(d_int) <- "RNA"
      d_int <- Seurat::FindVariableFeatures(
        d_int,
        selection.method = "vst",
        nfeatures = 3000
      )
      d_int <- Seurat::NormalizeData(d_int)
      d_int <- Seurat::ScaleData(d_int)
      d_int <- Seurat::RunPCA(d_int)
      d_int <- Seurat::FindVariableFeatures(d_int, nfeatures = 3000)
      ## ATAC
      Seurat::DefaultAssay(d_int) <- "ufy.peaks"
      d_int <- Signac::FindTopFeatures(d_int, min.cutoff = 5)
      d_int <- Signac::RunTFIDF(d_int)
      d_int <- Signac::RunSVD(d_int)
      print("---- Step 6: Integrate layers ----")
      ## Integrate Layers
      Seurat::DefaultAssay(d_int) <- "SCT"
      library(Seurat)
      d_int <- Seurat::IntegrateLayers(
        object = d_int,
        method = CCAIntegration, # nolint
        assay = "SCT",
        orig.reduction = "pca",
        new.reduction = "int.SCT.CCA",
        verbose = TRUE,
        normalization.method = "SCT"
      )
      Seurat::DefaultAssay(d_int) <- "RNA"
      d_int <- SeuratObject::JoinLayers(d_int)
      print("---- Step 7: Harmony batch effect correction ----")
      ## Harmony batch effect correction for GEX and ATAC
      d_int <- harmony::RunHarmony(
        d_int,
        assay.use = "ufy.peaks",
        group.by.vars = "Code",
        reduction.use = "lsi",
        reduction.save = "ufy.peaks.corrected",
        project.dim = FALSE
      )
      d_int <- harmony::RunHarmony(
        d_int,
        assay.use = "RNA",
        group.by.vars = "Code",
        reduction.use = "int.SCT.CCA",
        reduction.save = "RNA.corrected",
        project.dim = FALSE
      )
      print("---- Step 8: WNN for final integration of multimodal data ----") # nolint
      ## WNN
      d_int <- Seurat::FindMultiModalNeighbors(
        d_int,
        reduction.list = list("RNA.corrected", "ufy.peaks.corrected"),
        dims.list = list(1:50, 2:50)
      )
      d_int <- Seurat::RunUMAP(
        d_int,
        nn.name = "weighted.nn",
        reduction.name = "wnn.umap",
        reduction.key = "wnnUMAP_",
        n.components = 3
      )
      d_int <- Seurat::FindClusters(
        d_int,
        graph.name = "wsnn",
        algorithm = 3,
        verbose = TRUE,
        resolution = res_clus
      )
      print("Integration complete!")
      return(d_int)
    }
  }
  # Standard scRNA-Seq data integration
  if(proc_mode == "scRNA-Seq") { # nolint
    print("Integrating samples using scRNA-Seq method...")
    # Integration for Seurat versions 5*
    if (unlist(packageVersion("Seurat"))[1] < 5) {
      print("Integrating data using Seurat v4 and older was deprecated in WASP 4.04; Upgrade Seurat version and re-run data integration.") # nolint
    }
    if(unlist(packageVersion("Seurat"))[1] == 5) { # nolint
      # Weighted nearest neighbors (Integrate RNA assays)
      print("---- Step 1: Pull processed objects from list ----")
      if (file.exists("processing/int_merged.rds")) {
        print("Merged Seurat object already exists; loading existing object...") # nolint
        ld_int <- readRDS("processing/int_merged.rds")
      }
      if(!file.exists("processing/int_merged.rds")) { # nolint
        ld_int <- l_so
        print("---- Step 2: Merge objects ----")
        if (length(ld_int) == 1) {
          print("Only a single object is present in the dataset; skipping merge step...") # nolint
          d_int <- ld_int[[1]]
          saveRDS(d_int, "processing/int_merged.rds")
        }
        if (length(ld_int) > 1) {
          d_int <- merge(
            x = ld_int[[1]],
            y = ld_int[2:length(ld_int)],
            add.cell.ids = names(ld_int)
          )
          saveRDS(d_int, "processing/int_merged.rds")
        }
        remove(ld_int)
        gc(reset = TRUE)
      }
      print("---- Step 3: Normalize data ----")
      # Normalization
      ## GEX
      Seurat::DefaultAssay(d_int) <- "RNA"
      d_int <- Seurat::FindVariableFeatures(
        d_int,
        selection.method = "vst",
        nfeatures = 3000
      )
      d_int <- Seurat::NormalizeData(d_int)
      d_int <- Seurat::ScaleData(d_int)
      d_int <- Seurat::RunPCA(d_int)
      d_int <- Seurat::FindVariableFeatures(d_int, nfeatures = 3000)
      print("---- Step 4: Integrate layers ----")
      ## Integrate Layers
      Seurat::DefaultAssay(d_int) <- "SCT"
      library(Seurat)
      d_int <- Seurat::IntegrateLayers(
        object = d_int,
        method = CCAIntegration, # nolint
        orig.reduction = "pca",
        new.reduction = "int.SCT.CCA",
        verbose = TRUE,
        normalization.method = "SCT"
      )
      Seurat::DefaultAssay(d_int) <- "RNA"
      d_int <- SeuratObject::JoinLayers(d_int)
      print("---- Step 5: Harmony batch effect correction ----")
      ## Harmony batch effect correction for GEX
      d_int <- harmony::RunHarmony(
        d_int,
        assay.use = "RNA",
        group.by.vars = "Code",
        reduction.use = "int.SCT.CCA",
        reduction.save = "RNA.corrected",
        project.dim = FALSE
      )
      print("---- Step 6: Find clusters and output data ----")
      d_int <- Seurat::RunUMAP(
        d_int,
        reduction = "int.SCT.CCA",
        reduction.name = "int.umap",
        reduction.key = "wnnUMAP_",
        n.components = 3
      )
      d_int <- Seurat::FindClusters(
        d_int,
        graph.name = "int.SCT",
        algorithm = 3,
        verbose = TRUE,
        resolution = res_clus
      )
      print("Integration complete!")
      return(d_int)
    }
  }
}
