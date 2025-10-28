#' Batch Integration of scRNA-Seq and Multiome Datasets
#'
#' Integrates a list of samples processed as Seurat objects.
#'
#' @param l_so List of processed Seurat objects to be integrated.
#' @param l_par List of processing parameters passed to function.
#' For multiome data integration, must have a gene annotation and
#' reference genome file.
#' @param proc_mode Processing mode to be used
#' (either "scRNA-Seq" or "multiome").
#' @param parl Logical indicating whether processing should be run
#' in parallel. Set to FALSE if running sequentially.
#' @param core_perc Proportion (as a numeric value) of available cores
#' to use if running in parallel. Set to 1 if running sequentially.
#' @return A Seurat object containing integrated data for all samples
#' present in a scRNA-Seq or multiome experiment.
#' @examples
#'
#' # d_integrated <- sc_integrate_data(
#' #   l_so = list_data,
#' #   l_par = l_params,
#' #   proc_mode = "multiome",
#' #   parl = TRUE,
#' #   core_perc = 0.5
#' # )
#'
#' @export
sc_integrate_data <- function( # nolint
  l_so,
  l_par,
  proc_mode,
  parl,
  core_perc
) {

  if(proc_mode == "multiome") { # nolint
    if(unlist(packageVersion("Seurat"))[1] == 5) { # nolint
      if(!file.exists("processed/data.processed.unified.rds")) { # nolint
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
        ufy_peaks <- Signac::UnifyPeaks(
          object.list = ld_int,
          mode = "reduce"
        )
        ## Filter low quality peaks
        ufy_width <- BSgenome::width(ufy_peaks)
        range(ufy_width)
        quantile(ufy_width)
        median(ufy_width)
        # Create new peak assays using unified list
        ld_int <- setNames(
          parallel::mclapply(
            mc.cores = ceiling(
              parallel::detectCores() *
                core_perc
            ),
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
        saveRDS(ld_int, "processed/data.processed.unified.rds")
        # Compute WNN for each sample, then integrate layers
        ## Run Harmony to correct for batch effects introduced by Code
        ### Remove previous objects prior to integration steps
        remove(list = ls())
        gc(reset = TRUE)
      }
      if(file.exists("processed/data.processed.unified.rds")) { # nolint
        ld_int <- readRDS("processed/data.processed.unified.rds")
      }

      # Weighted nearest neighbors (Integrate RNA and unified ATAC peak assays)
      if(!file.exists("processed/data.processed.merge.rds")) { # nolint
        ld_int <- setNames(
          parallel::mclapply(
            mc.cores = ceiling(
              parallel::detectCores() *
                0.1
            ),
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
        gc(reset = TRUE)
        library(Seurat)
        d_int <- merge(
          x = ld_int[[1]],
          y = ld_int[2:length(ld_int)],
          add.cell.ids = names(ld_int)
        )
        saveRDS(d_int, "processed/data.processed.merge.rds")
        remove(list = ls())
        gc(reset = TRUE)
      }
      if(file.exists("processed/data.processed.merge.rds")) { # nolint
        d_int <- readRDS("processed/data.processed.merge.rds")
      }

      if(file.exists("analysis/data.integrated.rds")) { # nolint
        print(
          "Integrated data file already exists!
          Use existing file or remove to re-run integration
          step from merged data"
        )
      }
      if(!file.exists("analysis/data.integrated.rds")) { # nolint
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

        ## ATAC Peaks
        Seurat::DefaultAssay(d_int) <- "ufy.peaks"
        d_int <- Signac::FindTopFeatures(d_int, min.cutoff = 5)
        d_int <- Signac::RunTFIDF(d_int)
        d_int <- Signac::RunSVD(d_int)
        Signac::DepthCor(d_int)

        ## Integrate Layers
        library(Seurat)
        library(SeuratWrappers)
        Seurat::DefaultAssay(d_int) <- "RNA"

        d_int <- Seurat::IntegrateLayers(
          object = d_int,
          method = CCAIntegration, # nolint
          orig.reduction = "pca",
          new.reduction = "int.SCT.CCA",
          verbose = TRUE
        )
        d_int <- SeuratObject::JoinLayers(d_int)

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
          resolution = 0.5
        )

        ## Visualize Clusters
        ggplot2::ggsave(
          "analysis/plot.umap.panel.WNN.png",
          sc_umap_panel( # nolint
            d_int,
            c("Group", "Code", "seurat_clusters"),
            "wnn.umap"
          ),
          height = 8,
          width = 24,
          dpi = 600
        )

        ## Integration Performance
        d1qc_pre <- Seurat::VlnPlot(
          object = d_int,
          features = c(
            "nFeature_RNA",
            "nCount_RNA",
            "percent.mt",
            "nCount_ATAC",
            "TSS.enrichment",
            "nucleosome_signal"
          ),
          layer = "counts",
          ncol = 4,
          pt.size = 0.2
        )
        ggplot2::ggsave(
          "analysis/plot.integration.qc.png",
          d1qc_pre,
          width = 24,
          height = 8,
          dpi = 800
        )
        ## Save
        saveRDS(d_int, "analysis/data.integrated.rds")
      }
    }
  }

  if(proc_mode == "scRNA-Seq") { # nolint
    if(unlist(packageVersion("Seurat"))[1] == 5) { # nolint
      # Weighted nearest neighbors (Integrate RNA assays)
      if(!file.exists("processed/data.processed.merge.rds")) { # nolint
        ld_int <- l_so
        library(Seurat)
        d_int <- merge(
          x = ld_int[[1]],
          y = ld_int[2:length(ld_int)],
          add.cell.ids = names(ld_int)
        )
        saveRDS(d_int, "processed/data.processed.merge.rds")
        remove(list = ls())
        gc(reset = TRUE)
      }
      if(file.exists("processed/data.processed.merge.rds")) { # nolint
        d_int <- readRDS("processed/data.processed.merge.rds")
      }

      if(!file.exists("analysis/data.integrated.rds")) { # nolint
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

        ## Integrate Layers
        library(Seurat)
        library(SeuratWrappers)
        Seurat::DefaultAssay(d_int) <- "RNA"

        d_int <- Seurat::IntegrateLayers(
          object = d_int,
          method = CCAIntegration, # nolint
          orig.reduction = "pca",
          new.reduction = "int.SCT.CCA",
          verbose = TRUE
        )
        d_int <- SeuratObject::JoinLayers(d_int)

        ## Harmony batch effect correction for GEX
        d_int <- harmony::RunHarmony(
          d_int,
          assay.use = "RNA",
          group.by.vars = "Code",
          reduction.use = "int.SCT.CCA",
          reduction.save = "RNA.corrected",
          project.dim = FALSE
        )

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
          resolution = 0.5
        )

        ## Visualize Clusters
        ggplot2::ggsave(
          "analysis/plot.umap.panel.WNN.png",
          sc_umap_panel( # nolint
            d_int,
            c("Group", "Code", "seurat_clusters"),
            "int.umap"
          ),
          height = 8,
          width = 24,
          dpi = 600
        )

        ## Integration Performance
        d1qc_pre <- Seurat::VlnPlot(
          object = d_int,
          features = c(
            "nFeature_RNA",
            "nCount_RNA",
            "percent.mt"
          ),
          layer = "counts",
          ncol = 4,
          pt.size = 0.2
        )
        ggplot2::ggsave(
          "analysis/plot.integration.qc.png",
          d1qc_pre,
          width = 24,
          height = 8,
          dpi = 800
        )
        ## Save
        saveRDS(d_int, "analysis/data.integrated.rds")
      }
      if(file.exists("analysis/data.integrated.rds")) { # nolint
        print(
          "Integrated data file already exists!
          Use existing file or remove to re-run integration
          step from merged data"
        )
      }
    }
  }

  if(proc_mode == "scRNA-Seq") { # nolint
    if(unlist(packageVersion("Seurat"))[1] < 5) { # nolint
      if(Sys.info()[["sysname"]] != "Windows" && parl == TRUE) { # nolint
        fun_int <- function(list_d_subset, future_size) {
          ## Find integration anchors
          future::plan(
            "multisession",
            workers = ceiling(parallel::detectCores() * core_perc)
          )
          options(future.globals.maxSize = future_size * 1024^2)
          d_anchor <- Seurat::FindIntegrationAnchors(
            object.list = list_d_subset,
            anchor.features = 4000,
            dims = 1:50
          )
          ## Integrate data
          future::plan("sequential")
          options(future.globals.maxSize = 500 * 1024^2)
          d_merged <- Seurat::IntegrateData(
            anchorset = d_anchor,
            dims = 1:50
          )
          return(d_merged) # nolint
        }
        if(length(l_so) <= 12 & length(l_so) > 9) { # nolint
          print(
            paste(
              length(l_so),
              "samples present: dividing integration into 4 batches...",
              sep = " "
            )
          )
          # Batch 1
          d_merged1 <- fun_int(l_so[1:3], 5000)
          saveRDS(
            d_merged1,
            "processed/data_integrated_batch_1.rds"
          )
          # Batch 2
          d_merged2 <- fun_int(l_so[4:6], 5000)
          saveRDS(
            d_merged2,
            "processed/data_integrated_batch_2.rds"
          )
          # Batch 3
          d_merged3 <- fun_int(l_so[7:9], 5000)
          saveRDS(
            d_merged3,
            "processed/data_integrated_batch_3.rds"
          )
          # Batch 4
          d_merged4 <- fun_int(l_so[10:length(l_so)], 5000)
          saveRDS(
            d_merged4,
            "processed/data_integrated_batch_4.rds"
          )
          ## Combine batch 1&2 then 3&4
          d_merged12 <- fun_int(list(d_merged1, d_merged2), 10000)
          saveRDS(
            d_merged12,
            "processed/data_integrated_batch_1and2.rds"
          )
          d_merged34 <- fun_int(list(d_merged3, d_merged4), 10000)
          saveRDS(
            d_merged34,
            "processed/data_integrated_batch_3and4.rds"
          )
          remove(
            d_merged1,
            d_merged2,
            d_merged3,
            d_merged4
          )
          ## Combine into final batch and save
          d_integrated <- fun_int(list(d_merged12, d_merged34), 15000)
          saveRDS(
            d_integrated,
            "processed/data_integrated.rds"
          )
          remove(d_merged12, d_merged34)
          file.remove("processed/data_integrated_batch_1.rds")
          file.remove("processed/data_integrated_batch_2.rds")
          file.remove("processed/data_integrated_batch_3.rds")
          file.remove("processed/data_integrated_batch_4.rds")
          file.remove("processed/data_integrated_batch_1and2.rds")
          file.remove("processed/data_integrated_batch_3and4.rds")
        }

        if(length(l_so) <= 9 & length(l_so) > 6) { # nolint
          print(
            paste(
              length(l_so),
              "samples present: dividing integration into 3 batches...",
              sep = " "
            )
          )
          # Batch 1
          d_merged1 <- fun_int(l_so[1:3], 5000)
          saveRDS(
            d_merged1,
            "processed/data_integrated_batch_1.rds"
          )
          # Batch 2
          d_merged2 <- fun_int(l_so[4:6], 5000)
          saveRDS(
            d_merged2,
            "processed/data_integrated_batch_2.rds"
          )
          # Batch 3
          d_merged3 <- fun_int(l_so[7:length(l_so)], 5000)
          saveRDS(
            d_merged3,
            "processed/data_integrated_batch_3.rds"
          )
          ## Combine batch 1&2 then 3
          d_merged12 <- fun_int(list(d_merged1, d_merged2), 10000)
          saveRDS(
            d_merged12,
            "processed/data_integrated_batch_1and2.rds"
          )
          remove(d_merged1, d_merged2)
          ## Combine into final batch and save
          d_integrated <- fun_int(list(d_merged12, d_merged3), 15000)
          saveRDS(
            d_integrated,
            "processed/data_integrated.rds"
          )
          remove(d_merged12, d_merged3)
          file.remove("processed/data_integrated_batch_1.rds")
          file.remove("processed/data_integrated_batch_2.rds")
          file.remove("processed/data_integrated_batch_3.rds")
          file.remove("processed/data_integrated_batch_1and2.rds")
        }

        if(length(l_so) <= 6 && length(l_so) > 3) { # nolint
          print(
            paste(
              length(l_so),
              "samples present: dividing integration into 2 batches...",
              sep = " "
            )
          )
          # Batch 1
          d_merged1 <- fun_int(l_so[1:3], 5000)
          saveRDS(
            d_merged1,
            "processed/data_integrated_batch_1.rds"
          )
          # Batch 2
          d_merged2 <- fun_int(l_so[4:6], 5000)
          saveRDS(
            d_merged2,
            "processed/data_integrated_batch_2.rds"
          )
          ## Combine batch 1&2
          d_integrated <- fun_int(list(d_merged1, d_merged2), 10000)
          saveRDS(
            d_integrated,
            "processed/data_integrated.rds"
          )
          remove(d_merged1, d_merged2)
          file.remove("processed/data_integrated_batch_1.rds")
          file.remove("processed/data_integrated_batch_2.rds")
        }

        if(length(l_so) <= 3) { # nolint
          print(
            paste(
              length(l_so),
              "samples present: integrating all samples as 1 batch...",
              sep = " "
            )
          )
          # Batch 1
          d_integrated <- fun_int(l_so[1:length(l_so)], 5000) # nolint
          saveRDS(
            d_integrated,
            "processed/data_integrated.rds"
          )
        }
      }
    }

    if(Sys.info()[["sysname"]] == "Windows") { # nolint
      fun_int <- function(list_d_subset, future_size) {
        ## Find integration anchors
        d_anchor <- Seurat::FindIntegrationAnchors(
          object.list = list_d_subset,
          anchor.features = 4000,
          dims = 1:50
        )
        ## Integrate data
        d_merged <- Seurat::IntegrateData(
          anchorset = d_anchor,
          dims = 1:50
        )
        return(d_merged) # nolint
      }

      if(length(l_so) <= 12 && length(l_so) > 9) { # nolint
        print(
          paste(
            length(l_so),
            "samples present: dividing integration into 4 batches...",
            sep = " "
          )
        )
        # Batch 1
        d_merged1 <- fun_int(l_so[1:3], 5000)
        saveRDS(
          d_merged1,
          "processed/data_integrated_batch_1.rds"
        )
        # Batch 2
        d_merged2 <- fun_int(l_so[4:6], 5000)
        saveRDS(
          d_merged2,
          "processed/data_integrated_batch_2.rds"
        )
        # Batch 3
        d_merged3 <- fun_int(l_so[7:9], 5000)
        saveRDS(
          d_merged3,
          "processed/data_integrated_batch_3.rds"
        )
        # Batch 4
        d_merged4 <- fun_int(l_so[10:length(l_so)], 5000)
        saveRDS(
          d_merged4,
          "processed/data_integrated_batch_4.rds"
        )
        ## Combine batch 1&2 then 3&4
        d_merged12 <- fun_int(list(d_merged1, d_merged2), 10000)
        saveRDS(
          d_merged12,
          "processed/data_integrated_batch_1and2.rds"
        )
        d_merged34 <- fun_int(list(d_merged3, d_merged4), 10000)
        saveRDS(
          d_merged34,
          "processed/data_integrated_batch_3and4.rds"
        )
        remove(d_merged1, d_merged2, d_merged3, d_merged4)
        ## Combine into final batch and save
        d_integrated <- fun_int(list(d_merged12, d_merged34), 15000)
        saveRDS(
          d_integrated,
          "processed/data_integrated.rds"
        )
        remove(d_merged12, d_merged34)
        file.remove("processed/data_integrated_batch_1.rds")
        file.remove("processed/data_integrated_batch_2.rds")
        file.remove("processed/data_integrated_batch_3.rds")
        file.remove("processed/data_integrated_batch_4.rds")
        file.remove("processed/data_integrated_batch_1and2.rds")
        file.remove("processed/data_integrated_batch_3and4.rds")
      }

      if(length(l_so) <= 9 && length(l_so) > 6) { # nolint
        print(
          paste(
            length(l_so),
            "samples present: dividing integration into 3 batches...",
            sep = " "
          )
        )
        # Batch 1
        d_merged1 <- fun_int(l_so[1:3], 5000)
        saveRDS(
          d_merged1,
          "processed/data_integrated_batch_1.rds"
        )
        # Batch 2
        d_merged2 <- fun_int(l_so[4:6], 5000)
        saveRDS(
          d_merged2,
          "processed/data_integrated_batch_2.rds"
        )
        # Batch 3
        d_merged3 <- fun_int(l_so[7:length(l_so)], 5000)
        saveRDS(
          d_merged3,
          "processed/data_integrated_batch_3.rds"
        )
        ## Combine batch 1&2 then 3
        d_merged12 <- fun_int(list(d_merged1, d_merged2), 10000)
        saveRDS(
          d_merged12,
          "processed/data_integrated_batch_1and2.rds"
        )
        remove(d_merged1, d_merged2)
        ## Combine into final batch and save
        d_integrated <- fun_int(list(d_merged12, d_merged3), 15000)
        saveRDS(
          d_integrated,
          "processed/data_integrated.rds"
        )
        remove(d_merged12, d_merged3)
        file.remove("processed/data_integrated_batch_1.rds")
        file.remove("processed/data_integrated_batch_2.rds")
        file.remove("processed/data_integrated_batch_3.rds")
        file.remove("processed/data_integrated_batch_1and2.rds")
      }

      if(length(l_so) <= 6 && length(l_so) > 3) { # nolint
        print(
          paste(
            length(l_so),
            "samples present: dividing integration into 2 batches...",
            sep = " "
          )
        )
        # Batch 1
        d_merged1 <- fun_int(l_so[1:3], 5000)
        saveRDS(
          d_merged1,
          "processed/data_integrated_batch_1.rds"
        )
        # Batch 2
        d_merged2 <- fun_int(l_so[4:6], 5000)
        saveRDS(
          d_merged2,
          "processed/data_integrated_batch_2.rds"
        )
        ## Combine batch 1&2
        d_integrated <- fun_int(list(d_merged1, d_merged2), 10000)
        saveRDS(
          d_integrated,
          "processed/data_integrated.rds"
        )
        remove(d_merged1, d_merged2)
        file.remove("processed/data_integrated_batch_1.rds")
        file.remove("processed/data_integrated_batch_2.rds")
      }

      if(length(l_so) <= 3) { # nolint
        print(
          paste(
            length(l_so),
            "samples present: integrating all samples as 1 batch...",
            sep = " "
          )
        )
        # Batch 1
        d_integrated <- fun_int(l_so[1:length(l_so)], 5000) # nolint
        saveRDS(
          d_integrated,
          "processed/data_integrated.rds"
        )
      }
    }
    remove(l_so)
    return(d_integrated) # nolint
  }
}


#' scRNA-Seq Integration Quality Control
#'
#' Displays the feature number, average counts, and percentage of reads
#' derived from the mitochondrial genome given an integrated Seurat object.
#'
#' @param so Integrated Seurat object.
#' @param cl_var Clustering variable provided as a
#' character string (generally use "seurat_clusters").
#' @return A panel of violin plots providing QC
#' measures for an integrated Seurat object.
#' @examples
#'
#' # d.integrated <- sc_integrate_data(d.integrated)
#'
#' @export
sc_integration_qc <- function(
  so,
  cl_var
) {
  d <- so
  df <- d@meta.data[, c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    cl_var
  )]
  p_qc <- ggpubr::ggarrange(
    plotlist = lapply(
      names(
        dplyr::select(
          df,
          -.data[[cl_var]] # nolint
        )
      ),
      function(x) {
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            x = .data[[cl_var]], # nolint
            y = .data[[x]],
            fill = .data[[cl_var]]
          )
        ) +
          ggplot2::scale_fill_manual(
            name = "Cluster",
            values = col_univ() # nolint
          ) +
          # Add violin plot and dotplot
          ggplot2::geom_violin(
            trim = FALSE
          ) +
          ggplot2::geom_jitter(
            ggplot2::aes(
              alpha = 0.2
            ),
            shape = 16,
            size = 0.1,
            position = ggplot2::position_jitter(
              width = 0.4
            ),
            show.legend = FALSE
          ) +
          # Add Theme
          sc_theme1() + # nolint
          ggplot2::labs(y = x) +
          ggplot2::theme(
            plot.margin = ggplot2::unit(
              c(0.1, 0.1, 0.1, 0.1),
              "cm"
            )
          )
      }
    ),
    common.legend = TRUE,
    legend = "none",
    nrow = 1,
    ncol = 3
  )
  return(p_qc) # nolint
}