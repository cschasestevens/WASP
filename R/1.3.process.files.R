#' Process CellRanger Data Files
#'
#' Processes a single sample into a Seurat object for data integration.
#'
#' @param df_par Data frame of parameters to use for data processing.
#' @param i Numeric value for sample ID.
#' @param rho_adj Numeric value indicating a proportion to scale
#' rho values calculated during the ambient RNA removal step. 10% (0.1)
#' is generally suitable for most data sets.
#' @param m_cell Threshold for including features in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @param m_feat Threshold for including cells in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @return A processed sample file converted into a Seurat object
#' with a summary list of QC and processing details.
#' @examples
#'
#' # proc_data <- sc_process_file(
#' #   df_par = list_params,
#' #   i = 1
#' # )
#'
#' @export
sc_process_file <- function(
  df_par,
  i,
  rho_adj = 0.1,
  m_cell = 5,
  m_feat = 200
) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1234)
  # Select file from chosen input parameter df
  df_p <- df_par
  d <- df_p[i, ]
  # Estimate contamination fraction
  d <- SoupX::autoEstCont(
    SoupX::load10X(
      df_p[i, "Path"]
    ),
    doPlot = FALSE
  )
  df_p[["rho"]] <- as.vector(
    unlist(
      unique(
        d[["metaData"]][["rho"]]
      )
    )
  )
  ## Adjusted rho
  df_p[["adj.rho"]] <- df_p[["rho"]] +
    rho_adj *
      df_p[["rho"]]
  ## Cluster and cell numbers
  df_p[["Clusters"]] <- length(
    unique(
      (d)$metaData$clusters
    )
  )
  df_p[["Cell.No"]] <- nrow(d$metaData)
  # Adjust contamination fraction by selected adj.rho value
  d <- SoupX::adjustCounts(
    SoupX::setContaminationFraction(
      d,
      df_p[["adj.rho"]]
    )
  )
  # Doublet removal/Create single cell experiment
  d <- scDblFinder::scDblFinder(
    SingleCellExperiment::SingleCellExperiment(
      list(
        counts = as.matrix(d)
      )
    )
  )
  d <- SummarizedExperiment::assay(
    d[, d$scDblFinder.class == "singlet"],
    "counts"
  )
  # Create Seurat object
  d <- Seurat::CreateSeuratObject(
    counts = d,
    project = paste(df_p[i, "File.ID"]),
    min.cells = m_cell,
    min.features = m_feat
  )
  d_md <- df_p[i, 5:ncol(df_p)]
  d_md <- setNames(
    as.data.frame(
      lapply(
        seq.int(
          1,
          ncol(d_md),
          1
        ),
        function(x) {
          rep(
            d_md[1, x],
            nrow(d@meta.data)
          )
        }
      )
    ),
    names(d_md)
  )
  d <- Seurat::AddMetaData(
    object = d,
    metadata = d_md
  )
  # Add gene names
  d_g <- setNames(
    read.table(
      df_p[i,
           c("Path.feat")],
      sep = "\t",
      header = FALSE
    ),
    c("GENEID", "GENE", "TYPE")
  )
  d_g <- d_g[!duplicated(d_g[["GENE"]]), ]
  rownames(d_g) <- d_g[["GENE"]]
  d[["RNA"]] <- Seurat::AddMetaData(
    object = d[["RNA"]],
    metadata = d_g
  )
  # Subset Seurat object
  d[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object = d,
    pattern = "^MT-"
  )
  plot_pre_qc <- Seurat::VlnPlot(
    object = d,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt"),
    ncol = 3,
    pt.size = 0.2
  )
  d_filt <- BiocGenerics::subset(
    d,
    subset = nFeature_RNA > 300 & # nolint
      nFeature_RNA < 7000 &  # nolint
      percent.mt < 20 # nolint
  )
  plot_pos_qc <- Seurat::VlnPlot(
    object = d_filt,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt"),
    ncol = 3,
    pt.size = 0.2
  )
  # Log normalize data
  d_norm <- Seurat::NormalizeData(
    object = d_filt,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  d_norm <- Seurat::FindVariableFeatures(
    object = d_norm,
    selection.method = "vst",
    nfeatures = 2000
  )
  d_norm_sum <- head(
    Seurat::VariableFeatures(d_norm),
    25
  )
  plot_var_feat <- Seurat::VariableFeaturePlot(
    d_norm,
    assay = "RNA"
  )
  plot_var_feat_out <- Seurat::LabelPoints(
    plot_var_feat,
    points = d_norm_sum,
    repel = TRUE
  )
  df_p <- df_p[i, ]

  return(
    list(
      "Seurat.obj" = d_norm,
      "Params" = df_p,
      "QC.pre" = plot_pre_qc,
      "QC.post" = plot_pos_qc,
      "var.feat.sum" = d_norm_sum,
      "var.feat.plot" = plot_var_feat_out
    )
  )
}

#' Batch Processing of CellRanger Files
#'
#' Processes a list of scRNA-Seq samples for data integration.
#'
#' @param df_par Data frame of parameters to use for data processing.
#' @param rho_adj Numeric value indicating a proportion to scale
#' rho values calculated during the ambient RNA removal step. 10% (0.1)
#' is generally suitable for most data sets.
#' @param m_cell Threshold for including features in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @param m_feat Threshold for including cells in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @param parl Logical indicating whether processing should be run in
#' parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core_perc Percentage of available cores to use if running
#' in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A processed list of sample files converted into Seurat
#' objects with a summary list of QC and processing details.
#' @examples
#'
#' # list_data <- sc_process_batch(
#' # # parameter list
#' # list_params,
#' # # adj.rho proportion
#' # 0.1,
#' # # minimum cells per feature
#' # 5,
#' # # minimum features per cell
#' # 200
#' # )
#'
#' @export
sc_process_batch <- function(
  df_par,
  rho_adj,
  m_cell,
  m_feat,
  parl,
  core_perc
) {
  if( # nolint
    Sys.info()[["sysname"]] != "Windows" && parl == TRUE
  ) {
    list_data <- setNames(
      parallel::mclapply(
        mc.cores = ceiling(
          parallel::detectCores() *
            core_perc
        ),
        seq.int(
          1,
          nrow(df_par),
          1
        ),
        function(x) {
          sc_process_file(
            # parameter list
            list_params, # nolint
            # sample ID
            x,
            # adj.rho proportion
            rho_adj,
            # minimum cells per feature
            m_cell,
            # minimum features per cell
            m_feat
          )
        }
      ),
      df_par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] != "Windows" && parl == FALSE) { # nolint
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(df_par),
        1
      ),
      function(x) {
        sc_process_file(
          # parameter list
          list_params, # nolint
          # sample ID
          x,
          # adj.rho proportion
          rho_adj,
          # minimum cells per feature
          m_cell,
          # minimum features per cell
          m_feat
        )
      }
    ),
    df_par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" && parl == TRUE) { # nolint
    print(
      "Windows OS does not support parallel sample processing, 
      defaulting to sequential processing..."
    )
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(df_par),
        1
      ),
      function(x) {
        sc_process_file(
          # parameter list
          list_params, # nolint
          # sample ID
          x,
          # adj.rho proportion
          rho_adj,
          # minimum cells per feature
          m_cell,
          # minimum features per cell
          m_feat
        )
      }
    ),
    df_par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" && parl == FALSE) { # nolint
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(df_par),
        1
      ),
      function(x) {
        sc_process_file(
          # parameter list
          list_params, # nolint
          # sample ID
          x,
          # adj.rho proportion
          rho_adj,
          # minimum cells per feature
          m_cell,
          # minimum features per cell
          m_feat
        )
      }
    ),
    df_par[["File.ID"]]
    )
  }
  return(list_data)
}

#' Plot Sample Contamination Fraction
#'
#' Generates a scatter plot indicating the individual
#' and average contamination fractions for all samples.
#'
#' @param df_par Data frame of updated parameters used for data processing.
#' @return A scatter plot indicating the individual
#' and average contamination fractions for all samples.
#' @examples
#'
#' # sc_plot_rho(list_params)
#'
#' @export
sc_plot_rho <- function(
  df_par
) {
  p_rho <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df_par,
      ggplot2::aes(
        x = .data[["File.ID"]], # nolint
        y = .data[["rho"]]
      ),
      shape = 21,
      size = 3,
      alpha = 0.25,
      fill = col_univ()[2] # nolint
    ) +
    ggplot2::geom_point(
      data = df_par,
      ggplot2::aes(
        x = .data[["File.ID"]],
        y = .data[["adj.rho"]]
      ),
      shape = 21,
      size = 3,
      alpha = 1,
      fill = col_univ()[1]
    ) +
    ggplot2::geom_hline(
      yintercept = mean(
        df_par[["rho"]]
      ),
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = mean(
        df_par[["adj.rho"]]
      ),
      linetype = "dashed",
      color = "firebrick2"
    ) +
    sc_theme1() # nolint
  ggplot2::ggsave(
    "processed/data.ambientRNA.cont.png",
    p_rho,
    width = 10,
    height = 10,
    dpi = 600
  )
}

#' Process CellRanger Multiome Data Files
#'
#' Processes a single multiome sample for data integration.
#'
#' @param df_par Data frame of parameters to use for data processing.
#' @param i Numeric value for sample ID.
#' @param p_name Character string providing the project name.
#' @param m_cell Minimum cells per feature.
#' @param m_feat Minimum features per cell.
#' @param atac_high Upper limit for per cell ATAC counts.
#' @param rna_high Upper limit for per cell RNA counts.
#' @param rna_low Lower limit for per cell RNA counts.
#' @param perc_mt Upper limit for percentage of mitochondrial reads per cell.
#' @param nsm_high Upper limit for nucleosome signal. Used to filter
#' potential doublets, potential artifacts, etc.
#' @param tss_low Lower limit for TSS enrichment score. Used to filter
#' low quality cells.
#' @return A list containing the processed multimodal Seurat object,
#' processing parameters, and QC plots.
#' @examples
#'
#' # proc_d_multiome <- sc_multiome_process(
#' #   # parameter list
#' #   l_params,
#' #   # sample ID
#' #   1,
#' #   # Project name
#' #   "Hiro.CF.multiome",
#' #   # Minimum cells per feature
#' #   5,
#' #   # Minimum features per cell
#' #   200,
#' #   # ATAC counts upper limit
#' #   100000,
#' #   # RNA counts upper limit
#' #   25000,
#' #   # RNA counts lower limit
#' #   1000,
#' #   # Percentage mitochondrial reads upper limit
#' #   20,
#' #   # Nucleosome signal upper limit
#' #   2.5,
#' #   # TSS enrichment lower limit
#' #   2
#' # )
#'
#' @export
sc_multiome_process <- function(
  df_par,
  i,
  p_name,
  m_cell,
  m_feat,
  atac_high,
  rna_high,
  rna_low,
  perc_mt,
  nsm_high,
  tss_low
) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1234)
  # Select file from chosen input parameter df
  dfp <- df_par
  df_p <- dfp[["param"]]
  d <- df_p[i, ]
  ## IMPORTANT: run "conda deactivate" in terminal then
  ## install hdf5r package if not present on machine

  # Process RNA
  ## Ambient contamination removal
  d_cnt <- Seurat::Read10X_h5(d[["Path.count"]])
  d_raw <- Seurat::Read10X_h5(d[["Path.raw"]])
  d_cnt2 <- decontX::decontX(
    SingleCellExperiment::SingleCellExperiment(
      list(
        counts = d_cnt$`Gene Expression`
      )
    ),
    background = SingleCellExperiment::SingleCellExperiment(
      list(
        counts = d_raw$`Gene Expression`
      )
    )
  )

  ## Doublet removal/Create single cell experiment
  dx_cnt <- scDblFinder::scDblFinder(d_cnt2, clusters = TRUE)
  dx_cnt <- round(
    SummarizedExperiment::assay(
      dx_cnt[, dx_cnt$scDblFinder.class == "singlet"],
      "decontXcounts"
    ),
    0
  )

  ## Create Seurat Object
  d1 <- SeuratObject::CreateSeuratObject(
    counts = dx_cnt,
    assay = "RNA",
    project = "samp1",
    min.cells = 5,
    min.features = 200
  )

  # Process ATAC
  d1_atac <- SeuratObject::CreateSeuratObject(
    counts = d_cnt$`Peaks`,
    assay = "ATAC"
  )
  d1_atac <- d1_atac[, colnames(d1_atac) %in% colnames(d1)]
  d1_atac <- Signac::CreateChromatinAssay(
    counts = SeuratObject::LayerData(
      d1_atac,
      layer = "counts"
    ),
    sep = c(":", "-"),
    fragments = d[["Path.frag"]],
    annotation = dfp[["ref.gtf"]]
  )
  d1[["ATAC"]] <- d1_atac
  remove(d1_atac)

  # Add metadata
  Seurat::DefaultAssay(d1) <- "RNA"
  d_md <- d[, 8:ncol(d)]
  d_md <- setNames(
    as.data.frame(
      lapply(
        seq.int(
          1,
          ncol(d_md),
          1
        ),
        function(x) {
          rep(
            d_md[1, x],
            nrow(d1@meta.data)
          )
        }
      )
    ),
    names(d_md)
  )
  d1 <- Seurat::AddMetaData(
    object = d1,
    metadata = d_md
  )

  ## Percentage mitochondrial reads
  d1[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object = d1,
    pattern = "^MT-"
  )

  # Add nucleosome signal and TSS enrichment values
  Seurat::DefaultAssay(d1) <- "ATAC"
  d1 <- Signac::NucleosomeSignal(d1)
  d1 <- Signac::TSSEnrichment(d1)
  # QC plot
  names(d1@meta.data)
  d1qc_pre <- Seurat::VlnPlot(
    object = d1,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "nCount_ATAC",
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    layer = "counts",
    ncol = 3,
    pt.size = 0.2
  )

  d1f <- d1[
    ,
    d1[["nCount_ATAC"]] < atac_high &
      d1[["nCount_RNA"]] < rna_high &
      d1[["percent.mt"]] < perc_mt &
      d1[["nCount_RNA"]] > rna_low &
      d1[["nucleosome_signal"]] < nsm_high &
      d1[["TSS.enrichment"]] > tss_low
  ]
  remove(d1)
  gc(reset = TRUE)
  d1qc_pst <- Seurat::VlnPlot(
    object = d1f,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "nCount_ATAC",
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    layer = "counts",
    ncol = 3,
    pt.size = 0.2
  )

  # Call ATAC peaks using default parameters
  ## NOTE: Ensure conda is activated otherwise MACS2
  ## will not be found!!
  d1pk <- Signac::CallPeaks(d1f)
  d1pk2 <- GenomeInfoDb::keepStandardChromosomes(
    d1pk,
    pruning.mode = "coarse"
  )

  d1pk_cnts <- Signac::FeatureMatrix(
    fragments = Signac::Fragments(d1f),
    features = d1pk2,
    cells = colnames(d1f)
  )

  d1f[["peaks"]] <- Signac::CreateChromatinAssay(
    counts = d1pk_cnts,
    fragments = d[["Path.frag"]],
    annotation = dfp[["ref.gtf"]]
  )

  # Normalization
  ## GEX
  Seurat::DefaultAssay(d1f) <- "RNA"
  d1f <- Seurat::SCTransform(d1f)
  d1f <- Seurat::RunPCA(d1f)

  ## ATAC Peaks
  Seurat::DefaultAssay(d1f) <- "peaks"
  d1f <- Signac::FindTopFeatures(d1f, min.cutoff = 5)
  d1f <- Signac::RunTFIDF(d1f)
  d1f <- Signac::RunSVD(d1f)

  d1f_depth <- Signac::DepthCor(d1f)

  Seurat::DefaultAssay(d1f) <- "SCT"
  d1f_sum <- head(
    Seurat::VariableFeatures(d1f),
    25
  )
  plot_var_feat <- Seurat::VariableFeaturePlot(
    d1f,
    assay = "SCT"
  )
  plot_var_feat_out <- Seurat::LabelPoints(
    plot_var_feat,
    points = d1f_sum,
    repel = TRUE
  )

  return(
    list(
      "Seurat.obj" = d1f,
      "Params" = df_p,
      "QC.pre" = d1qc_pre,
      "QC.post" = d1qc_pst,
      "var.feat.sum" = d1f_sum,
      "var.feat.plot" = plot_var_feat_out,
      "QC.depth" = d1f_depth
    )
  )
}

#' Batch Processing of CellRanger Multiome Files
#'
#' Processes a list of 10X multiome samples for downstream integration.
#'
#' @param df_par Data frame of parameters to use for data processing.
#' @param p_name Character string providing the project name.
#' @param m_cell Minimum cells per feature.
#' @param m_feat Minimum features per cell.
#' @param atac_high Upper limit for per cell ATAC counts.
#' @param rna_high Upper limit for per cell RNA counts.
#' @param rna_low Lower limit for per cell RNA counts.
#' @param perc_mt Upper limit for percentage of mitochondrial reads per cell.
#' @param nsm_high Upper limit for nucleosome signal. Used to filter
#' potential doublets, potential artifacts, etc.
#' @param tss_low Lower limit for TSS enrichment score. Used to filter
#' low quality cells.
#' @param parl Logical indicating whether processing should be run in
#' parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core_perc Percentage of available cores to use if running
#' in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A processed list of 10X multiome sample files converted
#' into Seurat objects with a summary list of QC and processing details.
#' @examples
#'
#' # list_data <- sc_process_batch_multiome(
#' #  # parameter list
#' #  l_params[["param"]],
#' #  # Project name
#' #  "Hiro.CF.multiome",
#' #  # Minimum cells per feature
#' #  5,
#' #  # Minimum features per cell
#' #  200,
#' #  # ATAC counts upper limit
#' #  100000,
#' #  # RNA counts upper limit
#' #  25000,
#' #  # RNA counts lower limit
#' #  1000,
#' #  # Percentage mitochondrial reads upper limit
#' #  20,
#' #  # Nucleosome signal upper limit
#' #  2.5,
#' #  # TSS enrichment lower limit
#' #  2,
#' #  # Run in parallel?
#' #  TRUE,
#' #  # Core percentage
#' #  0.1
#' # )
#'
#' @export
sc_process_batch_multiome <- function(
  df_par,
  p_name,
  m_cell,
  m_feat,
  atac_high,
  rna_high,
  rna_low,
  perc_mt,
  nsm_high,
  tss_low,
  parl,
  core_perc
) {
  dfp <- df_par
  if( # nolint
    Sys.info()[["sysname"]] != "Windows" && parl == TRUE
  ) {
    list_data <- setNames(
      parallel::mclapply(
        mc.cores = ceiling(
          parallel::detectCores() *
            core_perc
        ),
        seq.int(
          1,
          nrow(dfp[["param"]]),
          1
        ),
        function(x) {
          sc_multiome_process(
            # parameter list
            dfp,
            # sample ID
            x,
            # Project name
            p_name,
            # Minimum cells per feature
            m_cell,
            # Minimum features per cell
            m_feat,
            # ATAC counts upper limit
            atac_high,
            # RNA feature counts upper limit
            rna_high,
            # RNA feature counts lower limit
            rna_low,
            # Percentage mitochondrial reads
            perc_mt,
            # Nucleosome signal upper limit
            nsm_high,
            # TSS enrichment lower limit
            tss_low
          )
        }
      ),
      dfp[["param"]][["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] != "Windows" && parl == FALSE) { # nolint
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(dfp[["param"]]),
        1
      ),
      function(x) {
        sc_multiome_process(
          # parameter list
          dfp,
          # sample ID
          x,
          # Project name
          p_name,
          # Minimum cells per feature
          m_cell,
          # Minimum features per cell
          m_feat,
          # ATAC counts upper limit
          atac_high,
          # RNA feature counts upper limit
          rna_high,
          # RNA feature counts lower limit
          rna_low,
          # Percentage mitochondrial reads
          perc_mt,
          # Nucleosome signal upper limit
          nsm_high,
          # TSS enrichment lower limit
          tss_low
        )
      }
    ),
    dfp[["param"]][["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" && parl == TRUE) { # nolint
    print(
      "Windows OS does not support parallel sample processing, 
      defaulting to sequential processing..."
    )
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(dfp[["param"]]),
        1
      ),
      function(x) {
        sc_multiome_process(
          # parameter list
          dfp,
          # sample ID
          x,
          # Project name
          p_name,
          # Minimum cells per feature
          m_cell,
          # Minimum features per cell
          m_feat,
          # ATAC counts upper limit
          atac_high,
          # RNA feature counts upper limit
          rna_high,
          # RNA feature counts lower limit
          rna_low,
          # Percentage mitochondrial reads
          perc_mt,
          # Nucleosome signal upper limit
          nsm_high,
          # TSS enrichment lower limit
          tss_low
        )
      }
    ),
    dfp[["param"]][["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" && parl == FALSE) { # nolint
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(dfp[["param"]]),
        1
      ),
      function(x) {
        sc_multiome_process(
          # parameter list
          dfp,
          # sample ID
          x,
          # Project name
          p_name,
          # Minimum cells per feature
          m_cell,
          # Minimum features per cell
          m_feat,
          # ATAC counts upper limit
          atac_high,
          # RNA feature counts upper limit
          rna_high,
          # RNA feature counts lower limit
          rna_low,
          # Percentage mitochondrial reads
          perc_mt,
          # Nucleosome signal upper limit
          nsm_high,
          # TSS enrichment lower limit
          tss_low
        )
      }
    ),
    dfp[["param"]][["File.ID"]]
    )
  }
  return(list_data)
}
