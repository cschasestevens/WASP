#' Reclustering Analysis
#'
#' Performs reclustering and basic visualization
#' of a subsetted Seurat object subset.
#' Subsets are determined based on cell groups.
#'
#' @param so A Seurat object.
#' @param rc_type Data set type being reclustered (either "GEX" or "Mult").
#' @param ct_col Column containing cell type information.
#' @param slot1 Assay from Seurat object subset to use (ex. "RNA").
#' @param slot2 Secondary assay to use for reclustering
#' (only if rc_type is "Mult").
#' @param red Original dimension reduction to use for reclustering.
#' @param sid Sample ID column.
#' @param md_list A vector of character strings indicating metadata
#' columns for overlaying on a loadings plot.
#' @param dim1 Number of dimensions to use in PCA.
#' @param batch_cor Should Harmony based batch correction be performed?
#' This parameter is automatically set to TRUE for multiome datasets.
#' @param cor_col Variables to account for in batch correction.
#' @param res1 Clustering resolution to use for final cluster determination.
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param cl_c Cluster columns?
#' @param cl_r Cluster rows?
#' @param rot_c Rotation of column names.
#' @param col1 Gradient color scheme to use
#' (must be exactly 4 colors in length).
#' @param fsize Future size for cell type predictions (in GB).
#' @param mdir Directory for saving marker gene table.
#' @return A reclustered Seurat Object with summary plots, proportion tables,
#' and marker gene lists.
#' @examples
#'
#' # d_recluster <- sc_recluster(
#' #   so = d_gene,
#' #   md_list = c("Code", "Airway", "Batch"),
#' #   cor_col = c("Batch", "Code")
#' # )
#'
#' @export
sc_recluster <- function(
  so,
  rc_type = "GEX",
  ct_col = "CellType",
  slot1 = "RNA",
  slot2 = "ufy.peaks.corrected",
  red = "wnn.umap",
  sid = "Code",
  md_list,
  dim1 = 50,
  batch_cor = TRUE,
  cor_col,
  res1 = 0.5,
  h_w = 38,
  h_h = 18,
  fs_c = 6,
  fs_r = 9,
  cl_c = FALSE,
  cl_r = FALSE,
  rot_c = 45,
  col1 = col_grad(scm = 3),
  mdir = "analyze/",
  fsize = 45000
) {
  # Load objects
  d <- so
  # Reclustering of traditional scRNA-Seq data
  if(rc_type == "GEX") { # nolint
    print("Performing RNA reclustering...")
    SeuratObject::DefaultAssay(d) <- slot1
    # Re-run dimension reduction on object subset
    options(future.globals.maxSize = fsize * 1024 ^ 2)
    print(paste("---- Step 1: Normalize subsetted data ----"))
    d <- Seurat::SCTransform(d)
    # Run PCA
    d <- Seurat::RunPCA(d)
    ## Run Harmony if batch effect exists
    print(paste("---- Step 2: Batch correction ----"))
    if(batch_cor == TRUE) { # nolint
      d <- harmony::RunHarmony(
        d,
        assay.use = "SCT",
        group.by.vars = cor_col,
        reduction.use = "pca",
        reduction.save = "pca.cor",
        project.dim = FALSE
      )
      # Run UMAP
      Seurat::Idents(d) <- ct_col
      d <- Seurat::RunUMAP(
        object = d,
        reduction = "pca.cor",
        dims = 1:dim1,
        n.components = 3
      )
    }
    if(batch_cor == FALSE) { # nolint
      # Run UMAP
      Seurat::Idents(d) <- ct_col
      d <- Seurat::RunUMAP(
        object = d,
        reduction = "pca",
        dims = 1:dim1,
        n.components = 3
      )
    }
    SeuratObject::DefaultAssay(d) <- "SCT"
    print(paste("---- Step 3: Recluster ----"))
    d <- Seurat::FindNeighbors(
      d,
      reduction = "umap",
      dims = 1:3,
      verbose = TRUE
    )
    d <- Seurat::FindClusters(
      d,
      cluster.name = "recluster",
      resolution = res1
    )
    ## Reclustering Performance
    print(paste("---- Step 4: QC and marker gene calculation ----"))
    d1qc_pre <- sc_violin(
      so = d,
      ct_col = "recluster",
      gene_col = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt"
      )
    )
    # Marker gene heatmap
    d_hmap <- sc_heatmap( # nolint
      so = d,
      asy = "SCT",
      cl_var = "recluster",
      h_w = h_w,
      h_h = h_h,
      fs_c = fs_c,
      fs_r = fs_r,
      cl_c = cl_c,
      cl_r = cl_r,
      rot_c = rot_c,
      col1 = col1
    )
    # Cell type predictions
    print(paste("---- Step 5: Predict cell types ----"))
    Seurat::DefaultAssay(d) <- "RNA"
    d <- sc_predict(
      so = d,
      md_list = c(sid),
      parl = FALSE,
      f_size = fsize,
      cl_var = "recluster"
    )
    ## Visualize Clusters
    d_umap1 <- sc_umap_panel( # nolint
      d[["data"]],
      c(
        ct_col, md_list, "recluster", "predicted.id"
      ),
      "re_wnn_umap"
    )
  }
  if(rc_type == "Mult") { # nolint
    print("Performing multiome reclustering...")
    # Re-run SCTransform and batch correct GEX data;
    # leave corrected ATAC data as-is
    # Normalization
    ## GEX
    Seurat::DefaultAssay(d) <- slot1
    options(future.globals.maxSize = fsize * 1024 ^ 2)
    print(paste("---- Step 1: Normalize subsetted data ----"))
    d <- Seurat::FindVariableFeatures(
      d,
      selection.method = "vst",
      nfeatures = 3000
    )
    d <- Seurat::SCTransform(d)
    d <- Seurat::RunPCA(d)
    names(d@meta.data)
    print(paste("---- Step 2: Batch correction ----"))
    d <- harmony::RunHarmony(
      d,
      assay.use = "SCT",
      group.by.vars = c(sid),
      reduction.use = "pca",
      reduction.save = "pca_cor",
      project.dim = FALSE
    )
    ## WNN
    print(paste("---- Step 3: WNN ----"))
    d <- Seurat::FindMultiModalNeighbors(
      d,
      reduction.list = list(
        "pca_cor",
        slot2
      ),
      dims.list = list(1:50, 2:50)
    )
    d <- Seurat::RunUMAP(
      d,
      nn.name = "weighted.nn",
      reduction.name = "re_wnn_umap",
      reduction.key = "wnnUMAP_",
      n.components = 3
    )
    d <- Seurat::FindNeighbors(
      d,
      reduction = "re_wnn_umap",
      dims = 1:3,
      verbose = TRUE
    )
    print(paste("---- Step 4: Recluster ----"))
    d <- Seurat::FindClusters(
      d,
      graph.name = "wsnn",
      cluster.name = "recluster",
      algorithm = 3,
      verbose = TRUE,
      resolution = res1
    )
    d@meta.data[["recluster"]] <- factor(
      as.numeric(d@meta.data[["recluster"]]),
      levels = seq.int(
        1,
        length(unique(as.numeric(d@meta.data[["recluster"]]))),
        1
      )
    )
    ## Reclustering Performance
    print(paste("---- Step 5: QC and marker gene calculation ----"))
    d1qc_pre <- sc_violin(
      so = d,
      ct_col = "recluster",
      gene_col = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt",
        "nCount_ATAC",
        "TSS.enrichment",
        "nucleosome_signal"
      )
    )
    # Marker gene heatmap
    d_hmap <- sc_heatmap( # nolint
      so = d,
      asy = "SCT",
      cl_var = "recluster",
      h_w = h_w,
      h_h = h_h,
      fs_c = fs_c,
      fs_r = fs_r,
      cl_c = cl_c,
      cl_r = cl_r,
      rot_c = rot_c,
      col1 = col1,
      mark_dir = mdir
    )
    # Cell type predictions
    print(paste("---- Step 6: Cell type predictions ----"))
    Seurat::DefaultAssay(d) <- "RNA"
    d <- sc_predict(
      so = d,
      md_list = c(sid),
      parl = FALSE,
      f_size = fsize,
      cl_var = "recluster"
    )
    ## Visualize Clusters
    d_umap1 <- sc_umap_panel( # nolint
      d[["data"]],
      c(
        ct_col, md_list, "recluster", "predicted.id"
      ),
      "re_wnn_umap"
    )
  }
  print("Reclustering completed!")
  return(
    list = c(
      "data" = d,
      "umap_panel" = d_umap1,
      "hmap_top10" = d_hmap,
      "qc_vio" = d1qc_pre
    )
  )
}
