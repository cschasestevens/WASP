#' DGEA and DA
#'
#' Performs differential gene expression/accessibility analysis
#' per cell type from a Seurat Object.
#'
#' @param so A Seurat object
#' @param asy Character string providing the name of the assay
#' to use for differential analysis.
#' @param slot1 Character string selecting the slot from the
#' Seurat object to pull from the chosen assay.
#' @param ct CellType column name.
#' @param mast_comp Character string indicating the name of the MAST
#' group comparison for conducting DGEA. MAST names are comprised of the chosen
#' variable name and the leading factor level within that variable.
#' @param mast_name User-defined name of a DGEA comparison,
#' provided as a character string.
#' @param form1 Formula to use for MAST generalized linear model. Usually
#' consists of two terms, the first of which is the treatment group column
#' name and the second of which is the column indicating the number of features
#' present per cell.
#' @param cores Number of available cores to use if running
#' in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A list of DGEA results per cell type for the chosen group comparison,
#' including genes missing fold changes and cell type DGEA results
#' with errors.
#' @examples
#'
#' # scdgea <- sc_diff(
#' #   so = d,
#' #   mast_comp = "AirwaySAE",
#' #   mast_name = "SAE vs. LAE",
#' #   form1 = c("Group", "nFeature_SCT")
#' # )
#'
#' @export
sc_diff <- function( # nolint
  so,
  asy = "SCT",
  slot1 = "scale.data",
  ct = "CellType",
  mast_comp,
  mast_name,
  form1,
  cores = 1
) {
  # Seurat object
  d <- so
  # Cell type column
  c <- ct
  # Slot name
  slt <- slot1
  # Assay name
  assy <- asy
  # MAST formula terms
  form <- form1
  # MAST Comparison (combines column name and leading factor level for name)
  mc <- mast_comp
  # MAST Comparison name
  mn <- mast_name
  ## Input
  deg_mat <- as.matrix(
    SeuratObject::GetAssayData(d, layer = slt, assay = assy)
  )
  ## Set correct rownames if performing DA analysis
  if (assy == "chromvar") { # nolint
    rownames(deg_mat) <- rownames(d@assays[[assy]])
  }
  if (assy == "ATAC") { # nolint
    rownames(deg_mat) <- paste(
      d@assays[["chromvar"]]@meta.features$nearestGene,
      seq.int(1, nrow(d@assays[[assy]]@meta.features), 1),
      sep = "."
    )
  }
  deg_cols <- data.frame(
    d@meta.data[, c(form, c)]
  )
  ## Format input as DGEA/DA object
  if (assy == "chromvar" || assy == "ATAC") {
    dgea_sc <- MAST::FromMatrix(
      deg_mat,
      cData = deg_cols,
      check_sanity = FALSE
    )
  }
  if (assy != "chromvar" && assy != "ATAC") {
    dgea_sc <- MAST::FromMatrix(
      deg_mat,
      cData = deg_cols
    )
  }
  dgea_celltype <- unique(
    as.character(SingleCellExperiment::colData(dgea_sc)[[c]])
  )
  list_dgea <- list("SCE" = dgea_sc, "CellType" = dgea_celltype)
  remove(d)
  ## Subset data prior to differential analysis
  list_dgea_sub <- setNames(lapply(
    seq.int(1, length(list_dgea[["CellType"]]), 1),
    function(i) {
      # Subset data
      ## split cell types
      s1 <- list_dgea[[1]][ , SingleCellExperiment::colData(list_dgea[["SCE"]])[[c]] == list_dgea[["CellType"]][[i]]] #nolint
      ## filter genes not expressed in cell type
      s1_sum <- rowSums(SummarizedExperiment::assay(s1) > 0)
      ## remove NA
      s1_sum[is.na(s1_sum)] <- 0
      ## retain genes expressed in > 5% of cells per type
      s1 <- s1[s1_sum / ncol(s1) >= 0.05, ]
      return(s1) # nolint
    }
  ), list_dgea[["CellType"]])
  ## Return list of variable genes removed from each cell type
  list_miss <- dplyr::bind_rows(
    lapply(
      seq.int(1, length(list_dgea_sub), 1),
      function(i) {
        data.frame(
          "CellType" = rep(
            unique(
              as.character(
                SingleCellExperiment::colData(list_dgea_sub[[i]])[[c]]
              )
            ),
            length(
              rownames(
                list_dgea[["SCE"]]
              )[
                rownames(list_dgea[["SCE"]]) %in%
                  rownames(list_dgea_sub[[i]]) == FALSE
              ]
            )
          ),
          "GENE" = rownames(
            list_dgea[["SCE"]]
          )[
            rownames(list_dgea[["SCE"]]) %in%
              rownames(list_dgea_sub[[i]]) == FALSE
          ]
        )
      }
    )
  )
  list_miss[["Note"]] <- "Removed based on filtering"
  remove(list_dgea)
  ## Define output function
  d_mast_sum_fun <- function(
    glm_fit,
    comp1,
    ct2,
    comp1_name
  ) {
    s1_res <- MAST::summary(
      glm_fit, # nolint
      doLRT = comp1,
      logFC = TRUE,
      parallel = FALSE
    )
    ### make dfs to display summary results by comp
    s1_dt <- reshape2::melt(
      dplyr::select(
        dplyr::filter(
          s1_res$datatable,
          contrast == comp1 & # nolint
            component != "S" # nolint
        ),
        -c("contrast")
      ),
      id.vars = c("primerid", "component")
    )
    s1_dt[["vars"]] <- paste(
      s1_dt$component,
      s1_dt$variable,
      sep = "."
    )
    s1_dt <- dplyr::select(
      dplyr::mutate(
        reshape2::dcast(
          dplyr::select(
            dplyr::filter(
              s1_dt,
              vars != "logFC.Pr(>Chisq)" & # nolint
                vars != "H.ci.hi" &
                vars != "H.ci.lo" &
                vars != "H.coef" &
                vars != "H.z"
            ),
            -c(
              "component",
              "variable"
            )
          ),
          primerid ~ vars
        ),
        "CellType" = ct2,
        "Comparison" = comp1_name
      ),
      c(
        "CellType", "Comparison", "primerid",
        "logFC.coef", "H.Pr(>Chisq)", "C.Pr(>Chisq)",
        "D.Pr(>Chisq)", everything() # nolint
      )
    )
    names(s1_dt) <- c(
      "CellType", "Comparison", "GENE",
      "logFC", "H.pval", "C.pval",
      "D.pval",
      names(
        s1_dt[8:ncol(
          s1_dt
        )]
      )
    )
    return(s1_dt) # nolint
  }
  # Run DGEA/DA
  if (Sys.info()[["sysname"]] != "Windows" && cores > 1) {
    list_dgea_res <- setNames(parallel::mclapply(
      mc.cores = 4,
      seq.int(1, length(list_dgea_sub), 1),
      function(x) {
        tryCatch(
          {
            ### create glm (generalized linear model for each variable)
            s1_fit <- MAST::zlm( # nolint
              formula = as.formula(
                paste(
                  "~", form[[1]], "+",
                  form[[2]],
                  sep = " "
                )
              ),
              list_dgea_sub[[x]],
              method = "glm",
              ebayes = FALSE,
              parallel = FALSE
            )
            d1 <- d_mast_sum_fun(
              s1_fit,
              mc,
              names(list_dgea_sub)[[x]],
              mn
            )
            print(
              paste(
                "Differential analysis completed for",
                names(list_dgea_sub)[[x]]
              )
            )
            return(d1)
          },
          error = function(e) {
            print(
              paste(
                "Differential analysis unsuccessful for",
                names(list_dgea_sub)[[x]]
              )
            )
          }
        )
      }
    ), names(list_dgea_sub))
  }
  # Sequential processing
  if (Sys.info()[["sysname"]] == "Windows" || cores == 1) {
    list_dgea_res <- setNames(lapply(
      seq.int(1, length(list_dgea_sub), 1),
      function(x) {
        tryCatch(
          {
            ### create glm (generalized linear model for each variable)
            s1_fit <- MAST::zlm( # nolint
              formula = as.formula(
                paste(
                  "~", form[[1]], "+",
                  form[[2]],
                  sep = " "
                )
              ),
              list_dgea_sub[[x]],
              method = "glm",
              ebayes = FALSE,
              parallel = FALSE
            )
            d1 <- d_mast_sum_fun(
              s1_fit,
              mc,
              names(list_dgea_sub)[[x]],
              mn
            )
            print(
              paste(
                "Differential analysis completed for",
                names(list_dgea_sub)[[x]]
              )
            )
            return(d1)
          },
          error = function(e) {
            print(
              paste(
                "Differential analysis unsuccessful for",
                names(list_dgea_sub)[[x]]
              )
            )
          }
        )
      }
    ), names(list_dgea_sub))
  }
  # Combine results
  dgea_comb <- dplyr::bind_rows(list_dgea_res[lengths(list_dgea_res) > 1])
  dgea_res <- dgea_comb[!is.na(dgea_comb[["logFC"]]), ]
  dgea_res[["CellType"]] <- factor(
    dgea_res[["CellType"]],
    levels = gtools::mixedsort(
      unique(dgea_res[["CellType"]])
    )
  )
  dgea_res <- dgea_res[
    order(dgea_res[["CellType"]], dgea_res[["GENE"]]),
  ]
  ## isolate DGEA results with errors
  dgea_error <- unlist(list_dgea_res[lengths(list_dgea_res) <= 1])
  ## return genes for each result with missing logFC
  ## and combine with filtered list
  dgea_miss <- dgea_comb[
    is.na(dgea_comb[["logFC"]]),
    c("CellType", "GENE")
  ]
  dgea_miss[["Note"]] <- "Only expressed in MAST comparison group"
  dgea_miss <- dplyr::bind_rows(
    dgea_miss,
    list_miss
  )
  dgea_miss[["CellType"]] <- factor(
    dgea_miss[["CellType"]],
    levels = gtools::mixedsort(
      unique(dgea_miss[["CellType"]])
    )
  )
  dgea_miss <- dgea_miss[
    order(dgea_miss[["CellType"]], dgea_miss[["GENE"]]),
  ]
  # Output final results
  dgea_sum <- list(
    "results" = dplyr::bind_rows(dgea_res),
    "missing" = dplyr::bind_rows(dgea_miss),
    "errors" = dplyr::bind_rows(dgea_error)
  )
  dgea_sum[["results"]][["H.qval"]] <- p.adjust(
    dgea_sum[["results"]][["H.pval"]],
    method = "BH"
  )
  dgea_sum[["results"]][["log2FC"]] <- log2(
    exp(dgea_sum[["results"]][["logFC"]])
  )
  dgea_sum[["results"]] <- dplyr::select(
    dgea_sum[["results"]],
    1:4,
    H.qval, # nolint
    log2FC, # nolint
    everything() # nolint
  )
  ## Add TF names for differential accessibility analysis
  ## of chromvar transcription factor activity scores
  if (assy == "chromvar") {
    if (Sys.info()[["sysname"]] == "Windows" || cores == 1) {
      list_tf <- dplyr::bind_rows(
        lapply(
          seq.int(
            1,
            length(
              unique(
                c(dgea_sum[["results"]][["GENE"]], dgea_sum[["missing"]][["GENE"]]) # nolint
              )
            ),
            1
          ),
          function(i) {
            ltf <- unique(
              c(
                dgea_sum[["results"]][["GENE"]],
                dgea_sum[["missing"]][["GENE"]]
              )
            )
            ltf1 <- data.frame(
              "GENE" = ltf[[i]],
              "TF" = name(TFBSTools::getMatrixByID(JASPAR2020, ID = ltf[[i]])) # nolint
            )
            return(ltf1) # nolint
          }
        )
      )
    }
    if (Sys.info()[["sysname"]] != "Windows" && cores > 1) {
      list_tf <- dplyr::bind_rows(
        parallel::mclapply(
          mc.cores = 12,
          seq.int(
            1,
            length(
              unique(
                c(dgea_sum[["results"]][["GENE"]], dgea_sum[["missing"]][["GENE"]]) # nolint
              )
            ),
            1
          ),
          function(i) {
            ltf <- unique(
              c(
                dgea_sum[["results"]][["GENE"]],
                dgea_sum[["missing"]][["GENE"]]
              )
            )
            ltf1 <- data.frame(
              "GENE" = ltf[[i]],
              "TF" = name(TFBSTools::getMatrixByID(JASPAR2020, ID = ltf[[i]])) # nolint
            )
            return(ltf1) # nolint
          }
        )
      )
    }
    dgea_sum[["results"]] <- dplyr::select(
      dplyr::left_join(
        dgea_sum[["results"]],
        list_tf,
        by = "GENE"
      ),
      "CellType", "Comparison", "GENE", "TF", everything()
    )
    dgea_sum[["missing"]] <- dplyr::select(
      dplyr::left_join(
        dgea_sum[["missing"]],
        list_tf,
        by = "GENE"
      ),
      "CellType", "GENE", "TF", "Note"
    )
  }
  return(dgea_sum)
}
