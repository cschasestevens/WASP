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
#' @param var_id Sample ID column name.
#' @param mtd DGEA/DA method to use (either "MAST" or "nebula").
#' Set to nebula by default.
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
  slot1 = "counts",
  ct = "CellType",
  var_id = "Code",
  mtd = "nebula",
  form1,
  cores = 1
) {
  #---- Define parameters ----
  # method
  mtd1 <- mtd
  # Seurat object
  d <- so
  # Cell type column
  c <- ct
  # Sample ID column
  v <- var_id
  # Slot name
  slt <- slot1
  # Assay name
  assy <- asy
  # MAST formula terms
  form <- form1
  #---- Function: input ----
  fun.input <- function(dso, asy1, slt1) { # nolint
    d1 <- dso
    ## Input
    deg_mat <- as.matrix(
      SeuratObject::GetAssayData(d1, layer = slt1, assay = asy1)
    )
    ## Set correct rownames if performing DA analysis
    if (assy == "chromvar") { # nolint
      rownames(deg_mat) <- rownames(d1@assays[[asy1]])
    }
    if (assy == "ATAC") { # nolint
      rownames(deg_mat) <- paste(
        d1@assays[["chromvar"]]@meta.features$nearestGene,
        seq.int(1, nrow(d1@assays[[asy1]]@meta.features), 1),
        sep = "."
      )
    }
    # Select metadata columns
    deg_cols <- data.frame(
      d1@meta.data[, c(form, c, v)]
    )
    # Method-specific inputs
    if (mtd1 == "nebula") {
      list_dgea <- list(
        "counts" = deg_mat,
        "meta" = deg_cols,
        "formula" = model.matrix(
          as.formula(
            paste(
              "~", form[[1]], "+",
              form[[2]],
              sep = " "
            )
          ),
          data = deg_cols
        ),
        "CellType" = gtools::mixedsort(unique(as.character(deg_cols[[c]])))
      )
    }
    if (mtd1 == "MAST") {
      ## Format input as DGEA/DA object
      if (assy == "chromvar" || assy == "ATAC") {
        dgea_sc <- MAST::FromMatrix(
          deg_mat,
          cData = deg_cols,
          check_sanity = FALSE
        )
      } else {
        dgea_sc <- MAST::FromMatrix(
          deg_mat,
          cData = deg_cols
        )
      }
      dgea_celltype <- gtools::mixedsort(unique(
        as.character(SingleCellExperiment::colData(dgea_sc)[[c]])
      ))
      dgea_form <- as.formula(
        paste(
          "~", form[[1]], "+",
          form[[2]],
          sep = " "
        )
      )
      list_dgea <- list(
        "SCE" = dgea_sc,
        "CellType" = dgea_celltype,
        "Formula" = dgea_form
      )
    }
    return(list_dgea) # nolint
  }
  #---- Function: subset ----
  fun.subset <- function(ldgea) { # nolint
    ld1 <- ldgea
    list_dgea_sub <- setNames(lapply(
      seq.int(1, length(ld1[["CellType"]]), 1),
      function(i) {
        # Subset data
        if (mtd1 == "MAST") {
          ## split cell types
          s1 <- ld1[["SCE"]][ , SingleCellExperiment::colData(ld1[["SCE"]])[[c]] == ld1[["CellType"]][[i]]] #nolint
          ## filter genes not expressed in cell type
          s1_sum <- rowSums(SummarizedExperiment::assay(s1) > 0)
          ## remove NA
          s1_sum[is.na(s1_sum)] <- 0
          ## retain genes expressed in > 5% of cells per type
          s1 <- s1[s1_sum / ncol(s1) >= 0.05, ]
        }
        if (mtd1 == "nebula") {
          # Subset counts matrix
          s1 <- ld1[["counts"]][ , ld1[["meta"]][[c]] == ld1[["CellType"]][[i]]] # nolint
          # Subset id column
          s1_id <- ld1[["meta"]][ ld1[["meta"]][[c]] == ld1[["CellType"]][[i]], ] # nolint
          s1 <- list(
            "counts" = s1,
            "id" = s1_id[[v]],
            "formula" = model.matrix(
              as.formula(
                paste(
                  "~", form[[1]], "+",
                  form[[2]],
                  sep = " "
                )
              ),
              data = s1_id
            )
          )
        }
        return(s1) # nolint
      }
    ), ld1[["CellType"]])
    return(list_dgea_sub) # nolint
  }
  #---- Function: run MAST----
  fun.run.mast <- function( # nolint
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
  #---- Function: run nebula ----
  fun.run.nebula <- function(ob1, ct2) { # nolint
    dgea_neb <- nebula::nebula(
      count  = ob1[["counts"]],
      id     = ob1[["id"]],
      pred   = ob1[["formula"]],
      offset = log(colSums(ob1[["counts"]]) + 1),
      model  = "NBLMM",   # Negative Binomial Linear Mixed Model
      ncore = 1
    )
    dgea_neb[["summary"]][["FDR"]] <- p.adjust(
      dgea_neb[["summary"]][[8]], method = "fdr"
    )
    dgea_neb[["summary"]][["log2FC"]] <- log2(
      exp(dgea_neb[["summary"]][[2]])
    )
    dgea_neb[["summary"]][["CellType"]] <- ct2
    dgea_neb <- setNames(
      dplyr::select(dgea_neb[["summary"]], "CellType", "gene", "log2FC", 5, 8, "FDR"), # nolint
      c(
        "CellType",
        "GENE",
        paste(
          "log2FC:",
          unique(dgea1[["meta"]][[form[[1]]]])[[1]],
          "_",
          unique(dgea1[["meta"]][[form[[1]]]])[[2]],
          sep = ""
        ),
        "std_error",
        "pval",
        "FDR"
      )
    )
    return(dgea_neb) # nolint
  }
  #---- Function: format results ----
  fun.format <- function(list_ob1) { # nolint
    d1 <- list_ob1
    if (mtd1 == "MAST") {
      dgea_comb <- dplyr::bind_rows(d1[lengths(d1) > 1])
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
      # Output final results
      dgea_sum <- dgea_res
      dgea_sum[["H.qval"]] <- p.adjust(
        dgea_sum[["H.pval"]],
        method = "BH"
      )
      dgea_sum[["log2FC"]] <- log2(
        exp(dgea_sum[["logFC"]])
      )
      dgea_sum <- dplyr::select(
        dgea_sum,
        1:4,
        H.qval, # nolint
        log2FC, # nolint
        everything() # nolint
      )
    }
    if (mtd1 == "nebula") {
      dgea_comb <- dplyr::bind_rows(d1[lengths(d1) > 1])
      dgea_comb[["CellType"]] <- factor(
        dgea_comb[["CellType"]],
        levels = gtools::mixedsort(
          unique(dgea_comb[["CellType"]])
        )
      )
      dgea_sum <- dgea_comb[
        order(dgea_comb[["CellType"]], dgea_comb[["GENE"]]),
      ]
    }
    return(dgea_sum) # nolint
  }
  #---- Function: add transcription factor motif names ----
  fun.add.tfs <- function(ob1, cores1 = 1) { # nolint
    d1 <- dgea4
    tf.fun <- function(gname) { # nolint
      library(TFBSTools) # nolint
      library(JASPAR2020) # nolint
      ltf <- unique(d1[["GENE"]])
      ltf1 <- data.frame(
        "GENE" = ltf[[gname]],
        "TF" = name(TFBSTools::getMatrixByID(JASPAR2020, ID = ltf[[gname]])) # nolint
      )
      return(ltf1) # nolint
    }
    if (Sys.info()[["sysname"]] != "Windows" && cores1 > 1) {
      list_tf <- dplyr::bind_rows(
        parallel::mclapply(
          mc.cores = cores1,
          seq.int(1, length(unique(d1[["GENE"]])), 1),
          function(j) tf.fun(j)
        )
      )
    } else {
      list_tf <- dplyr::bind_rows(
        lapply(
          seq.int(1, length(unique(d1[["GENE"]])), 1),
          function(j) tf.fun(j)
        )
      )
    }
    d1 <- dplyr::select(
      dplyr::left_join(d1, list_tf, by = "GENE"),
      "CellType", "GENE", "TF", everything() # nolint
    )
    return(d1) # nolint
  }
  #---- Run test and output results ----
  # Input
  dgea1 <- fun.input(d, assy, slt)
  # Subsetting
  dgea2 <- fun.subset(dgea1)
  # Run test
  if (Sys.info()[["sysname"]] != "Windows" && cores > 1) {
    dgea3 <- setNames(parallel::mclapply(
      mc.cores = cores,
      seq.int(1, length(dgea2), 1),
      function(i) {
        if (mtd1 == "MAST") {
          tryCatch({
            cat("MAST started for:", names(dgea2)[[i]], "\n")
            ### create glm (generalized linear model for each variable)
            s1_fit <- MAST::zlm( # nolint
              formula = dgea1[["Formula"]],
              dgea2[[i]],
              method = "glm",
              ebayes = FALSE,
              parallel = FALSE
            )
            d1 <- fun.run.mast(
              s1_fit,
              paste(form[[1]], unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[1]], sep = ""), # nolint
              names(dgea2)[[i]],
              paste(
                unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[1]], # nolint
                "vs.",
                unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[2]] # nolint
              )
            )
            cat("MAST successful for:", names(dgea2)[[i]], "\n")
          }, error = function(e) {
            cat("MAST failed for:", names(dgea2)[[i]], "\n")
          })
        }
        if (mtd1 == "nebula") {
          tryCatch({
            cat("NEBULA started for:", names(dgea2)[[i]], "\n")
            d1 <- fun.run.nebula(dgea2[[i]], names(dgea2)[[i]])
            cat("NEBULA successful for:", names(dgea2)[[i]], "\n")
          }, error = function(e) {
            cat("NEBULA failed for:", names(dgea2)[[i]], "\n")
          })
        }
        return(d1) # nolint
      }
    ), names(dgea2))
  } else {
    dgea3 <- setNames(lapply(
      seq.int(1, length(dgea2), 1),
      function(i) {
        if (mtd1 == "MAST") {
          tryCatch({
            cat("MAST started for:", names(dgea2)[[i]], "\n")
            ### create glm (generalized linear model for each variable)
            s1_fit <- MAST::zlm( # nolint
              formula = dgea1[["Formula"]],
              dgea2[[i]],
              method = "glm",
              ebayes = FALSE,
              parallel = FALSE
            )
            d1 <- fun.run.mast(
              s1_fit,
              paste(form[[1]], unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[1]], sep = ""), # nolint
              names(dgea2)[[i]],
              paste(
                unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[1]], # nolint
                "vs.",
                unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[2]] # nolint
              )
            )
            cat("MAST successful for:", names(dgea2)[[i]], "\n")
          }, error = function(e) {
            cat("MAST failed for:", names(dgea2)[[i]], "\n")
          })
        }
        if (mtd1 == "nebula") {
          tryCatch({
            cat("NEBULA started for:", names(dgea2)[[i]], "\n")
            d1 <- fun.run.nebula(dgea2[[i]], names(dgea2)[[i]])
            cat("NEBULA successful for:", names(dgea2)[[i]], "\n")
          }, error = function(e) {
            cat("NEBULA failed for:", names(dgea2)[[i]], "\n")
          })
        }
        return(d1) # nolint
      }
    ), names(dgea2))
  }
  # Combine results for each cell type and save
  dgea4 <- fun.format(dgea3)
  if (assy == "chromvar") {
    dgea4 <- fun.add.tfs(dgea4, cores1 = 12) # nolint
  }
  return(dgea4) # nolint
}
