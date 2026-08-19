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
  #---- Function: linear mixed effects model (chromvar assays only) ----
  fun.lme <- function(dso) { # nolint
    d1 <- dso
    ## Input
    deg_mat <- as.matrix(
      SeuratObject::GetAssayData(d1, layer = "data", assay = "chromvar")
    )
    deg_mat <- as.data.frame(t(deg_mat))
    # Select metadata columns
    deg_cols <- data.frame(
      d1@meta.data[, c(form[[1]], c, v)]
    )
    deg_cols[[form[[1]]]] <- factor(
      deg_cols[[form[[1]]]],
      levels = gtools::mixedsort(
        unique(deg_cols[[form[[1]]]])
      )
    )
    rownames(deg_cols) <- rownames(d1@meta.data)
    cat("Checking for identical rownames in counts and metadata", "\n")
    stopifnot(identical(rownames(deg_mat), rownames(deg_cols)))
    list_dgea <- cbind(deg_mat, deg_cols)
    head(list_dgea)
    # Subset data by cell type
    list_dgea_sub <- setNames(lapply(
      seq.int(1, length(levels(list_dgea[[c]])), 1),
      function(i) {
        # Subset into chromvar scores and metadata
        s1 <- list_dgea[list_dgea[[c]] == list_dgea[[c]][[i]], ] # nolint
        return(s1) # nolint
      }
    ), levels(list_dgea[[c]]))
    # Fit lme to each motif per cell type
    list_motif <- names(deg_mat)
    fit.lme <- function(motif, data1) { # nolint
      # define formula
      def_form <- as.formula(
        paste0(
          "`", motif, "` ~ ",
          form[[1]], " + ",
          "(1 | ", v, ")"
        )
      )
      # Run model
      tryCatch({
        fit <- lmerTest::lmer(def_form, data = list_dgea_sub[[1]])
        coefs <- summary(fit)$coefficients
        # Extract the row(s) corresponding to your condition variable
        # This assumes a two-level factor; adjust if you have >2 levels
        coef_row <- coefs[
          grep(paste0("^", form[[1]]), rownames(coefs)),
          ,
          drop = FALSE
        ]
        dout <- data.frame(
          motif      = motif,
          contrast   = rownames(coef_row),
          estimate   = coef_row[, "Estimate"],
          std_error  = coef_row[, "Std. Error"],
          df = coef_row[, "df"],
          t_value    = coef_row[, "t value"],
          p_value    = coef_row[, "Pr(>|t|)"],
          row.names  = NULL
        )
      }, error = function(e) {
        # Return NAs for motifs where the model fails to converge
        message("Model failed for motif: ", motif, " — ", e$message)
      })
      return(dout) # nolint
    }
    # Fit linear model for each celltype
    list_lme <- dplyr::bind_rows(lapply(
      seq.int(1, length(list_dgea_sub), 1),
      function(i) {
        cat("Fitting lme for:", names(list_dgea_sub)[[i]], "\n")
        if (Sys.info()[["sysname"]] != "Windows" && cores > 1) {
          fit_lme <- parallel::mclapply(
            mc.cores = cores,
            seq.int(1, length(list_motif), 1),
            function(j) {
              fit.lme(list_motif[[j]], list_dgea_sub[[i]])
            }
          )
        } else {
          fit_lme <- lapply(
            seq.int(1, length(list_motif), 1),
            function(j) {
              fit.lme(list_motif[[j]], list_dgea_sub[[i]])
            }
          )
        }
        fit_lme <- dplyr::bind_rows(fit_lme)
        # calculate FDR adjusted p-values
        fit_lme[["FDR"]] <- p.adjust(fit_lme[["p_value"]], method = "BH")
        # assign cell type column
        fit_lme[["CellType"]] <- names(list_dgea_sub)[[i]]
        # format and output
        fit_lme <- dplyr::select(fit_lme, "CellType", everything())
        fit_lme <- fit_lme[order(fit_lme[["p_value"]], na.last = TRUE), ]
        return(fit_lme) # nolint
      }
    ))
    return(list_lme) # nolint
  }
  #---- Function: limma and empirical Bayes for chromVAR data ----
  fun.limma <- function(dso) { # nolint
    d1 <- dso
    ## Input
    deg_mat <- SeuratObject::GetAssayData(
      d1, layer = "data", assay = "chromvar"
    )
    # Select metadata columns
    deg_cols <- data.frame(
      d1@meta.data[, c(form[[1]], c, v)]
    )
    deg_cols[[form[[1]]]] <- factor(
      deg_cols[[form[[1]]]],
      levels = gtools::mixedsort(
        unique(deg_cols[[form[[1]]]])
      )
    )
    rownames(deg_cols) <- rownames(d1@meta.data)
    # compute mean z-scores for each motif per id, celltype, and group
    deg_cols2 <- dplyr::distinct(deg_cols)
    # subset and calculate rowMeans, then combine
    deg_mat <- dplyr::bind_rows(
      lapply(
        seq.int(1, nrow(deg_cols2), 1),
        function(j) {
          df3 <- rownames(deg_cols)[
            deg_cols[[form[[1]]]] == deg_cols2[j, ][[form[[1]]]] &
              deg_cols[[c]] == deg_cols2[j, ][[c]] &
              deg_cols[[v]] == deg_cols2[j, ][[v]]
          ]
          df3 <- setNames(
            as.data.frame(rowMeans(deg_mat[, df3, drop = FALSE])),
            "mean_zscore"
          )
          df3[["motif"]] <- rownames(df3)
          df3[["CellType"]] <- deg_cols2[j, ][[c]]
          df3[["Code"]] <- deg_cols2[j, ][[v]]
          df3[["Group"]] <- deg_cols2[j, ][[form[[1]]]]
          df3 <- dplyr::select(
            df3, "Group", "Code", "CellType", "motif", "mean_zscore"
          )
          return(df3) # nolint
        }
      )
    )
    # fit limma for each celltype
    ## define function
    fun.run.lim <- function(ob1) { # nolint
      d1 <- ob1
      d1 <- dplyr::select(
        magrittr::set_rownames(d1, d1[["motif"]]),
        -c("motif")
      )
      d1grp <- data.frame(
        "colname" = paste(
          dplyr::distinct(deg_mat, Group, Code)[["Code"]], # nolint
          dplyr::distinct(deg_mat, Group, Code)[["Group"]],
          sep = "_"
        ),
        "Group" = dplyr::distinct(deg_mat, Group, Code)[["Group"]]
      )
      d1grp <- d1grp[match(colnames(d1), d1grp[["colname"]]), ]
      lim_grp <- factor( # nolint
        d1grp[["Group"]],
        levels = unique(d1grp[["Group"]])
      )
      lim_dsg <- model.matrix(~ lim_grp)
      # fit model
      lim_fit <- limma::lmFit(d1, lim_dsg)
      # empirical bayes shrinkage for small sample size
      lim_fit <- limma::eBayes(lim_fit)
      # Return results df
      lim_res <- limma::topTable(
        lim_fit,
        coef = paste0("lim_grp", unique(d1grp[["Group"]])[[2]]),
        number = Inf,
        sort.by = "p"
      )
      lim_res[["Comparison"]] <- paste(
        unique(d1grp[["Group"]])[[2]],
        unique(d1grp[["Group"]])[[1]],
        sep = " vs. "
      )
      lim_res[["log2FC"]] <- log2(exp(lim_res[["logFC"]]))
      return(lim_res) # nolint
    }
    ## run for each cell type
    deg_fit <- dplyr::bind_rows(
      lapply(
        seq.int(1, length(unique(deg_mat[["CellType"]])), 1),
        function(k) {
          cat("Running limma for:", levels(deg_mat[["CellType"]])[[k]], "\n")
          # subset data and format
          d1 <- deg_mat[
            deg_mat[["CellType"]] == unique(deg_mat[["CellType"]])[[k]],
          ]
          d1 <- dplyr::select(d1, -c("CellType"))
          # cast by donor and group
          d2 <- reshape2::dcast(d1, motif ~ Code + Group, value.var = "mean_zscore") # nolint
          # run limma fit
          d3 <- fun.run.lim(d2)
          d3[["CellType"]] <- levels(deg_mat[["CellType"]])[[k]]
          d3[["motif"]] <- rownames(d3)
          d3 <- dplyr::select(
            d3, "Comparison", "CellType", "motif", "log2FC", everything()
          )
          return(d3) # nolint
        }
      )
    )
    return(deg_fit) # nolint
  }
  #---- Function: input ----
  fun.input <- function(dso, asy1, slt1) { # nolint
    d1 <- dso
    ## Input
    deg_mat <- as.matrix(
      SeuratObject::GetAssayData(d1, layer = slt1, assay = asy1)
    )
    ## Set correct rownames if performing DA analysis
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
      if (assy == "ATAC") {
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
          ## filter low expressed genes or TFs
          cat("removing the following genes/TFs with missing counts in cell type:", "\n") # nolint
          ### replace missing counts with 0
          s1[["counts"]][is.na(s1[["counts"]])] <- 0
          cat(rownames(s1[["counts"]])[rowSums(s1[["counts"]]) == 0], "\n") # nolint
          s1[["counts"]] <- s1[["counts"]][rowSums(s1[["counts"]]) > 0, ]
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
  fun.add.tfs <- function(ob1, cores1 = 1, motif_col = "motif") { # nolint
    d1 <- ob1
    if (motif_col != "motif") {
      d1[["motif"]] <- d1[[motif_col]]
      d1 <- dplyr::select(d1, -c(motif_col))
    }
    tf.fun <- function(gname) { # nolint
      library(TFBSTools) # nolint
      library(JASPAR2020) # nolint
      ltf <- unique(d1[["motif"]])
      ltf1 <- data.frame(
        "motif" = ltf[[gname]],
        "TF" = name(TFBSTools::getMatrixByID(JASPAR2020, ID = ltf[[gname]])) # nolint
      )
      return(ltf1) # nolint
    }
    if (Sys.info()[["sysname"]] != "Windows" && cores1 > 1) {
      list_tf <- dplyr::bind_rows(
        parallel::mclapply(
          mc.cores = cores1,
          seq.int(1, length(unique(d1[["motif"]])), 1),
          function(j) tf.fun(j)
        )
      )
    } else {
      list_tf <- dplyr::bind_rows(
        lapply(
          seq.int(1, length(unique(d1[["motif"]])), 1),
          function(j) tf.fun(j)
        )
      )
    }
    d1 <- dplyr::select(
      dplyr::left_join(d1, list_tf, by = "motif"),
      "CellType", "motif", "TF", everything() # nolint
    )
    return(d1) # nolint
  }
  #---- Run test and output results ----
  if (mtd1 == "MAST" || mtd1 == "nebula") {
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
            cat("Selected method is MAST; running MAST differential analysis...", "\n") # nolint
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
                paste(form[[1]], sort(unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]]))[[2]], sep = ""), # nolint
                names(dgea2)[[i]],
                paste(
                  unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[2]], # nolint
                  "vs.",
                  unique(SingleCellExperiment::colData(dgea1[[1]])[[form[[1]]]])[[1]] # nolint
                )
              )
              cat("MAST successful for:", names(dgea2)[[i]], "\n")
            }, error = function(e) {
              cat("MAST failed for:", names(dgea2)[[i]], "\n")
            })
          }
          if (mtd1 == "nebula") {
            cat("Selected method is nebula; running nebula differential analysis...", "\n") # nolint
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
            cat("Selected method is MAST; running MAST differential analysis...", "\n") # nolint
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
            cat("Selected method is nebula; running nebula differential analysis...", "\n") # nolint
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
      dgea4 <- fun.add.tfs(dgea4, cores1 = 12, motif_col = "GENE")
    }
  }
  if (mtd1 == "limma" && assy == "chromvar") {
    cat("Selected assay is `chromvar`; running limma model...", "\n")
    # Fit linear model
    dgea4 <- fun.limma(d)
    # Assign motif names
    dgea4 <- fun.add.tfs(dgea4, cores1 = 12) # nolint
  }
  if (mtd1 == "lme" && assy == "chromvar") {
    cat("Selected assay is `chromvar`; running linear model...", "\n")
    # Fit linear model
    dgea4 <- fun.lme(d)
    # Assign motif names
    dgea4 <- fun.add.tfs(dgea4, cores1 = 12) # nolint
  }
  return(dgea4) # nolint
}
