#' Cell Type Prediction (HLCA/Azimuth)
#'
#' Cell type predictions of a Seurat Object using the built-in
#' Azimuth Human Lung Cell Atlas reference annotation.
#'
#' @param so A Seurat object.
#' @param ref1 Azimuth reference annotation ("lungref" by default).
#' @param f_size Future size (in MB), provided as a numeric value.
#' For larger datasets, set to 5000 or higher.
#' @param cl_var Name of the cluster variable used for cell type predictions.
#' @param md_list A vector selecting metadata columns
#' stratifying the prediction score output table.
#' @param parl Should predictions be run in parallel? (TRUE/FALSE)
#' @param core_num Number of cores to use if parl is TRUE.
#' @return A list containing a Seurat object with predicted clusters and
#' QC tables/plots to evaluate prediction performance.
#' @examples
#'
#' # pred1 <- sc_predict(
#' #   so = d,
#' #   md_list = c("Code")
#' # )
#'
#' @export
sc_predict <- function(
  so,
  ref1 = "lungref",
  f_size = 10000,
  cl_var = "seurat_clusters",
  md_list,
  parl = TRUE,
  core_num = 2
) {
  # Load data
  d1 <- so
  r1 <- ref1
  fs1 <- f_size
  # Future size and parallel settings
  if(parl == TRUE && Sys.info()[["sysname"]] != "Windows") { # nolint
    future::plan(
      "multicore",
      workers = core_num
    )
  }
  options(future.globals.maxSize = fs1 * 1024^2)
  # Prediction
  pred1 <- Azimuth::RunAzimuth(
    d1,
    reference = r1
  )
  # Reset future size and parallel settings
  if(parl == TRUE && Sys.info()[["sysname"]] != "Windows") { # nolint
    future::plan("sequential")
  }
  options(future.globals.maxSize = 500 * 1024^2)
  # Cell Type Prediction Score Distribution
  fun_dist_score <- function(x) {
    p_score_dist <- ggplot2::ggplot( # nolint
      x,
      ggplot2::aes(x = .data[["predicted.ann_finest_level.score"]]) # nolint
    ) +
      ggplot2::geom_density(
        color = "black",
        fill = col_univ()[[2]] # nolint
      ) +
      ggplot2::labs(
        x = "Prediction Score",
        y = "Density",
        title = "Prediction Score Distribution"
      ) +
      sc_theme1() # nolint
    return(p_score_dist) # nolint
  }
  p1 <- fun_dist_score(pred1@meta.data)
  print(p1)
  # Predicted Cell Type Proportions for Each Cluster
  ## Counting function
  fun_predict_prop <- function(
    x,
    c,
    cls
  ) {
    # Count cell numbers
    data_pred_prop2 <- setNames(
      dplyr::count(
        x,
        .data[[cls]],
        .data[[c]] # nolint
      ),
      c(cls, c, "Total Cells")
    )
    # Calculate proportions relative to cluster numbers
    data_pred_prop2[["Proportion"]] <- unlist(
      lapply(
        seq.int(1, nrow(data_pred_prop2), 1),
        function(i) {
          tot_clus <- sum(data_pred_prop2[
            data_pred_prop2[[cls]] ==
              data_pred_prop2[i, ][[cls]],
          ][["Total Cells"]])
          prop_type <- round(
            data_pred_prop2[i, ][["Total Cells"]] / tot_clus,
            digits = 3
          )
          return(prop_type) # nolint
        }
      )
    )
    return(data_pred_prop2) # nolint
  }
  ## Cell Type Proportion Summary and Consensus Type
  t1 <- fun_predict_prop(
    pred1@meta.data,
    "predicted.ann_finest_level",
    cl_var
  )
  ## Assign final predictions
  t1[["predicted.id"]] <- unlist(
    lapply(
      seq.int(1, nrow(t1), 1),
      function(i) {
        tot_clus <- t1[
          t1[[cl_var]] == t1[i, ][[cl_var]],
        ]
        # Assign identities
        ## Consensus ID
        if (max(tot_clus[["Proportion"]]) > 0.75) {
          pred_id <- tot_clus[
            tot_clus[["Proportion"]] == max(tot_clus[["Proportion"]]),
          ]
          pred_id <- paste(
            pred_id[[cl_var]],
            ".",
            pred_id[["predicted.ann_finest_level"]],
            "(",
            (pred_id[["Proportion"]]) * 100,
            "%)",
            sep = ""
          )
        }
        ## Mixed IDs
        if (max(tot_clus[["Proportion"]]) < 0.75) {
          pred_id <- dplyr::slice_max(tot_clus, .data[["Proportion"]], n = 3)
          pred_id <- paste(
            unique(pred_id[[cl_var]]),
            ".",
            paste(
              paste(
                pred_id[["predicted.ann_finest_level"]],
                "(", pred_id[["Proportion"]] * 100, "%)",
                sep = ""
              ),
              collapse = "_"
            ),
            sep = ""
          )
        }
        return(pred_id) # nolint
      }
    )
  )
  ### Change grouping columns to factors
  t1[["predicted.id"]] <- factor(
    t1[["predicted.id"]],
    levels = c(gtools::mixedsort(unique(t1[["predicted.id"]])))
  )
  ## Add predicted id column to Seurat object
  t2 <- dplyr::left_join(
    pred1@meta.data,
    t1[
      !duplicated(t1[[cl_var]]),
      c(cl_var, "predicted.id")
    ],
    by = cl_var
  )
  pred1 <- Seurat::AddMetaData(
    pred1,
    t2[["predicted.id"]],
    col.name = "predicted.id"
  )
  list_d <- list(
    "data" = pred1,
    "predictions" = t1
  )
  return(list_d)
}

#' Cell Type Prediction (Custom)
#'
#' Predicts cell type identities of clusters in a Seurat object
#' from a previously annotated Seurat object.
#'
#' @param so A Seurat object.
#' @param ref1 Reference Seurat object containing annotated cell types.
#' @param f_size Future size (in MB), provided as a numeric value.
#' @param cl_var Cluster variable for cell type predictions.
#' @param ct_col Column to use from the reference dataset for annotating
#' the query dataset.
#' @param asy1 Assay to use for comparing the reference and query dataset.
#' @param md_list Vector of metadata variables for stratifying proportion table.
#' @param parl Should predictions be run in parallel? (TRUE/FALSE)
#' @param core_perc Percentage of cores to use if parl is TRUE.
#' @return A list containing a Seurat object with predicted clusters and
#' QC tables/plots to evaluate prediction performance.
#' @examples
#'
#' # pred1 <- sc_predict2(
#' #   so = d,
#' #   ref1 = d_ref,
#' #   md_list = c("Code"),
#' # )
#'
#' @export
sc_predict2 <- function(
  so,
  ref1,
  f_size = 10000,
  cl_var = "seurat_clusters",
  ct_col = "CellType",
  asy1 = "RNA",
  md_list,
  parl = TRUE,
  core_perc = 0.25
) {
  # Load data
  dr <- ref1
  d1 <- so
  # Change assay and active identity
  SeuratObject::DefaultAssay(dr) <- asy1
  SeuratObject::DefaultAssay(d1) <- asy1
  d1 <- SeuratObject::SetIdent(
    d1,
    value = d1@meta.data[[cl_var]]
  )
  # Future size and parallel settings
  options(future.globals.maxSize = f_size * 1024^2)
  if(Sys.info()[["sysname"]] != "Windows" && parl == TRUE) { # nolint
    future::plan(
      "multisession",
      workers = parallel::detectCores() * core_perc
    )
  }
  ### Reference features
  data_ref <- Seurat::FindVariableFeatures(
    dr,
    selection.method = "vst",
    nfeatures = 4000,
    assay = asy1
  )
  ### Query features
  data_qry <- Seurat::FindVariableFeatures(
    d1,
    selection.method = "vst",
    nfeatures = 4000,
    assay = asy1
  )
  ## Transfer anchors
  data_transfer_anchors <- Seurat::FindTransferAnchors(
    reference = data_ref,
    query = data_qry,
    features = Seurat::VariableFeatures(object = data_ref),
    reference.assay = asy1,
    query.assay = asy1,
    reduction = "pcaproject"
  )
  ## Transfer cell type predictions to query set
  data_predicted <- Seurat::TransferData(
    anchorset = data_transfer_anchors,
    refdata = data_ref@meta.data[[ct_col]],
    weight.reduction = "pcaproject"
  )
  # reset future size and parallel settings
  if(Sys.info()[["sysname"]] != "Windows" && parl == TRUE) { # nolint
    future::plan("sequential")
  }
  options(future.globals.maxSize = 500 * 1024^2)
  ## add metadata and return the query dataset with cell type predictions
  data_qry <- Seurat::AddMetaData(
    data_qry,
    metadata = data_predicted
  )
  # Combine marker genes and Seurat object as list
  list_d <- list(
    "Predicted Clusters" = data_qry
  )
  ## Save predictions as table and Seurat object with predictions
  list_d[["Prediction Scores"]] <- dplyr::select(
    list_d$`Predicted Clusters`@meta.data,
    c(
      md_list, cl_var,
      names(
        list_d$`Predicted Clusters`@meta.data[
          grepl(
            "predicted|prediction",
            names(list_d$`Predicted Clusters`@meta.data)
          )
        ]
      )
    )
  )
  print(
    table(
      list_d$`Predicted Clusters`$prediction.score.max >
        0.5
    )
  )
  ### Cell Type Prediction Score Distribution
  fun_dist_score <- function(x) {
    p_score_dist <- ggplot2::ggplot(
      x,
      ggplot2::aes(
        x = .data[["prediction.score.max"]] # nolint
      )
    ) +
      ggplot2::geom_density(
        color = "black",
        fill = col_univ()[[2]] # nolint
      ) +
      ggplot2::labs(
        x = "Prediction Score",
        y = "Density",
        title = "Prediction Score Distribution"
      ) +
      sc_theme1() # nolint
    return(p_score_dist) # nolint
  }
  list_d[["predicted_dist"]] <- fun_dist_score(
    list_d$`Prediction Scores`
  )
  ### Predicted Cell Type Proportions for Each Cluster
  # Counting function
  fun_predict_prop_alt <- function(
    x,
    c,
    md
  ) {
    data_pred_prop2 <- setNames(
      dplyr::count(
        x,
        .data[[c]] # nolint
      ),
      c(c, "Total Cells")
    )
    data_pred_prop <- dplyr::count(
      x,
      x[, c(md, c)]
    )
    data_pred_prop <- dplyr::left_join(
      data_pred_prop,
      data_pred_prop2,
      by = c
    )
    data_pred_prop[["Proportion"]] <- round(
      data_pred_prop$n /
        data_pred_prop$`Total Cells`,
      digits = 3
    )
    return(data_pred_prop) # nolint
  }
  ## Cell Type Proportion Summary and Consensus Type
  list_d[["Predicted Proportions"]] <- dplyr::filter(
    fun_predict_prop_alt(
      list_d[["Predicted Clusters"]]@meta.data,
      "predicted.id",
      c(md_list, cl_var)
    ),
      .data[["Proportion"]] > # nolint
        0.001
  )
  list_d[["pred.prop.summary"]] <- setNames(
    aggregate(
      list_d[["Predicted Proportions"]][["Proportion"]],
      list(
        list_d[["Predicted Proportions"]][[cl_var]],
        list_d[["Predicted Proportions"]][["predicted.id"]]
      ),
      FUN = sum
    ),
    c(cl_var, "predicted.id", "Proportion")
  )
  ## Return cell type predictions and assign highest ranked type to each cell
  list_d[["cluster.proportions"]] <- dplyr::select(
    unique(
      dplyr::left_join(
        setNames(
          aggregate(
            list_d$pred.prop.summary[["Proportion"]],
            list(
              list_d$pred.prop.summary[[cl_var]]
            ),
            FUN = max
          ),
          c(cl_var,
            "Proportion"
          )
        ),
        list_d$pred.prop.summary[, c("predicted.id", "Proportion")],
        by = c("Proportion")
      )
    ),
    c(cl_var, "predicted.id", "Proportion")
  )
  list_d[["cluster.proportions"]] <- list_d[["cluster.proportions"]][
    !duplicated(list_d[["cluster.proportions"]][[cl_var]]),
  ]
  list_d[["cluster.assignments"]] <- dplyr::mutate(
    list_d$cluster.proportions,
    "predicted.id" = paste(
      seq(1:nrow(list_d$cluster.proportions)), # nolint
      gsub(
        "FOXN4", "",
        gsub(
          "\\.", "",
          gsub(
            "^*..", "",
            list_d$cluster.proportions$predicted.id
          )
        )
      ),
      sep = "."
    ),
    "CellGroup" = as.factor(
      gsub(
        "FOXN4", "",
        gsub(
          "\\.", "",
          gsub(
            "^*..", "",
            list_d$cluster.proportions$predicted.id
          )
        )
      )
    )
  )
  ### Change grouping columns to factors
  list_d[["cluster.assignments"]][["predicted.id"]] <- factor(
    list_d[["cluster.assignments"]][["predicted.id"]],
    levels = c(
      gtools::mixedsort(list_d[["cluster.assignments"]][["predicted.id"]])
    )
  )
  list_d[["cluster.assignments"]][["CellGroup"]] <- factor(
    list_d[["cluster.assignments"]][["CellGroup"]],
    levels = sort(unique(list_d[["cluster.assignments"]][["CellGroup"]]))
  )
  if("CellGroup" %in% names(list_d[["Predicted Clusters"]]@meta.data) == TRUE) { # nolint
    list_d[["cluster.assignments"]] <- dplyr::select(
      dplyr::left_join(
        list_d$`Predicted Clusters`@meta.data,
        list_d$cluster.assignments,
        by = cl_var
      ),
      c("predicted.id.y", "CellGroup.y")
    )
    list_d$`Predicted Clusters` <- Seurat::AddMetaData(
      list_d$`Predicted Clusters`,
      list_d$cluster.assignments[["CellGroup.y"]],
      col.name = "CellGroup"
    )
  }
  if("CellGroup" %in% names(list_d[["Predicted Clusters"]]@meta.data) == FALSE) { # nolint
    list_d[["cluster.assignments"]] <- dplyr::select(
      dplyr::left_join(
        list_d$`Predicted Clusters`@meta.data,
        list_d$cluster.assignments,
        by = cl_var
      ),
      c("predicted.id.y", "CellGroup")
    )
    list_d$`Predicted Clusters` <- Seurat::AddMetaData(
      list_d$`Predicted Clusters`,
      list_d$cluster.assignments[["CellGroup"]],
      col.name = "CellGroup"
    )
  }
  ### Add Cell Type and Cell Group columns to seurat object
  list_d$`Predicted Clusters` <- Seurat::AddMetaData(
    list_d$`Predicted Clusters`,
    list_d$cluster.assignments[["predicted.id.y"]],
    col.name = "CellType"
  )
  list_d <- list(
    "data" = list_d[["Predicted Clusters"]],
    "predicted_dist" = list_d[["predicted_dist"]],
    "predicted_all" = list_d[["Predicted Proportions"]],
    "predicted_sum" = list_d[["pred.prop.summary"]],
    "assigned_types" = list_d[["cluster.assignments"]]
  )
  return(list_d)
}

#' Save Predicted Data
#'
#' Saves a predicted list object in the chosen directory.
#'
#' @param so_pred A list generated by either sc_predict() or sc_predict2().
#' @param file1 File name for saving predicted results list
#' @param dir1 Directory for saving files
#' (relative to the current working directory).
#' @return Individual RDS and prediction QC files.
#' @examples
#'
#' # sc_save_pred(
#' #   so_pred = pred1,
#' #   file1 = "test",
#' #   dir1 = "analysis/"
#' # )
#'
#' @export
sc_save_pred <- function(so_pred, file1, dir1) {
  d1 <- so_pred
  # Predicted Seurat object
  saveRDS(d1, paste(dir1, file1, sep = "_"))
  # Score distribution
  ggplot2::ggsave(
    paste(dir1, file1, "distribution.png", sep = "_"),
    d1[["predicted_dist"]],
    width = 8,
    height = 8,
    dpi = 300
  )
  # Predicted types (all cells)
  write.table(
    d1[["predicted_all"]],
    paste(dir1, file1, "predicted_all.txt", sep = "_"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  # Predicted types (summary)
  write.table(
    d1[["predicted_sum"]],
    paste(dir1, file1, "predicted_sum.txt", sep = "_"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  # Predicted types (assigned)
  write.table(
    d1[["assigned_types"]],
    paste(dir1, file1, "assigned_types.txt", sep = "_"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
}
