#' Top-Level Visualization Function
#'
#' General utility function used to toggle between plot types. All necessary
#' R objects (datasets, gene lists, etc.) must be in the global
#' environment to work properly and are defined
#'
#' @param ptype Plot type, selected from one of the following:
#' "umap_gex_list", "umap_atac_list", "umap_std", "vol_std" ...
#' @param p_params List of plot-specific parameters supplied to the function.
#' Note that function defaults are used when input parameters are missing.
#' @return Plot(s) corresponding to the chosen plot type.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom parallel mclapply
#' @import SeuratObject
#' @examples
#'
#' # p1 <- sc_visualize(
#' #   ptype = "umap_gex_list"
#' # )
#'
#' @export
sc_visualize <- function(
  ptype,
  p_params
) {
  # Input parameters specific to plot type
  par1 <- p_params
  #---- Volcano Plot ----
  if(ptype == "vol_std") { # nolint
    if(is.null(par1[["opt"]])) { # nolint
      p_vol <- sc_volcano( # nolint
        l_deg = par1[["obj_list"]],
        comp_name = par1[["title"]],
        gene_name = par1[["var_lab"]],
        diff_col = par1[["var_meas"]],
        p_col = par1[["var_sig"]],
        p_cut = 0.05,
        f_cut = 0.25,
        f_lim = par1[["limx"]],
        y_limit = par1[["limy"]],
        x_title = par1[["axis_lab"]]
      )
      ggplot2::ggsave(
        paste(
          "analysis/",
          "p_vol_",
          par1[["file"]],
          ".png",
          sep = ""
        ),
        p_vol,
        height = 10,
        width = 10,
        dpi = 700
      )
    }
    if(!is.null(par1[["opt"]])) { # nolint
      p_vol <- sc_volcano( # nolint
        l_deg = par1[["obj_list"]],
        comp_name = par1[["title"]],
        gene_name = par1[["var_lab"]],
        diff_col = par1[["var_meas"]],
        p_col = par1[["var_sig"]],
        p_cut = par1[["opt"]][[1]],
        f_cut = par1[["opt"]][[2]],
        f_lim = par1[["limx"]],
        y_limit = par1[["limy"]],
        x_title = par1[["axis_lab"]]
      )
      ggplot2::ggsave(
        paste(
          "analysis/",
          par1[["file"]],
          ".png",
          sep = ""
        ),
        p_vol,
        height = par1[["opt"]][[3]],
        width = par1[["opt"]][[4]],
        dpi = par1[["opt"]][[5]]
      )
    }
  }

  #---- Standard UMAP ----
  if(ptype == "umap_std") { # nolint
    if(is.null(par1[["opt"]])) { # nolint
      ggplot2::ggsave(
        paste(
          "analysis/UMAPs/p_umap_standard_",
          par1[["file"]],
          ".png",
          sep = ""
        ),
        sc_umap_standard( # nolint
          # Seurat object
          so = par1[["obj_list"]],
          # metadata column
          md_var = par1[["var"]],
          # reduction to plot
          slot1 = par1[["dimr"]]
        ),
        height = 8,
        width = 8,
        dpi = 900
      )
    }
    if(!is.null(par1[["opt"]])) { # nolint
      ggplot2::ggsave(
        paste(
          "analysis/UMAPs/p_umap_standard_",
          par1[["file"]],
          ".png",
          sep = ""
        ),
        sc_umap_standard( # nolint
          # Seurat object
          so = par1[["obj_list"]],
          # metadata column
          md_var = par1[["var"]],
          # reduction to plot
          slot1 = par1[["dimr"]],
          # Plot as 2D or 3D?
          dims1 = par1[["opt"]][[1]],
          # color scheme
          col1 = par1[["opt"]][[2]],
          # legend position
          pos_leg = par1[["opt"]][[3]]
        ),
        height = par1[["opt"]][[4]],
        width = par1[["opt"]][[5]],
        dpi = par1[["opt"]][[6]]
      )
    }
  }
  #---- Gene Expression UMAPs from Gene List ----
  if(ptype == "umap_gex_list") { # nolint
    out_umap_sum <- dplyr::bind_rows(parallel::mclapply(
      mc.cores = par1[["cores"]],
      seq.int(1, nrow(par1[["lg"]]), 1),
      function(x) {
        gex_umap <- dplyr::bind_rows(lapply(
          seq.int(1, length(par1[["obj_list"]]), 1),
          function(y) {
            tryCatch(
              {
                ggplot2::ggsave(
                  paste(
                    "analysis/UMAPs/p_umap_GEX",
                    "_", names(par1[["obj_list"]])[[y]],
                    "_", par1[["lg"]][x, "Gene"],
                    "_", par1[["lg"]][x, "abbv"],
                    ".png",
                    sep = ""
                  ),
                  sc_umap_panel_gene( # nolint
                    so = par1[["obj_list"]][[1]],
                    asy1 = par1[["asy"]],
                    md_var = par1[["var"]],
                    g_name = par1[["lg"]][x, "Gene"],
                    col_scheme = col_univ()[ # nolint
                      1:length( # nolint
                        levels(
                          par1[["obj_list"]][[y]]@meta.data[[par1[["var"]]]]
                        )
                      )
                    ],
                    col_names = c(
                      levels(par1[["obj_list"]][[y]]@meta.data[[par1[["var"]]]])
                    ),
                    slot1 = par1[["dimr"]],
                    col1 = par1[["col"]],
                    leg_x = par1[["leg_x"]],
                    leg_y = par1[["leg_y"]]
                  ),
                  width = 8,
                  height = 8,
                  dpi = 300
                )
              },
              error = function(e) {
                paste(par1[["lg"]][x, "Gene"], "is not present in the dataset!")
              }
            )
            SeuratObject::DefaultAssay(par1[["obj_list"]][[y]]) <- par1[["asy"]]
            chk_gene <- data.frame(
              "Gene" = par1[["lg"]][x, "Gene"],
              "Set" = names(par1[["obj_list"]])[[y]],
              "Present" = ifelse(
                par1[["lg"]][x, "Gene"] %in%
                  rownames(par1[["obj_list"]][[y]]),
                "Y",
                "N"
              ),
              "Filename" = ifelse(
                par1[["lg"]][x, "Gene"] %in%
                  rownames(par1[["obj_list"]][[y]]),
                paste(
                  "analysis/UMAPs/p_umap_GEX",
                  "_", names(par1[["obj_list"]])[[y]],
                  "_", par1[["lg"]][x, "Gene"],
                  "_", par1[["lg"]][x, "abbv"],
                  ".png",
                  sep = ""
                ),
                "NA"
              )
            )
            if( # nolint
              file.exists(
                paste(
                  "analysis/UMAPs/p_umap_GEX",
                  "_", names(par1[["obj_list"]])[[y]],
                  "_", par1[["lg"]][x, "Gene"],
                  "_", par1[["lg"]][x, "abbv"],
                  ".png",
                  sep = ""
                )
              ) &&
                chk_gene[["Present"]] == "N"
            ) {
              file.remove(
                paste(
                  "analysis/UMAPs/p_umap_GEX",
                  "_", names(par1[["obj_list"]])[[y]],
                  "_", par1[["lg"]][x, "Gene"],
                  "_", par1[["lg"]][x, "abbv"],
                  ".png",
                  sep = ""
                )
              )
            }
            return(chk_gene) # nolint
          }
        ))
        return(gex_umap) # nolint
      }
    ))
    write.table(
      out_umap_sum,
      file = "analysis/UMAPs/1_umap_GEX_output_summary.txt",
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )
  }

  #---- Accessibility UMAPs from Gene List ----
  if(ptype == "umap_atac_list") { # nolint
    out_umap_sum <- dplyr::bind_rows(parallel::mclapply(
      mc.cores = par1[["cores"]],
      seq.int(1, nrow(par1[["lg"]]), 1),
      function(x) {
        gex_umap <- dplyr::bind_rows(lapply(
          seq.int(1, length(par1[["obj_list"]]), 1),
          function(y) {
            tryCatch(
              {
                ggplot2::ggsave(
                  paste(
                    "analysis/UMAPs/p_umap_ATAC",
                    "_", names(par1[["obj_list"]])[[y]],
                    "_",
                    gsub("\\:|\\-|\\.|\\(|\\)", "", par1[["lg"]][x, "Motif"]),
                    ".png",
                    sep = ""
                  ),
                  sc_umap_panel_gene( # nolint
                    so = par1[["obj_list"]][[y]],
                    asy1 = par1[["asy"]],
                    md_var = par1[["var"]],
                    g_name = par1[["lg"]][x, "ID"],
                    col_scheme = col_univ()[ # nolint
                      1:length( # nolint
                        levels(
                          par1[["obj_list"]][[y]]@meta.data[[par1[["var"]]]]
                        )
                      )
                    ],
                    col_names = c(
                      levels(par1[["obj_list"]][[y]]@meta.data[[par1[["var"]]]])
                    ),
                    slot1 = par1[["dimr"]],
                    col1 = par1[["col"]],
                    leg_x = par1[["leg_x"]],
                    leg_y = par1[["leg_y"]]
                  ),
                  width = 8,
                  height = 8,
                  dpi = 300
                )
              },
              error = function(e) {
                paste(
                  par1[["lg"]][x, "Motif"], "is not present in the dataset!"
                )
              }
            )
            SeuratObject::DefaultAssay(par1[["obj_list"]][[y]]) <- "chromvar"
            chk_gene <- data.frame(
              "Gene" = par1[["lg"]][x, "Motif"],
              "Set" = names(par1[["obj_list"]])[[y]],
              "Present" = ifelse(
                par1[["lg"]][x, "ID"] %in%
                  rownames(par1[["obj_list"]][[y]]),
                "Y",
                "N"
              ),
              "Filename" = ifelse(
                par1[["lg"]][x, "ID"] %in%
                  rownames(par1[["obj_list"]][[y]]),
                paste(
                  "analysis/UMAPs/p_umap_ATAC",
                  "_", names(par1[["obj_list"]])[[y]],
                  "_",
                  gsub("\\:|\\-|\\.|\\(|\\)", "", par1[["lg"]][x, "Motif"]),
                  ".png",
                  sep = ""
                ),
                "NA"
              )
            )
            if( # nolint
              file.exists(
                paste(
                  "analysis/UMAPs/p_umap_ATAC",
                  "_", names(par1[["obj_list"]])[[y]],
                  "_",
                  gsub("\\:|\\-|\\.|\\(|\\)", "", par1[["lg"]][x, "Motif"]),
                  ".png",
                  sep = ""
                )
              ) &&
                chk_gene[["Present"]] == "N"
            ) {
              file.remove(
                paste(
                  "analysis/UMAPs/p_umap_ATAC",
                  "_", names(par1[["obj_list"]])[[y]],
                  "_",
                  gsub("\\:|\\-|\\.|\\(|\\)", "", par1[["lg"]][x, "Motif"]),
                  ".png",
                  sep = ""
                )
              )
            }
            return(chk_gene) # nolint
          }
        ))
        return(gex_umap) # nolint
      }
    ))
    write.table(
      out_umap_sum,
      file = "analysis/UMAPs/1_umap_ATAC_output_summary.txt",
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )
  }
}
