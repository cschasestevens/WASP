#' scATAC-Seq Gene Track
#'
#' Generates a gene track from a Signac ChromatinAssay.
#' Requires scATAC-Seq peak information and a reference GRanges object
#' for plotting gene positions.
#'
#' @param so An object of class Seurat. Must contain an ATAC assay.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @param g_name1 Gene to plot, provided as a character string.
#' @param bp_window Numeric value indicating the number of base pairs to
#' extend the plotting window on each end of the selected gene's location.
#' Useful for visualizing peaks corresponding to neighboring genes.
#' @return A sequence plot including the specified gene track and all other
#' genes present within the specified window.
#' @examples
#'
#' # p_cov <- sc_seqanno_plot(
#' #   d,
#' #   rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' #   ),
#' #   "TMEM45A",
#' #   20000
#' # )
#'
#' @export
sc_seqanno_plot <- function(
  so,
  dref,
  g_name1,
  bp_window
) {
  d <- d
  ref_gene <- dref
  g_name <- "SFTPB"

  # Format reference gene annotation file
  ref_gene <- dplyr::as_tibble(ref_gene)
  # Extract detected genes from reference gene list
  g_data <- ref_gene[
    ref_gene[["gene_name"]] %in% rownames(d@assays$RNA$counts),
  ]

  # Extract individual gene location from Seurat object and map
  # against reference
  g_loc <- g_data[
    g_data$gene_name == g_name &
      g_data$type == "gene",
    c("start", "end", "width", "seqnames")
  ]

  p_pos <- g_data[
    g_data$start >= g_loc$start - 60000 &
      g_data$end <= g_loc$end + 60000,
  ]

  p_pos <- p_pos[p_pos[["seqnames"]] == as.character(g_loc[["seqnames"]]), ]
  p_pos <- p_pos[p_pos[["type"]] == "exon", ]
  p_pos <- p_pos[p_pos[["seqnames"]] == g_loc[["seqnames"]], ]
  p_peak <- p_pos[p_pos[["nearestGene"]] %in% unique(p_pos[["gene_name"]]), ]
  p_peak <- p_peak[!is.na(p_peak[["seqnames"]]), ]
  p_g_rng <- dplyr::bind_rows(
    setNames(
      lapply(
        seq.int(
          1,
          nrow(unique(p_pos[, c("gene_name", "gene_id")])),
          1
        ),
        function(x) {
          r <- unique(p_pos[, c("gene_name", "gene_id")])
          d <- p_pos[
            p_pos[["gene_name"]] == r[x, ][["gene_name"]] &
              p_pos[["gene_id"]] == r[x, ][["gene_id"]],
          ]
          d <- data.frame(
            "gene_name" = unique(d[["gene_name"]]),
            "gene_id" = unique(d[["gene_id"]]),
            "start" = min(d[["start"]]),
            "end" = max(d[["end"]])
          )
          return(d)
        }
      ),
      unique(p_pos[, c("gene_name", "gene_id")][["gene_id"]])
    )
  )
  # Plot
  p_seq <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = p_peak,
      ggplot2::aes(
        x = .data[["start"]], # nolint
        xend = .data[["end"]],
        y = 0
      ),
      linewidth = 12,
      alpha = 0.6,
      color = col_univ()[[20]] # nolint
    ) +
    ggplot2::geom_segment(
      data = p_pos[
        p_pos[["type"]] == "exon",
      ],
      ggplot2::aes(
        x = start,
        xend = end,
        color = g_name,
        y = 0
      ),
      linewidth = 8,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = p_g_rng,
      ggplot2::aes(
        x = start,
        xend = end,
        y = 0
      ),
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = p_g_rng,
      ggplot2::aes(
        x = start,
        y = -0.1,
        label = g_name,
        color = g_name
      ),
      bg.color = "white",
      show.legend = FALSE,
      size = 5
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.2, 0.2)
    ) +
    ggplot2::scale_color_manual(values = col_univ()) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(face = "bold", size = 14),
      axis.title.y = ggplot2::element_text(
        size = 14,
        angle = 90,
        color = "white"
      ),
      axis.line.x.bottom = ggplot2::element_line(color = "black"),
      axis.text.x = ggplot2::element_text(
        color = "grey40",
        size = 14,
        face = "bold"
      ),
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
      axis.ticks.length.x = grid::unit(0.2, "cm"),
      axis.ticks.x = ggplot2::element_line(color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Position (bp)", y = "Seq")
  return(p_seq)
}

#' scATAC-Seq Coverage Plot 2
#'
#' Generates a coverage plot from a Signac ChromatinAssay.
#' Requires scATAC-Seq peak information and a reference GRanges object
#' for plotting gene positions.
#'
#' @param so An object of class Seurat. Must contain an ATAC assay.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @param asy1 ATAC peaks assay to use.
#' @param g_name Gene to plot, provided as a character string.
#' @param bp_window Numeric value indicating the number of base pairs to
#' extend the plotting window on each end of the selected gene's location.
#' Useful for visualizing peaks corresponding to neighboring genes.
#' @param md_list1 Character vector of up to 2 metadata variables for
#' stratifying peak data.
#' @return A coverage plot including the specified gene track and all other
#' genes present within the specified window.
#' @examples
#'
#' # p_cov <- sc_coverage_plot(
#' #   readRDS("analysis/data.annotated.withTFs.rds"),
#' #   rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' #   ),
#' #   "TMEM45A",
#' #   20000,
#' #   c("CellType", "Airway")
#' # )
#'
#' @export
sc_coverage_plot2 <- function(
  so,
  dref,
  asy1,
  g_name,
  bp_window,
  md_list1
) {
  d <- d
  ref_gene <- dref
  md_list <- c("CellType")

  # Format reference gene annotation file
  ref_gene <- dplyr::as_tibble(ref_gene)
  # Extract detected genes from reference gene list
  g_data <- ref_gene[
    ref_gene[["gene_name"]] %in% rownames(d@assays$RNA$counts),
  ]
  head(g_data)

  # Extract individual gene location from Seurat object and map
  # against reference
  g_loc <- g_data[
    g_data$gene_name == "SFTPB" &
      g_data$type == "gene",
    c("start", "end", "width", "seqnames")
  ]
  g_loc
  p <- d@assays[["ufy.peaks"]]@meta.features
  p[["ID"]] <- seq.int(1, nrow(p), 1)
  head(p)
  pos1 <- setNames(
    as.data.frame(
      stringr::str_split_fixed(rownames(p), "-", 3)
    ),
    c("seqnames", "start", "end")
  )
  p <- cbind(p, pos1)
  g_range <- p[
    p$start >= g_loc$start - 20000 &
      p$end <= g_loc$end + 20000,
  ]
  g_range <- g_range[g_range$seqnames == as.character(g_loc$seqnames), ]
  g_range <- g_range[
    as.character(g_range[["seqnames"]]) == g_loc[["seqnames"]],
  ]

  # Extract unified peak counts from ATAC peaks assay
  if(length(md_list) == 1) { # nolint
    p2 <- setNames(data.frame(
      "col1" = d@meta.data[[md_list[[1]]]],
      setNames(
        as.data.frame(
          t(
            as.matrix(d@assays[["ufy.peaks"]]$counts[g_range$ID, ])
          )
        ),
        c(g_range$ID)
      )
    ), c(md_list[[1]], g_range[["ID"]]))
  }
  if(length(md_list) == 2) { # nolint
    p2 <- setNames(data.frame(
      "col1" = d@meta.data[[md_list[[1]]]],
      "col2" = d@meta.data[[md_list[[2]]]],
      setNames(
        as.data.frame(t(as.matrix(d@assays[[asy1]]$counts[g_range$ID, ]))),
        c(g_range$ID)
      )
    ), c(md_list[[1]], md_list[[2]], g_range[["ID"]]))
  }

  ## Calculate sum of reads per fragment per cell type
  ## Normalize reads using the following: (f.raw/n.cells)*mean(n.reads)
  if(length(md_list) == 1) { # nolint
    p2 <- dplyr::group_by(
      p2,
      .data[[md_list[[1]]]] # nolint
    )
    # raw signal
    d_raw <- data.frame(
      "col1" = levels(p2[[md_list[[1]]]]),
      as.data.frame(
        lapply(
          p2[3:ncol(p2)],
          function(y) {
            aggregate(
              y,
              list(p2[[md_list[[1]]]]),
              function(x) sum(x)
            )[[2]]
          }
        )
      )
    )
    head(d_raw)

    ## group scaling factor
    sc_norm_atac <- function(
      # raw data frame
      df_raw,
      # raw data frame (sum frequency per group)
      df_sum,
      # raw signal
      f_raw,
      # number of metadata columns
      md
    ) {
      (
        # raw signal
        f_raw /
          (
            # total cells per group
            aggregate(
              df_raw[[(md + 1)]],
              list(df_raw[[md_list[[1]]]]),
              function(x) length(x)
            )[[2]]
          )
      ) *
        (
          # average sequencing depth per group
          rowMeans(
            df_sum[2:ncol(df_sum)]
          )
        )
    }

    # normalized signal
    d_norm <- setNames(
      data.frame(
        "col1" = d_raw[["col1"]],
        as.data.frame(
          lapply(
            seq.int(2, ncol(d_raw), 1),
            function(y) {
              sc_norm_atac(
                p2,
                d_raw,
                d_raw[[y]],
                2
              )
            }
          )
        )
      ),
      c(md_list[[1]], names(p2[3:ncol(p2)]))
    )

    # format for plotting
    d_norm <- setNames(
      reshape2::melt(
        d_norm,
        id.vars = 1
      ),
      c(md_list[[1]], "ID", "freq.norm")
    )
    g_range[["ID"]] <- as.factor(g_range[["ID"]])
    g_range[["start"]] <- as.numeric(g_range[["start"]])
    g_range[["end"]] <- as.numeric(g_range[["end"]])

    d_norm <- dplyr::left_join(
      d_norm,
      g_range[, c("ID", "start", "end")],
      by = "ID"
    )

    ## 1.Split by CellType;
    ## 2.Merge each with bp interval
    ## (in steps of 500 bp [this is the fragment length]);
    ## 3.Smooth each track
    ## 4.Bind rows and plot
    d_norm <- dplyr::bind_rows(
      setNames(
        lapply(
          unique(d_norm[[md_list[[1]]]]),
          function(x) {
            d <- d_norm[d_norm[[md_list[[1]]]] == x, ]
            p2_int <- data.frame(
              "start" = seq.int(
                min(d[["start"]]) - 5000,
                max(d[["start"]] + 5000),
                by = 501
              )
            )
            d <- dplyr::full_join(
              p2_int,
              d,
              by = "start"
            )
            d <- d[order(d[["start"]]), ]
            d[is.na(d[["freq.norm"]]), "freq.norm"] <- 0
            d[is.na(d[[md_list[[1]]]]), md_list[[1]]] <- x
            ## Perform loess smoothing of tracks
            d[["freq.loess"]] <- lowess(
              x = d[["start"]],
              y = d[["freq.norm"]],
              f = 0.05,
              iter = 5,
              delta = 0
            )[[2]]
            d[d[["freq.loess"]] < 0, "freq.loess"] <- 0
            return(d)
          }
        ),
        unique(d_norm[[md_list[[1]]]])
      )
    )

    d_norm[[md_list[[1]]]] <- factor(
      d_norm[[md_list[[1]]]],
      levels = c(
        gtools::mixedsort(
          unique(d_norm[[md_list[[1]]]])
        )
      )
    )
  }

  if(length(md_list) == 2) { # nolint
    p2 <- dplyr::group_by(
      p2,
      .data[[md_list[[1]]]], # nolint
      .data[[md_list[[2]]]]
    )
    # raw signal
    d_raw <- data.frame(
      "col1" = levels(p2[[md_list[[1]]]]),
      "col2" = levels(p2[[md_list[[2]]]]),
      as.data.frame(
        lapply(
          p2[3:ncol(p2)],
          function(y) {
            aggregate(
              y,
              list(
                p2[[md_list[[1]]]],
                p2[md_list[[2]]]
              ),
              function(x) sum(x)
            )[[2]]
          }
        )
      )
    )

    ## group scaling factor
    sc_norm_atac <- function(
      # raw data frame
      df_raw,
      # raw data frame (sum frequency per group)
      df_sum,
      # raw signal
      f_raw,
      # number of metadata columns
      md
    ) {
      (
        # raw signal
        f_raw /
          (
            # total cells per group
            aggregate(
              df_raw[[(md + 1)]],
              list(
                df_raw[[md_list[[1]]]],
                df_raw[[md_list[[2]]]]
              ),
              function(x) length(x)
            )[[2]]
          )
      ) *
        (
          # average sequencing depth per group
          rowMeans(
            df_sum[3:ncol(df_sum)]
          )
        )
    }

    # normalized signal
    d_norm <- setNames(
      data.frame(
        "col1" = d_raw[["col1"]],
        "col2" = d_raw[["col2"]],
        as.data.frame(
          lapply(
            seq.int(3, ncol(d_raw), 1),
            function(y) {
              sc_norm_atac(
                p2,
                d_raw,
                d_raw[[y]],
                3
              )
            }
          )
        )
      ),
      c(md_list[[1]], md_list[[2]], names(p2[3:ncol(p2)]))
    )

    # format for plotting
    d_norm <- setNames(
      reshape2::melt(
        d_norm,
        id.vars = 1
      ),
      c(md_list[[1]], md_list[[2]], "p.ID", "freq.norm")
    )
    g_range[["ID"]] <- as.factor(g_range[["ID"]])
    g_range[["start"]] <- as.numeric(g_range[["start"]])
    g_range[["end"]] <- as.numeric(g_range[["end"]])

    d_norm <- dplyr::left_join(
      d_norm,
      g_range[, c("p.ID", "start", "end")],
      by = "p.ID"
    )

    ## 1.Split by CellType;
    ## 2.Merge each with bp interval
    ## (in steps of 500 bp [this is the fragment length]);
    ## 3.Smooth each track
    ## 4.Bind rows and plot
    d_norm <- dplyr::bind_rows(
      setNames(
        lapply(
          seq.int(1, unique(d_norm[, md_list]), 1),
          function(x) {
            d1 <- unique(d_norm[, md_list])
            d <- d_norm[
              d_norm[[md_list[[1]]]] == d1[x, md_list[[1]]] &
                d_norm[[md_list[[2]]]] == d1[x, md_list[[2]]],
            ]
            p2_int <- data.frame(
              "start" = seq.int(
                min(d[["start"]]) - 5000,
                max(d[["start"]] + 5000),
                by = 501
              )
            )
            d <- dplyr::full_join(
              p2_int,
              d,
              by = "start"
            )
            d <- d[order(d[["start"]]), ]
            d[is.na(d[["freq.norm"]]), "freq.norm"] <- 0
            d[is.na(d[[md_list[[1]]]]), md_list[[1]]] <- d1[x, md_list[[1]]]
            d[is.na(d[[md_list[[2]]]]), md_list[[2]]] <- d1[x, md_list[[2]]]

            ## Perform loess smoothing of tracks
            d[["freq.loess"]] <- lowess(
              x = d[["start"]],
              y = d[["freq.norm"]],
              f = 0.05,
              iter = 5,
              delta = 0
            )[[2]]
            d[d[["freq.loess"]] < 0, "freq.loess"] <- 0
            return(d)
          }
        ),
        paste(
          unique(d_norm[, md_list])[[1]],
          unique(d_norm[, md_list])[[2]],
          sep = "."
        )
      )
    )

    d_norm[[md_list[[1]]]] <- factor(
      d_norm[[md_list[[1]]]],
      levels = c(
        gtools::mixedsort(
          unique(d_norm[[md_list[[1]]]])
        )
      )
    )

    d_norm[[md_list[[2]]]] <- factor(
      d_norm[[md_list[[2]]]],
      levels = c(
        gtools::mixedsort(
          unique(d_norm[[md_list[[2]]]])
        )
      )
    )
  }

  # Plot
  p_cov <- ggplot2::ggplot(
    data = d_norm
  ) +
    ggplot2::geom_area(
      ggplot2::aes(
        x = .data[["start"]], # nolint
        y = .data[["freq.loess"]],
        fill = .data[[md_list[[1]]]]
      ),
      alpha = 0.7,
      show.legend = FALSE
    ) +
    sc_theme1() + # nolint
    ggplot2::scale_fill_manual(
      values = col_univ() # nolint
    ) +
    ggplot2::scale_x_continuous(
      name = "Position (bp)",
      limits = c(
        (min(d_norm[["start"]]) - 5000),
        max(d_norm[["start"]] + 5000)
      )
    ) +
    ggplot2::facet_wrap(
      ~ .data[[md_list[[1]]]],
      ncol = 1,
      strip.position = "top"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      strip.text.x.top = ggplot2::element_text(
        angle = 0,
        margin = ggplot2::margin(0, 0, 0, 0, "cm"),
        size = 10
      ),
      strip.background.x = ggplot2::element_rect(fill = "grey95"),
      panel.spacing.y = grid::unit(0, "cm"),
      axis.text.y = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(color = "black"),
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
      axis.title.x.bottom = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      y = paste(
        "LOESS-smoothed Normalized Signal",
        paste(
          "(Range ",
          round(min(d_norm[["freq.loess"]]), digits = 2),
          " - ",
          round(max(d_norm[["freq.loess"]]), digits = 2), ")", sep = ""
        ),
        sep = " "
      )
    ) +
    ggplot2::ggtitle(
      paste(
        g_loc[["seqnames"]], ":",
        min(d_norm[["start"]] - 5000),
        "-",
        max(d_norm[["start"]] + 5000),
        " (",
        "SFTPB",
        ")",
        sep = ""
      )
    )

  p_cov
  p_pos <- g_data[
    g_data$start >= g_loc$start - 20000 &
      g_data$end <= g_loc$end + 20000,
  ]

  p_pos <- p_pos[p_pos[["seqnames"]] == as.character(g_loc[["seqnames"]]), ]

  p_pos <- p_pos[p_pos[["type"]] == "exon", ]

  p_pos <- p_pos[p_pos[["seqnames"]] == g_loc[["seqnames"]], ]
  p_peak <- p[p[["nearestGene"]] %in% unique(p_pos[["gene_name"]]), ]
  p_peak <- p_peak[!is.na(p_peak[["seqnames"]]), ]
  p_g_rng <- dplyr::bind_rows(
    setNames(
      lapply(
        seq.int(
          1,
          nrow(unique(p_pos[, c("gene_name", "gene_id")])),
          1
        ),
        function(x) {
          r <- unique(p_pos[, c("gene_name", "gene_id")])
          d <- p_pos[
            p_pos[["gene_name"]] == r[x, ][["gene_name"]] &
              p_pos[["gene_id"]] == r[x, ][["gene_id"]],
          ]
          d <- data.frame(
            "gene_name" = unique(d[["gene_name"]]),
            "gene_id" = unique(d[["gene_id"]]),
            "start" = min(d[["start"]]),
            "end" = max(d[["end"]])
          )
          return(d)
        }
      ),
      unique(p_pos[, c("gene_name", "gene_id")][["gene_id"]])
    )
  )
  # Plot
  p_seq <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = p_peak,
      ggplot2::aes(
        x = .data[["start"]], # nolint
        xend = .data[["end"]],
        y = 0
      ),
      linewidth = 12,
      alpha = 0.6,
      color = col_univ()[[20]] # nolint
    ) +
    ggplot2::geom_segment(
      data = p_pos[
        p_pos[["type"]] == "exon",
      ],
      ggplot2::aes(
        x = start,
        xend = end,
        color = "SFTPB",
        y = 0
      ),
      linewidth = 8,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = p_g_rng,
      ggplot2::aes(
        x = start,
        xend = end,
        y = 0
      ),
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = p_g_rng,
      ggplot2::aes(
        x = start,
        y = -0.1,
        label = "SFTPB",
        color = "SFTPB"
      ),
      bg.color = "white",
      show.legend = FALSE,
      size = 5
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.2, 0.2)
    ) +
    ggplot2::scale_color_manual(values = col_univ()) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(face = "bold", size = 14),
      axis.title.y = ggplot2::element_text(
        size = 14,
        angle = 90,
        color = "white"
      ),
      axis.line.x.bottom = ggplot2::element_line(color = "black"),
      axis.text.x = ggplot2::element_text(
        color = "grey40",
        size = 14,
        face = "bold"
      ),
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
      axis.ticks.length.x = grid::unit(0.2, "cm"),
      axis.ticks.x = ggplot2::element_line(color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Position (bp)", y = "Seq") +
    ggplot2::scale_x_continuous(
      limits = c(
        (min(d_norm[["start"]]) - 5000),
        max(d_norm[["start"]] + 5000)
      )
    )

  p_seq

  p_cov_comb <- ggpubr::ggarrange(
    p_cov,
    p_seq,
    ncol = 1,
    nrow = 2,
    heights = c(0.9, 0.1)
  ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm")
    )

  return(p_cov_comb)
}

#' Add JASPAR TF Motifs
#'
#' Creates a ChromatinAssay from a scATAC-Seq counts assay
#' and annotates transcription factor motifs based on
#' the JASPAR database.
#'
#' @param so An object of class Seurat. Must contain an ATAC assay.
#' @param asy1 Assay to use for annotating transcription factors.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @return An annotated ChromatinAssay object.
#' @examples
#'
#' # d_chrmasy <- sc_atac_motifs(
#' #   d,
#' #   "ufy.peaks",
#' #   rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' #    )
#' # )
#'
#' @export
sc_atac_motifs <- function(
  so,
  asy1 = "ufy.peaks",
  dref
) {
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Hsapiens.UCSC.hg38)
  ## Extract counts from Seurat ATAC assay
  d <- so
  Seurat::DefaultAssay(d) <- asy1
  dc <- SeuratObject::GetAssayData(d, slot = "counts")
  ## Format gene locations and create ChromatinAssay
  ref1 <- dref
  gene_coords <- ref1[ref1$type == "gene"]
  gene_coords$gene_biotype <- gene_coords$gene_type
  GenomeInfoDb::seqlevelsStyle(gene_coords) <- "UCSC"
  gene_coords <- GenomeInfoDb::keepStandardChromosomes(
    gene_coords,
    pruning.mode = "coarse"
  )
  d1 <- Signac::CreateChromatinAssay(
    counts = dc,
    sep = c("-", "-"),
    verbose = TRUE,
    meta.data = d@meta.data
  )
  ## Use CORE, CNE, PBM, PBM_HLH, PHYLOFACTS, POLII, SPLICE, and PBM_HOMEO
  list_clt1 <- c(
    "CORE", "CNE", "PBM", "PBM_HLH",
    "PBM_HOMEO", "PHYLOFACTS", "POLII",
    "SPLICE"
  )
  list_clt <- setNames(
    lapply(
      list_clt1,
      function(x) {
        opt1 <- list(
          collection = x,
          tax_group = "vertebrates",
          all_versions = FALSE
        )
        return(opt1)
      }
    ),
    list_clt1
  )
  ## Query JASPAR and combine position frequency matrices
  list_pfm <- setNames(
    lapply(
      list_clt,
      function(y) {
        pfm1 <- TFBSTools::getMatrixSet(
          x = JASPAR2020, # nolint
          opts = y
        )
        return(pfm1)
      }
    ),
    list_clt1
  )
  list_pfm <- list_pfm[lengths(list_pfm) > 0]
  ## 4 JASPAR collections have PFMs
  list_pfm <- c(
    list_pfm[[1]],
    list_pfm[[2]],
    list_pfm[[3]],
    list_pfm[[4]]
  )
  ## Add motifs to chromatinassay
  d1 <- Signac::AddMotifs(
    object = d1,
    genome = BSgenome.Hsapiens.UCSC.hg38, # nolint
    pfm = list_pfm,
  )
  ## Save motif list
  tf_list <- Signac::GetMotifData(d1)
  tf_list <- data.frame(
    "ID" = colnames(tf_list),
    "name" = unlist(
      lapply(
        colnames(tf_list),
        function(x) {
          name(TFBSTools::getMatrixByID(JASPAR2020, ID = x)) # nolint
        }
      )
    )
  )
  write.table(
    tf_list,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    file = paste(
      "analysis/table_motifs_jaspar2020.txt"
    )
  )
  return(d1)
}

#' scATAC-Seq Coverage Plot
#'
#' Generates a coverage plot from a Signac ChromatinAssay.
#' Requires scATAC-Seq peak information and a reference GRanges object
#' for plotting gene positions. The Signac function CoveragePlot() is
#' used to generate the final plot.
#'
#' @param so Input Seurat Object.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @param asy1 ATAC peaks assay to use.
#' @param asy2 GEX assay to use.
#' @param g_name Gene or region to plot, provided as a character string.
#' @param f_name Gene to plot expression alongside gene track.
#' @param g_region (optional) Region within the coverage plot to highlight.
#' @param bp_window Vector containing the upstream and downstream limits
#' for extending the plotting window (in base pairs).
#' @param ct_col Cell type column name to use in the selected Seurat object.
#' Must be a factor variable.
#' @param ct_spl Logical indicating whether cell types should be stratified
#' by an additional variable.
#' @param ct_sel Pattern for selecting individual cell types if ct_spl == TRUE.
#' @param md_list1 Character vector of 1 metadata variable for
#' stratifying peak data.
#' @param p_type Final coverage plot type
#' (either "ct_only", "gn_only", or "gn_ct");
#' Note that "gn_only" and "gn_ct" use the variable
#' from md_list1 for plotting. This can be useful for summarizing differences
#' in accessbility between treatment groups.
#' @return A coverage plot including the specified gene track and all other
#' genes present within the specified window.
#' @examples
#'
#' # p_cov <- sc_coverage_plot2(
#' #   so = d,
#' #   dref = rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' #   ),
#' #   asy1 = "ufy.peaks",
#' #   g_name = "SCGB3A2",
#' #   bp_window = 5000,
#' #   ct_col = "CellType",
#' #   ct_spl = TRUE,
#' #   ct_sel = "Secretory",
#' #   md_list1 = c("Group"),
#' #   p_type = "gn_ct"
#' # )
#'
#' @export
sc_coverage_plot <- function( # nolint
  so,
  dref = NULL,
  asy1 = "ufy.peaks",
  asy2 = "sct",
  asy_dat = "scale.data",
  g_name,
  g_region = NULL,
  f_name,
  bp_window = c(5000, 5000),
  ct_col = NULL,
  ct_spl = TRUE,
  ct_filt = FALSE,
  ct_sel = NULL,
  md_list1 = NULL,
  p_type = "ct_only"
) {
  d <- so
  if(!is.null(g_region)) { # nolint
    g_reg <- Signac::StringToGRanges(g_region)
    g_reg$color <- "orange"
  }
  if(is.null(g_region)) { # nolint
    g_reg <- NULL
  }
  # Change assay and format reference annotation file
  Seurat::DefaultAssay(d) <- asy1
  if(!file.exists("ref/ref.gene.formatted.rds") == TRUE) { # nolint
    if(exists("dref1") == TRUE) { # nolint
      ref_gene1 <- dref1 # nolint
    }
    if(exists("dref1") == FALSE) { # nolint
      if(is.null(dref) == TRUE) { # nolint
        print("Error: No gene reference or invalid reference path provided!")
      }
      if(is.null(dref) == FALSE) { # nolint
        ref_gene1 <- dref
        ## Format gene locations and create ChromatinAssay
        ref_gene1$gene_biotype <- ref_gene1$gene_type
        ref_gene1$tx_id <- ref_gene1$transcript_id
        GenomeInfoDb::seqlevelsStyle(ref_gene1) <- "UCSC"
        ref_gene1 <- GenomeInfoDb::keepStandardChromosomes(
          ref_gene1,
          pruning.mode = "coarse"
        )
        saveRDS(ref_gene1, "ref/ref.gene.formatted.rds")
      }
    }
  }
  if(!file.exists("ref/ref.gene.formatted.rds") == FALSE) { # nolint
    if(exists("dref1") == TRUE) { # nolint
      ref_gene1 <- dref1 # nolint
    }
    if(exists("dref1") == FALSE) { # nolint
      print(
        "Using formatted reference gene annotation present in ref/ folder..."
      )
      ref_gene1 <- readRDS("ref/ref.gene.formatted.rds")
    }
  }
  # Reassign gene annotation
  Signac::Annotation(d[[asy1]]) <- ref_gene1
  # Calculate linkages for selected genes
  d <- Signac::RegionStats(
    d,
    genome = BSgenome.Hsapiens.UCSC.hg38 # nolint
  )
  d1 <- tryCatch(
    {
      d <- Signac::LinkPeaks(
        object = d,
        peak.assay = asy1,
        expression.assay = asy2,
        genes.use = c(
          g_name
        )
      )
    },
    error = function(e) {
      print("No linked peaks; selected gene is likely below limit of detection...") # nolint
    }
  )
  if(class(d1) == "SeuratObject") { # nolint
    d <- d1
    lnk <- TRUE
    remove(d1)
  }
  if(class(d1) == "character") { # nolint
    lnk <- FALSE
    remove(d1)
  }
  # Plot tracks for CellType, Airway, and split condition
  # Add column to split Cell Type and Group
  # Manually assign existing CellType column using consensus IDs
  if(ct_spl == TRUE) { # nolint
    if(p_type == "gn_ct") { # nolint
      d_spl <- data.frame(
        "CellType.split" = gtools::mixedsort(
          paste(
            unique(d@meta.data[, c(md_list1, ct_col)])[[ct_col]],
            unique(d@meta.data[, c(md_list1, ct_col)])[[md_list1]],
            sep = "."
          )
        )
      )
      d_spl[[ct_col]] <- factor(
        gsub(
          paste(
            ".", unique(d@meta.data[[md_list1]])[[1]],
            "|.", unique(d@meta.data[[md_list1]])[[2]],
            sep = ""
          ),
          "",
          d_spl[["CellType.split"]]
        ),
        levels = levels(d@meta.data[[ct_col]])
      )
      d_spl[[md_list1]] <- gsub(
        ".*\\.",
        "",
        d_spl[["CellType.split"]]
      )
      d_spl[["CellType.split"]] <- factor(
        d_spl[["CellType.split"]],
        levels = gtools::mixedsort(d_spl[["CellType.split"]])
      )
      d_ct2 <- dplyr::select(
        dplyr::left_join(
          d@meta.data[, c(ct_col, md_list1)],
          d_spl,
          by = c(ct_col, md_list1)
        ),
        .data[["CellType.split"]] # nolint
      )
      d <- SeuratObject::AddMetaData(
        d,
        d_ct2[["CellType.split"]],
        col.name = "CellType.split"
      )
      if(ct_filt == TRUE) { # nolint
        d2 <- d[, grepl(ct_sel, d@meta.data[[ct_col]])]
      }
      if(ct_filt == FALSE) { # nolint
        d2 <- d
      }
      ct_col2 <- "CellType.split"
      Seurat::Idents(d2) <- md_list1
      p_aw <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = sort(unique(d2@meta.data[[md_list1]])),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          strip.text = ggplot2::element_text(size = 6),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0, 0.25, 0),
            "cm"
          )
        )
      Seurat::Idents(d2) <- ct_col2
      p_ct <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = unique(d2@meta.data[[ct_col2]]),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0.05, 0.25, 0.25),
            "cm"
          )
        )
      p_comb <- p_aw +
        p_ct +
        plot_layout(ncol = 1, heights = c(0.75, 4)) + # nolint
        ggplot2::theme(
          plot.margin = ggplot2::unit(
            c(0.25, 0.0, 0.25, 0),
            "cm"
          )
        )
    }
    if(p_type == "gn_only") { # nolint
      d2 <- d
      Seurat::Idents(d2) <- md_list1
      p_aw <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = sort(unique(d2@meta.data[[md_list1]])),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          strip.text = ggplot2::element_text(size = 6),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0, 0.25, 0),
            "cm"
          )
        )
      p_comb <- p_aw
    }
    if(p_type == "ct_only") { # nolint
      d_spl <- data.frame(
        "CellType.split" = gtools::mixedsort(
          paste(
            unique(d@meta.data[, c(md_list1, ct_col)])[[ct_col]],
            unique(d@meta.data[, c(md_list1, ct_col)])[[md_list1]],
            sep = "."
          )
        )
      )
      d_spl[[ct_col]] <- factor(
        gsub(
          paste(
            ".", unique(d@meta.data[[md_list1]])[[1]],
            "|.", unique(d@meta.data[[md_list1]])[[2]],
            sep = ""
          ),
          "",
          d_spl[["CellType.split"]]
        ),
        levels = levels(d@meta.data[[ct_col]])
      )
      d_spl[[md_list1]] <- gsub(
        ".*\\.",
        "",
        d_spl[["CellType.split"]]
      )
      d_spl[["CellType.split"]] <- factor(
        d_spl[["CellType.split"]],
        levels = gtools::mixedsort(d_spl[["CellType.split"]])
      )
      d_ct2 <- dplyr::select(
        dplyr::left_join(
          d@meta.data[, c(ct_col, md_list1)],
          d_spl,
          by = c(ct_col, md_list1)
        ),
        .data[["CellType.split"]] # nolint
      )
      d <- SeuratObject::AddMetaData(
        d,
        d_ct2[["CellType.split"]],
        col.name = "CellType.split"
      )
      if(ct_filt == TRUE) { # nolint
        d2 <- d[, grepl(ct_sel, d@meta.data[[ct_col]])]
      }
      if(ct_filt == FALSE) { # nolint
        d2 <- d
      }
      ct_col2 <- "CellType.split"
      Seurat::Idents(d2) <- ct_col2
      p_ct <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = unique(d2@meta.data[[ct_col2]]),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0.05, 0.25, 0.25),
            "cm"
          )
        )
      p_comb <- p_ct
    }
  }
  if(ct_spl == FALSE) { # nolint
    if(p_type == "gn_ct") { # nolint
      if(ct_filt == TRUE) { # nolint
        d2 <- d[, grepl(ct_sel, d@meta.data[[ct_col]])]
      }
      if(ct_filt == FALSE) { # nolint
        d2 <- d
      }
      Seurat::Idents(d2) <- md_list1
      p_aw <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = sort(unique(d2@meta.data[[md_list1]])),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          strip.text = ggplot2::element_text(size = 6),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0, 0.25, 0),
            "cm"
          )
        )
      Seurat::Idents(d2) <- ct_col2
      p_ct <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = levels(d2@meta.data[[ct_col]]),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0.05, 0.25, 0.25),
            "cm"
          )
        )
      p_comb <- p_aw +
        p_ct +
        plot_layout(ncol = 1, heights = c(0.75, 4)) + # nolint
        ggplot2::theme(
          plot.margin = ggplot2::unit(
            c(0.25, 0.0, 0.25, 0),
            "cm"
          )
        )
    }
    if(p_type == "gn_only") { # nolint
      d2 <- d
      Seurat::Idents(d2) <- md_list1
      p_aw <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = sort(unique(d2@meta.data[[md_list1]])),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          strip.text = ggplot2::element_text(size = 6),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0, 0.25, 0),
            "cm"
          )
        )
      p_comb <- p_aw
    }
    if(p_type == "ct_only") { # nolint
      if(ct_filt == TRUE) { # nolint
        d2 <- d[, grepl(ct_sel, d@meta.data[[ct_col]])]
      }
      if(ct_filt == FALSE) { # nolint
        d2 <- d
      }
      Seurat::Idents(d2) <- ct_col
      p_ct <- Signac::CoveragePlot(
        object = d2,
        region = g_name,
        features = f_name,
        expression.assay = asy2,
        expression.slot = asy_dat,
        annotation = TRUE,
        peaks = TRUE,
        links = lnk,
        idents = levels(d2@meta.data[[ct_col]]),
        extend.upstream = bp_window[[1]],
        extend.downstream = bp_window[[2]],
        region.highlight = g_reg
      ) &
        ggplot2::scale_fill_manual(values = col_univ()) & # nolint
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 8),
          axis.title.y = ggplot2::element_text(size = 8),
          axis.text.x = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::unit(
            c(0.25, 0.05, 0.25, 0.25),
            "cm"
          )
        )
      p_comb <- p_ct
    }
  }
  return(p_comb)
}
