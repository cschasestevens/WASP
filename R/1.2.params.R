#' Define Processing Parameters
#'
#' Creates a data frame of processing parameters
#' used in scRNA-Seq data processing by Seurat.
#'
#' @param d_path Character string indicating the path
#' to a set of CellRanger files for data processing.
#' @param study_md A data frame containing sample
#' metadata variables specific to an individual study.
#' @return Data frame containing a list of parameters
#' to use for scRNA-Seq data processing by Seurat.
#' @examples
#'
#' # # Dataset input parameters
#' # list_params <- create_proc_params(
#' #   "data/",
#' #   data.frame(
#' #     # ID column name
#' #     # (splits by underscore, should be listed
#' #     # first in folder name as in 's01_1_KO_a')
#' #     "Code" = unlist(
#' #        lapply(
#' #          strsplit(
#' #            basename(list.files("data/")),
#' #            "_",
#' #            fixed = T
#' #          ),
#' #          "[",
#' #          1
#' #        )
#' #      ),
#' #     # 1st metadata column (include in CellRanger folder name)
#' #     "Knockout" = as.factor(
#' #                    ifelse(
#' #                      grepl("NG",basename(list.files("data/"))),
#' #                      "ctrl",
#' #                      "KO"
#' #                    )
#' #                  ),
#' #     # 2nd metadata column
#' #     "Region" = as.factor(
#' #                  ifelse(
#' #                    grepl("LAE",basename(list.files("data/"))),
#' #                    "1",
#' #                    "2"
#' #                  )
#' #                ),
#' #     # 3rd metadata column (add/remove columns as needed)
#' #     "Time" = as.factor(
#' #                ifelse(
#' #                  grepl("D28",basename(list.files("data/"))),
#' #                  "a",
#' #                  "b"
#' #                )
#' #              )
#' #   )
#' # )
#'
#' @export
create_proc_params <- function(d_path, study_md) {

  list_params <- data.frame(
    # Universal columns
    data.frame(
      # Sample Number
      Sample.No = seq(1:length( # nolint
        basename(
          list.files(d_path)
        )
      )
      ),
      # Individual file names (uses CellRanger folder name by default)
      File.ID = basename(list.files(d_path)),
      # Data file paths (Location of CellRanger files: 'Data/' by default)
      Path = paste(d_path, list.files(d_path), sep = ""),
      # Path to feature files
      Path.feat = paste(
        d_path,
        list.files(d_path),
        "/filtered_feature_bc_matrix/features.tsv.gz",
        sep = ""
      )
    ),
    # Dataset-specific columns
    study_md
  )
  return(list_params) # nolint
}

#' Define Processing Parameters
#'
#' Creates a data frame of processing parameters
#' used in scRNA-Seq data processing by Seurat.
#'
#' @param study_md A data frame containing sample
#' metadata variables specific to an individual study.
#' @param gtf_path Character string providing the path
#' to a GENCODE gene annotation file.
#' @param fa_path Character string providing the path
#' to a GENCODE FASTA genome file.
#' @return A list of parameters to use for processing
#' 10X multiome files.
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle
#' @importFrom rtracklayer import
#' @importFrom Rsamtools FaFile
#' @examples
#'
#' # # Dataset input parameters
#' # l_params <- sc_multiome_params(
#' #   data.frame(
#' #     # ID column name
#' #     "Code" = basename(list.files("data/")),
#' #     # 1st metadata column (include in CellRanger folder name)
#' #     "Group" = as.factor(
#' #       ifelse(
#' #         grepl("DD", basename(list.files("data/"))),
#' #         "NonCF",
#' #         "CF"
#' #       )
#' #     ),
#' #   ),
#' #   "ref/gencode.v45.primary_assembly.annotation.gtf",
#' #   "ref/GRCh38.primary_assembly.genome.fa"
#' # )
#'
#' @export
sc_multiome_params <- function(study_md, gtf_path, fa_path) {
  if(!file.exists("ref/ref.gencode45.formatted.rds") && !exists("fa_genome")) { # nolint
    ## Format gene annotation .gtf
    ref1 <- gtf_path
    ref_gene1 <- rtracklayer::import(ref1)
    ref_gene1$gene_biotype <- ref_gene1$gene_type
    GenomeInfoDb::seqlevelsStyle(ref_gene1) <- "UCSC"
    rformat <- GenomeInfoDb::keepStandardChromosomes(
      ref_gene1,
      pruning.mode = "coarse"
    )
    saveRDS(rformat, "ref/ref.gencode45.formatted.rds")
    ## Import genome FASTA
    fa_genome <- Rsamtools::FaFile(fa_path)
  }
  if(file.exists("ref/ref.gencode45.formatted.rds") && !exists("fa_genome")) { # nolint
    print("Formatted gene annotation already exists; using formatted file...")
    ## Import formatted gene annotation .gtf
    rformat <- readRDS("ref/ref.gencode45.formatted.rds")
    ## Import genome FASTA
    fa_genome <- Rsamtools::FaFile(fa_path)
  }

  list_params <- data.frame(
    # Universal columns
    data.frame(
      # Sample Number
      Sample.No = seq(1:length( # nolint
        basename(
          list.files("data/")
        )
      )
      ),
      # Individual file names (uses CellRanger folder name by default)
      File.ID = basename(list.files("data/")),
      # Data file paths (Location of CellRanger files: 'Data/' by default)
      Path = paste("data/", list.files("data/"), sep = ""),
      # Path to count files
      Path.count = paste(
        "data/",
        list.files("data/"),
        "/filtered_feature_bc_matrix.h5",
        sep = ""
      ),
      # Path to raw count files (for ambient contamination removal)
      Path.raw = paste(
        "data/",
        list.files("data/"),
        "/raw_feature_bc_matrix.h5",
        sep = ""
      ),
      # Path to fragment files
      Path.frag = paste(
        "data/",
        list.files("data/"),
        "/atac_fragments.tsv.gz",
        sep = ""
      ),
      # Path to feature files
      Path.feat = paste(
        "data/",
        list.files("data/"),
        "/filtered_feature_bc_matrix/features.tsv.gz",
        sep = ""
      )
    ),

    # Dataset-specific columns
    study_md
  )
  return( # nolint
    list(
      "param" = list_params,
      "ref.gtf" = rformat,
      "ref.fa" = fa_genome
    )
  )
}
