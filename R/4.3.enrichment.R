#### KS-based enrichment analysis (Gene Ontology or custom gene list) ####

#' scRNA-Seq Enrichment Analysis
#'
#' Performs enrichment analysis based on the Kolmogorov-Smirnov Test for a chosen list of DGEA results.
#'
#' @param l_deg A dgea.output data frame created by sc.DGEA(). See sc.DGEA function documentation for details.
#' @param en_type Should enrichment statistics be calculated for GO terms or for a custom gene set list?
#'  Type "GO" for GO enrichment or "cstm" to use a custom gene set list. If using a custom list, 
#'  provide a text file containing gene names in column 1 and set IDs in column 2.
#' @param cstm_list Only used when conducting enrichment analysis based on a custom gene set list. 
#' Provide the path (relative to the working directory) to a custom gene set list for calculating enrichment statistics.
#' @param parl Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core_perc Percentage of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A list of Gene Ontology or custom gene set enrichment results per cell type for all comparisons
#' @examples
#'
#' # enr <- sc_enrichment(
#' # # Data frame containing DGEA results (created from sc.DGEA)
#' # l.deg = readRDS("analysis/object.dgea.result.rds")[[1]],
#' # # Enrichment type (Either "GO" or "cstm")
#' # en.type = "cstm",
#' # # If enrichment type is 'cstm', 
#' # # provide a custom gene set list (ignored when performing GO enrichment)
#' # cstm.list = setNames(
#' #   read.table(
#' #     "ref/gene.list.lipids.txt",header = T,sep = "\t"),
#' #     c("GENE","Set")
#' #     ),
#' # # run in parallel? (Set to FALSE if on Windows)
#' # parl = TRUE,
#' # # core percentage to use
#' # core.perc = 0.75
#' # )
#'
#' @export
sc_enrichment <- function(
  l_deg, en_type, cstm_list, parl, core_perc
) {
  # If no custom gene list is provided:
  if(en_type == "GO") {
    # parameters
    ## l.deg
    d <- l_deg

    # map gene names to Ensembl
    gene.db <- setNames(biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol"),
      mart = biomaRt::useMart(
        "ensembl",
        dataset = "hsapiens_gene_ensembl")
      ),c("ID","GENE"))
    d <- dplyr::left_join(
      d,
      gene.db,
      by = "GENE"
      )
    d <- d[!is.na(d[["ID"]]),]
    
    # groupwise filtering by CellType
    enrich.filt <- unique(
      d[["CellType"]]
      )
    enrich.list <- lapply(
      enrich.filt,
      function(x)
        dplyr::filter(
          d,
          CellType == x
          )
      )
    # Remove duplicate values
    enrich.list <- setNames(lapply(
      enrich.list,
      function (x) {
        x[!duplicated(x[["GENE"]]),]
        }
      ),c(unique(d[["CellType"]])))
    
    # GO enrichment function
    ## Definitions:
    ### df - data frame to use for enrichment
    ### id1 - column with list of gene names
    ### var1 - variables to use for topGO creation (ensembl ID, fold-change, and p-value)
    ### ont1 - GO ontology to use (BP - biological process, MF - molecular function, or CC - cellular component)
    ### db1 - database to use *In most cases this will be human "org.Hs.eg.db"
    ### map1 - mapping to be used for annotating genes (use "ensembl")
    ### desc1 - character string with description of study
    fun.GO.enrich <- function(df,
                              id1,
                              var1,
                              ont1,
                              db1,
                              map1,
                              desc1) {
      library(topGO)
      # use gene name as row name
      rownames(df) <- df[[id1]]
      # Select gene, log2FC, and p-value columns
      d.go1 <- setNames(df[,var1],c(id1,"FC.log2","p"))
      
      # Create the topGO object
      ## genelist
      genelist <- setNames(d.go1[["p"]],c(d.go1[[id1]]))
      ## topGO object
      d.go2 <- new(
        "topGOdata",
        ontology = ont1,
        allGenes = genelist,
        geneSelectionFun = function(x) x < 0.10,
        annot = annFUN.org,
        mapping = db1,
        ID = map1,
        nodeSize = 5
        )
      # Check sig
      go.sig <- sample(
        topGO::usedGO(d.go2),
        10
        )
      go.terms <- topGO::termStat(
        d.go2,go.sig
        )
      
      # perform enrichment statistics
      ## Update study description
      description(d.go2) <- paste(
        description(d.go2),
        desc1
        )
      ## KS test (Kolmogorov-Smirnov Test)
      KS.enrich <- topGO::runTest(
        d.go2,
        algorithm = "weight01",
        statistic = "KS"
        )
      
      KS.tab <- topGO::GenTable(
        d.go2,
        `P-Value` = KS.enrich,
        topNodes = length(KS.enrich@score),
        numChar = 120
        )
      
      KS.tab2 <- dplyr::select(
        KS.tab,
        GO.ID,
        Term,
        Annotated,
        Significant,
        `P-Value`
        )
      ## calculate FDR adjusted p-values
      KS.tab.out <- KS.tab2
      ## keep pathways with raw p < 0.05 and at least 3 significant differentially expressed genes
      KS.tab.out <- dplyr::filter(
        KS.tab.out,
        `P-Value` < 0.05 &
          Significant >= 3
        )
      KS.tab.out[["Description"]] <- desc1
      ## Unlist topGO object to obtain genes contained within each pathway
      genes.in.term <- genesInTerm(d.go2)
      gene.terms <- KS.tab.out[["GO.ID"]]
      genes.in.term2 <- genes.in.term[c(gene.terms)]
      genes.in.term3 <- reshape2::melt(genes.in.term2)
      names(genes.in.term3) <- c(id1, "GO.ID")
      gene.sig <- sigGenes(d.go2)
      ## return list of GO terms for each gene
      genes.output <- merge(
        genes.in.term3,
        d.go1,
        by = id1
        )
      genes.output.sig <- dplyr::left_join(
        df,
        dplyr::filter(
          genes.output,
          p < 0.05
        ),
        by = id1
        )
      genes.output.sig[["Description"]] <- desc1
      genes.output.sig <- dplyr::left_join(
        genes.output.sig,
        KS.tab.out[
          ,
          c("GO.ID","Term")
        ],
        by = "GO.ID"
        )
      return(
        list(
          "GO.output" = KS.tab.out,
          "GO.input" = genes.output.sig
          )
        )
      }
    
    # Run for all celltypes and ontologies
    if(parl == TRUE && Sys.info()[["sysname"]] != "Windows"){
      # Groupwise GO enrichment by Cell Type function
      fun.GO.run.ct <- function(
        l1,
        z
        ) {
        library(org.Hs.eg.db)
        ggo <- parallel::mclapply(
          mc.cores = ceiling(
            parallel::detectCores()*core_perc
            ),
          l1,
          function (y)
            fun.GO.enrich(
              y,
              "ID",
              c("ID","log2FC","H.qval"),
              z,
              "org.Hs.eg.db",
              "ensembl",
              paste("GO",
                    z,
                    unique(y[["CellType"]]),
                    sep = "_")
              )
          )
          return(ggo)
          }
    }
    
    if(parl == FALSE & Sys.info()[["sysname"]] != "Windows"){
      # Groupwise GO enrichment by Cell Type function
      fun.GO.run.ct <- function(
    l1,
    z
      ) {
        library(org.Hs.eg.db)
        ggo <- lapply(
          l1,
          function (y)
            fun.GO.enrich(
              y,
              "ID",
              c("ID","log2FC","H.qval"),
              z,
              "org.Hs.eg.db",
              "ensembl",
              paste("GO",
                    z,
                    unique(y[["CellType"]]),
                    sep = "_")
            )
        )
        return(ggo)
      }
    }
    
    if(parl == TRUE & Sys.info()[["sysname"]] == "Windows"){
      print("Windows OS detected; defaulting to sequential processing...")
      # Groupwise GO enrichment by Cell Type function
      fun.GO.run.ct <- function(
    l1,
    z
      ) {
        library(org.Hs.eg.db)
        ggo <- lapply(
          l1,
          function (y)
            fun.GO.enrich(
              y,
              "ID",
              c("ID","log2FC","H.qval"),
              z,
              "org.Hs.eg.db",
              "ensembl",
              paste("GO",
                    z,
                    unique(y[["CellType"]]),
                    sep = "_")
            )
        )
        return(ggo)
      }
    }
    
    if(parl == FALSE & Sys.info()[["sysname"]] != "Windows"){
      # Groupwise GO enrichment by Cell Type function
      fun.GO.run.ct <- function(
    l1,
    z
      ) {
        library(org.Hs.eg.db)
        ggo <- lapply(
          l1,
          function (y)
            fun.GO.enrich(
              y,
              "ID",
              c("ID","log2FC","H.qval"),
              z,
              "org.Hs.eg.db",
              "ensembl",
              paste("GO",
                    z,
                    unique(y[["CellType"]]),
                    sep = "_")
            )
        )
        return(ggo)
      }
    }
      # Run with all ontologies
      GO.output <- list(
        ## Biological Process
        "Biological Process" = fun.GO.run.ct(
          enrich.list,
          "BP"
          ),
        ## Molecular Function
        "Molecular Function" = fun.GO.run.ct(
          enrich.list,
          "MF"
          ),
        ## Cell Component
        "Cell Component" = fun.GO.run.ct(
          enrich.list,
          "CC"
          )
        )
      
      # Format results
      go1 <- list(
        "Results" = dplyr::bind_rows(
          lapply(
            seq.int(1,length(GO.output),1),
            function(x) dplyr::bind_rows(
              lapply(
                seq.int(1,length(GO.output[[x]]),1),
                function(y) dplyr::bind_rows(
                  GO.output[[x]][[y]][["GO.output"]])
              )
            )
          )
        ),
        "Input" = dplyr::bind_rows(
          lapply(
            seq.int(1,length(GO.output),1),
            function(x) dplyr::bind_rows(
              setNames(lapply(
                seq.int(1,length(GO.output[[x]]),1),
                function(y) dplyr::bind_rows(
                  GO.output[[x]][[y]][["GO.input"]])
              ),c(names(GO.output[[x]]))),
              .id = "CellType"
            )
          )
        )
      )
      
      go1.dir <- dplyr::count(
        dplyr::mutate(
          go1[["Input"]],
          "inc" = go1[["Input"]][["log2FC"]] > 0),
        .data[["Description"]],
        .data[["Term"]],
        .data[["inc"]])
      go1.dir <- dplyr::full_join(
        go1.dir[go1.dir[["inc"]] == TRUE,][-3],
        go1.dir[go1.dir[["inc"]] == FALSE,][-3],
        by = c("Description","Term")
        )
      go1.dir[is.na(go1.dir[["n.x"]]),"n.x"] <- 0
      go1.dir[is.na(go1.dir[["n.y"]]),"n.y"] <- 0
      go1.dir <- go1.dir[!duplicated(go1.dir[,c("Description","Term")]),]
      go1.dir[["ratio"]] <- ifelse(
          go1.dir[["n.y"]] == 0,
          1,
          (go1.dir[["n.x"]])/(go1.dir[["n.x"]]+go1.dir[["n.y"]])
          )
      
      go1.comb <- data.frame(
        "Description" = unique(go1[["Input"]][,c("Description","Term")])[["Description"]],
        "Term" = unique(go1[["Input"]][,c("Description","Term")])[["Term"]]
        )
      go1.comb <- dplyr::left_join(
        go1.comb,
        go1.dir[,c("Description","Term","ratio")],
        by = c("Description","Term")
        )
      go1[["Results"]] <- dplyr::left_join(
        go1[["Results"]],
        go1.comb,
        by = c("Description","Term")
        )
      names(go1[["Results"]]) <- c("GO.ID","Term","cl.size","sig","p","Description","ratio")
      go1[["Results"]][["p"]] <- as.numeric(
        ifelse(
          go1[["Results"]][["p"]] == "< 1e-30",
          1e-30,
          go1[["Results"]][["p"]])
        )
      
      write.table(
        go1[["Results"]],
        "analysis/table.enrichment.GO.results.txt",
        sep = "\t",
        col.names = T,
        row.names = F
        )
      
      write.table(
        go1[["Input"]],
        "analysis/table.enrichment.GO.input.txt",
        sep = "\t",
        col.names = T,
        row.names = F
        )
      
      return(go1[["Results"]])
      
    }
  
  ## If a custom list is provided:
  if(en_type == "cstm" & missing(cstm_list)){print("Error: must provide a custom gene set list!")}
  if(en_type == "cstm" & !missing(cstm_list)) {
    # params
    ## l.deg
    l1 <- l_deg
    ## cstm.list
    l2 <- cstm_list
    
    # format df
    cr <- dplyr::bind_rows(
      lapply(
        levels(l1[["CellType"]]),
        function(x){
          d1 <- dplyr::left_join(
            l2,
            l1[l1[["CellType"]] == x,],
            by = c("GENE")
            )
          }
        )
      )
    cr2 <- list(
      "Input" = cr[!is.na(cr[["CellType"]]),
                   c("GENE","Set","CellType","Comparison","log2FC","H.qval")],
      "Missing" = cr[is.na(cr[["CellType"]]),c("GENE","Set","CellType","Comparison")])
    
    # Conduct KS test for each comparison and cell type
    if(parl == TRUE & Sys.info()[["sysname"]] != "Windows"){
      ks.vars <- list(
        "CellType" = tryCatch(
          {levels(cr2[["Input"]][["CellType"]])},
          error = function(e){print("Error: CellType variable must be a factor!")}),
        "Comparison" = unique(cr2[[1]][["Comparison"]]),
        "Set" = unique(cr2[[1]][["Set"]])
        )
      ks1 <- setNames(
        lapply(
          ks.vars[["Comparison"]],
          function(x){
            d1 <- setNames(
              parallel::mclapply(
                mc.cores = ceiling(
                  parallel::detectCores()*core_perc
                  ),
                ks.vars[["CellType"]],
                function(y){
                  d1 <- setNames(
                    lapply(
                      ks.vars[["Set"]],
                      function(z){
                        d <- tryCatch(
                          {cr2[["Input"]][
                            cr2[["Input"]][["CellType"]] == y &
                              cr2[["Input"]][["Comparison"]] == x &
                              cr2[["Input"]][["Set"]] == z,
                            c("GENE","Set","Comparison","CellType","H.qval","log2FC")]},
                          error = function(e){print("Warning: no genes have calculated fold changes for the selected set; skipping to next gene set...")})
                        d1 <- tryCatch(
                          {
                            d1 <- ks.test(
                              d[,5:6],
                              "punif",
                              alternative = "greater")
                            d1 <- data.frame(
                              "Comparison" = x,
                              "CellType" = y,
                              "Set" = z,
                              "genes" = paste(d[["GENE"]],collapse = ","),
                              "cl.size" = length(d[["GENE"]]),
                              "inc" = length(d[d[["log2FC"]] > 0,"GENE"]),
                              "dec" = length(d[d[["log2FC"]] < 0,"GENE"]),
                              "p.value" = d1[["p.value"]]
                              )
                            d1 <- cbind(
                              d1,
                              data.frame(
                                "ratio" = ifelse(
                                  d1[["dec"]] == 0,
                                  1,
                                  (d1[["inc"]])/(d1[["inc"]]+d1[["dec"]])
                                  )
                                )
                              )
                          },
                          error = function(e){print("Warning: KS test unsuccessful for selected set (likely due to insufficient genes present in a set); skipping to next gene set...")})
                        return(d1)
                        }
                      ),
                    c(ks.vars[["Set"]])
                    )
                  }),
              c(ks.vars[["CellType"]])
              )
          }),
        c(ks.vars[["Comparison"]])
        )
      }
    
    if(parl == TRUE & Sys.info()[["sysname"]] == "Windows"){
      print("Windows OS detected; defaulting to sequential processing...")
      ks.vars <- list(
        "CellType" = tryCatch(
          {levels(cr2[["Input"]][["CellType"]])},
          error = function(e){print("Error: CellType variable must be a factor!")}),
        "Comparison" = unique(cr2[[1]][["Comparison"]]),
        "Set" = unique(cr2[[1]][["Set"]])
      )
      ks1 <- setNames(
        lapply(
          ks.vars[["Comparison"]],
          function(x){
            d1 <- setNames(
              lapply(
                ks.vars[["CellType"]],
                function(y){
                  d1 <- setNames(
                    lapply(
                      ks.vars[["Set"]],
                      function(z){
                        d <- tryCatch(
                          {cr2[["Input"]][
                            cr2[["Input"]][["CellType"]] == y &
                              cr2[["Input"]][["Comparison"]] == x &
                              cr2[["Input"]][["Set"]] == z,
                            c("GENE","Set","Comparison","CellType","H.qval","log2FC")]},
                          error = function(e){print("Warning: no genes have calculated fold changes for the selected set; skipping to next gene set...")})
                        d1 <- tryCatch(
                          {
                            d1 <- ks.test(
                              d[,5:6],
                              "punif",
                              alternative = "greater")
                            d1 <- data.frame(
                              "Comparison" = x,
                              "CellType" = y,
                              "Set" = z,
                              "genes" = paste(d[["GENE"]],collapse = ","),
                              "cl.size" = length(d[["GENE"]]),
                              "inc" = length(d[d[["log2FC"]] > 0,"GENE"]),
                              "dec" = length(d[d[["log2FC"]] < 0,"GENE"]),
                              "p.value" = d1[["p.value"]]
                            )
                            d1 <- cbind(
                              d1,
                              data.frame(
                                "ratio" = ifelse(
                                  d1[["dec"]] == 0,
                                  1,
                                  (d1[["inc"]])/(d1[["inc"]]+d1[["dec"]])
                                )
                              )
                            )
                          },
                          error = function(e){print("Warning: KS test unsuccessful for selected set (likely due to insufficient genes present in a set); skipping to next gene set...")})
                        return(d1)
                      }
                    ),
                    c(ks.vars[["Set"]])
                  )
                }),
              c(ks.vars[["CellType"]])
            )
          }),
        c(ks.vars[["Comparison"]])
      )
    }
    
    if(parl == FALSE & Sys.info()[["sysname"]] != "Windows"){
      ks.vars <- list(
        "CellType" = tryCatch(
          {levels(cr2[["Input"]][["CellType"]])},
          error = function(e){print("Error: CellType variable must be a factor!")}),
        "Comparison" = unique(cr2[[1]][["Comparison"]]),
        "Set" = unique(cr2[[1]][["Set"]])
      )
      ks1 <- setNames(
        lapply(
          ks.vars[["Comparison"]],
          function(x){
            d1 <- setNames(
              lapply(
                ks.vars[["CellType"]],
                function(y){
                  d1 <- setNames(
                    lapply(
                      ks.vars[["Set"]],
                      function(z){
                        d <- tryCatch(
                          {cr2[["Input"]][
                            cr2[["Input"]][["CellType"]] == y &
                              cr2[["Input"]][["Comparison"]] == x &
                              cr2[["Input"]][["Set"]] == z,
                            c("GENE","Set","Comparison","CellType","H.qval","log2FC")]},
                          error = function(e){print("Warning: no genes have calculated fold changes for the selected set; skipping to next gene set...")})
                        d1 <- tryCatch(
                          {
                            d1 <- ks.test(
                              d[,5:6],
                              "punif",
                              alternative = "greater")
                            d1 <- data.frame(
                              "Comparison" = x,
                              "CellType" = y,
                              "Set" = z,
                              "genes" = paste(d[["GENE"]],collapse = ","),
                              "cl.size" = length(d[["GENE"]]),
                              "inc" = length(d[d[["log2FC"]] > 0,"GENE"]),
                              "dec" = length(d[d[["log2FC"]] < 0,"GENE"]),
                              "p.value" = d1[["p.value"]]
                            )
                            d1 <- cbind(
                              d1,
                              data.frame(
                                "ratio" = ifelse(
                                  d1[["dec"]] == 0,
                                  1,
                                  (d1[["inc"]])/(d1[["inc"]]+d1[["dec"]])
                                )
                              )
                            )
                          },
                          error = function(e){print("Warning: KS test unsuccessful for selected set (likely due to insufficient genes present in a set); skipping to next gene set...")})
                        return(d1)
                      }
                    ),
                    c(ks.vars[["Set"]])
                  )
                }),
              c(ks.vars[["CellType"]])
            )
          }),
        c(ks.vars[["Comparison"]])
      )
    }
    
    if(parl == FALSE & Sys.info()[["sysname"]] == "Windows"){
      ks.vars <- list(
        "CellType" = tryCatch(
          {levels(cr2[["Input"]][["CellType"]])},
          error = function(e){print("Error: CellType variable must be a factor!")}),
        "Comparison" = unique(cr2[[1]][["Comparison"]]),
        "Set" = unique(cr2[[1]][["Set"]])
      )
      ks1 <- setNames(
        lapply(
          ks.vars[["Comparison"]],
          function(x){
            d1 <- setNames(
              lapply(
                ks.vars[["CellType"]],
                function(y){
                  d1 <- setNames(
                    lapply(
                      ks.vars[["Set"]],
                      function(z){
                        d <- tryCatch(
                          {cr2[["Input"]][
                            cr2[["Input"]][["CellType"]] == y &
                              cr2[["Input"]][["Comparison"]] == x &
                              cr2[["Input"]][["Set"]] == z,
                            c("GENE","Set","Comparison","CellType","H.qval","log2FC")]},
                          error = function(e){print("Warning: no genes have calculated fold changes for the selected set; skipping to next gene set...")})
                        d1 <- tryCatch(
                          {
                            d1 <- ks.test(
                              d[,5:6],
                              "punif",
                              alternative = "greater")
                            d1 <- data.frame(
                              "Comparison" = x,
                              "CellType" = y,
                              "Set" = z,
                              "genes" = paste(d[["GENE"]],collapse = ","),
                              "cl.size" = length(d[["GENE"]]),
                              "inc" = length(d[d[["log2FC"]] > 0,"GENE"]),
                              "dec" = length(d[d[["log2FC"]] < 0,"GENE"]),
                              "p.value" = d1[["p.value"]]
                            )
                            d1 <- cbind(
                              d1,
                              data.frame(
                                "ratio" = ifelse(
                                  d1[["dec"]] == 0,
                                  1,
                                  (d1[["inc"]])/(d1[["inc"]]+d1[["dec"]])
                                )
                              )
                            )
                          },
                          error = function(e){print("Warning: KS test unsuccessful for selected set (likely due to insufficient genes present in a set); skipping to next gene set...")})
                        return(d1)
                      }
                    ),
                    c(ks.vars[["Set"]])
                  )
                }),
              c(ks.vars[["CellType"]])
            )
          }),
        c(ks.vars[["Comparison"]])
      )
    }
    
    # Format results
    ks2 <- list(
      "Input" = cr2[["Input"]],
      "Results" = dplyr::bind_rows(
        lapply(
          seq.int(1,length(ks1),1),
          function(x) dplyr::bind_rows(
            lapply(
              seq.int(1,length(ks1[[x]]),1),
              function(y) dplyr::bind_rows(
                ks1[[x]][[y]][lengths(ks1[[x]][[y]]) > 1])
              )
            )
          )
        ),
      "Missing" = dplyr::bind_rows(
        lapply(
          seq.int(1,length(ks1),1),
          function(x) dplyr::bind_rows(
            setNames(lapply(
              seq.int(1,length(ks1[[x]]),1),
              function(y) dplyr::bind_rows(
                ks1[[x]][[y]][lengths(ks1[[x]][[y]]) == 1])
              ),c(names(ks1[[x]]))),
            .id = "CellType"
            )
          )
        )
      )
    
    # Save tables
    write.table(
      ks2[["Results"]],
      "analysis/table.enrichment.custom.list.txt",
      col.names = T,
      row.names = F,
      sep = "\t"
      )
    
    return(ks2)
  }
  
}





#' Enrichment Plots
#'
#' Plots the top 10 terms from a provided set of enrichment results.
#'
#' @param l_enr A data frame containing enrichment results generated by sc.enrichment().
#' @param en_type Are the provided enrichment results based on the Gene Ontology database or a custom gene set? Type "GO"
#' for Gene Ontology Results and "cstm" for custom gene set enrichment results.


#' @return An enrichment plot depicting the most significant gene sets for a given comparison and cell type.
#' @examples
#'
#' # # Custom gene set example
#' # sc_plot_enrichment(
#' #   enr[["Results"]][enr[["Results"]][["CellType"]] == "2.Secretory",],
#' #   "cstm"
#' # )
#' # # Gene ontology example
#' # sc_plot_enrichment(
#' #   enr[enr[["Description"]] == "GO_BP_2.Secretory",],
#' #   "GO"
#' # )
#' # 
#'
#' @export
sc_plot_enrichment <- function(
  l_enr,
  en_type
  ) {
  
  if(en_type == "cstm"){
    l1 <- l_enr
    if(nrow(l1) > 10){
      l1 <- dplyr::slice_max(
        l1,
        order_by = -.data[["p.value"]],
        n = 10
        )
      }
    if(nrow(l1) < 10){
      l1 <- dplyr::slice_max(
        l1,
        order_by = -.data[["p.value"]],
        n = nrow(l1)
        )
      }
  # Main Plot
  cr.plot <- ggplot2::ggplot(
    l1,
    ggplot2::aes(
      x = Set,
      y = -log10(p.value),
      fill = ratio
    )
  ) +
    ## Connect lines to significantly altered clusters
    ggplot2::geom_linerange(
      data = l1[l1[["p.value"]] < 0.05,],
      ggplot2::aes(
        x = l1[l1[["p.value"]] < 0.05,"Set"],
        ymax = -log10(l1[l1[["p.value"]] < 0.05,"p.value"]),
        ymin = 0
      ),
      color = "black"
    ) +
    ## Plot points
    ggplot2::geom_point(
      shape = 21,
      col = "black",
      ggplot2::aes(
        size = .data[["cl.size"]]
      ),
      alpha = ifelse(
        l1[["p.value"]] < 0.05,
        1,
        0.5
      )
    ) +
    ## Significance line
    ggplot2::geom_hline(
      yintercept = 1.3,
      linetype = "dashed",
      color = "firebrick4"
    ) +
    ggplot2::geom_hline(
      yintercept = 3.3,
      linetype = "dashed",
      color = "firebrick1"
    ) +
    shadowtext::geom_shadowtext(
      ggplot2::aes(1,
      1.3),
      label = "p < 0.05",
      vjust = -0.5,
      bg.color = "firebrick4",
      color = "white",
      size = 3
      ) +
    shadowtext::geom_shadowtext(
      ggplot2::aes(1,
                   3.3),
      label = "p < 0.0005",
      vjust = -0.5,
      bg.color = "firebrick1",
      color = "white",
      size = 3
    ) +
    ## Theme
    sc_theme1() +
    ggplot2::theme(plot.margin = grid::unit(
      c(1,4,1,4),
      "cm"
      ),
      axis.text.x = ggplot2::element_text(
        size = 10)) +
    ## Axis labels
    ggplot2::xlab("Gene Set") +
    ggplot2::ylab("-log10(adjusted p-value)") +
    ## Gradient and size scaling
    ggplot2::scale_fill_gradientn(
      name = "Increased Ratio",
      colors = c(
        'midnightblue',
        'royalblue2',
        "white",
        'goldenrod1',
        "firebrick3"
      ),
      breaks = c(
        0,.25,.5,.75,1
      ),
      labels = c(
        "100% Decreased",
        "25% Increased:75% Decreased",
        "50% Increased:50% Decreased",
        "75% Increased:25% Decreased",
        "100% Increased"
      ),
      limits = c(0,1)
    ) +
    ggplot2::scale_size_area(
      name = "Cluster Size",
      max_size = 16
    )
  
  # Save
  ggplot2::ggsave(
    paste(
      "analysis/",
      "plot.enrichment.custom.",
      unique(l1[["Comparison"]]),".",
      unique(l1[["CellType"]]),
      ".png",
      sep = ""
    ),
    width = 16,
    height = 10,
    dpi = 1000
  )}
  
  if(en_type == "GO"){
    l1 <- l_enr
    l1 <- dplyr::slice_max(
      l1,
      order_by = -.data[["p"]],
      n = 10
      )
    # Main Plot
    cr.plot <- ggplot2::ggplot(
      l1,
      ggplot2::aes(
        x = Term,
        y = -log10(p),
        fill = ratio
      )
    ) +
      ## Connect lines to significantly altered clusters
      ggplot2::geom_linerange(
        data = l1[l1[["p"]] < 0.05,],
        ggplot2::aes(
          x = l1[l1[["p"]] < 0.05,"Term"],
          ymax = -log10(l1[l1[["p"]] < 0.05,"p"]),
          ymin = 0
        ),
        color = "black"
      ) +
      ## Plot points
      ggplot2::geom_point(
        shape = 21,
        col = "black",
        ggplot2::aes(
          size = .data[["cl.size"]]
        ),
        alpha = ifelse(
          l1[["p"]] < 0.05,
          1,
          0.5
        )
      ) +
      ## Significance line
      ggplot2::geom_hline(
        yintercept = 1.3,
        linetype = "dashed",
        color = "firebrick4"
      ) +
      ggplot2::geom_hline(
        yintercept = 3.3,
        linetype = "dashed",
        color = "firebrick1"
      ) +
      shadowtext::geom_shadowtext(
        ggplot2::aes(1,
                     1.3),
        label = "p < 0.05",
        vjust = -0.5,
        bg.color = "firebrick4",
        color = "white",
        size = 3
      ) +
      shadowtext::geom_shadowtext(
        ggplot2::aes(1,
                     3.3),
        label = "p < 0.0005",
        vjust = -0.5,
        bg.color = "firebrick1",
        color = "white",
        size = 3
      ) +
      ## Theme
      sc_theme1() +
      ggplot2::theme(plot.margin = grid::unit(
        c(1,4,1,4),
        "cm"
      ),
      axis.text.x = ggplot2::element_text(
        size = 10)) +
      ## Axis labels
      ggplot2::xlab("Gene Set (GO)") +
      ggplot2::ylab("-log10(adjusted p-value)") +
      ## Gradient and size scaling
      ggplot2::scale_fill_gradientn(
        name = "Increased Ratio",
        colors = c(
          'midnightblue',
          'royalblue2',
          "white",
          'goldenrod1',
          "firebrick3"
        ),
        breaks = c(
          0,.25,.5,.75,1
        ),
        labels = c(
          "100% Decreased",
          "25% Increased:75% Decreased",
          "50% Increased:50% Decreased",
          "75% Increased:25% Decreased",
          "100% Increased"
        ),
        limits = c(0,1)
      ) +
      ggplot2::scale_size_area(
        name = "Cluster Size",
        max_size = 16
      )
    
    # Save
    ggplot2::ggsave(
      paste(
        "analysis/",
        "plot.enrichment.GO.",
        unique(l1[["Description"]]),
        ".png",
        sep = ""
      ),
      width = 16,
      height = 10,
      dpi = 1000
    )
    
  }
}
