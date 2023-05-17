#' @title readGseaFile
#' @name readGseaFile
#' @param filePath the path of the GSEA software enrichment outputs, defalut is NULL.
#'
#' @return a list contains meta(intergated enrichment results),
#' glist(the ordered gene lists), gset(all background enrichment terms).
#' @export

globalVariables(c("term","gene","filePath","coreEnrichment"))

# define function
readGseaFile <- function(filePath = NULL){
  # ==================================================================================
  # 1.prepare meta data
  # ==================================================================================
  enrich.report <- list.files(filePath,
                              # pattern = "^gsea_report.*.tsv",
                              pattern = "^gsea_report.*.(tsv|xls|xlsx)$",
                              full.names = TRUE)

  # check suffix
  suffix <- unlist(strsplit(enrich.report[1],split = "\\."))
  suffix <- suffix[length(suffix)]

  # read to table
  purrr::map_df(enrich.report,function(x){
    tmp <- suppressMessages(readr::read_tsv(x,show_col_types = FALSE))
    tmp <- tmp[,c(-3,-12)]

    # check data
    if(nrow(tmp) > 0){
      return(tmp)
    }else{
      return(NULL)
    }
  }) -> enrich.meta

  colnames(enrich.meta) <- c("ID","Description","setSize","enrichmentScore","NES",
                             "pvalue","p.adjust","qvalue","rank","leading_edge")

  # load ranked data
  # x = 1
  purrr::map_df(seq_along(1:nrow(enrich.meta)),function(x){
    # check whether file exist
    dir.file <- paste(filePath,enrich.meta$ID[x],".",suffix,sep = "")
    if(file.exists(dir.file)){
      tmp <- suppressMessages(readr::read_tsv(dir.file,show_col_types = FALSE)) %>%
        dplyr::mutate(Description = enrich.meta$ID[x],
                      id = enrich.meta$ID[x])

      # tmp <- tmp[,c(-3,-8)]
      tmp <- cbind(tmp[,c(1:2)],tmp[,c("RANK IN GENE LIST","RANK METRIC SCORE",
                                       "RUNNING ES","CORE ENRICHMENT","Description","id")])
    }
    # return(tmp)
  }) -> enrich.rank

  colnames(enrich.rank)[1:6] <- c("name","symbol","rank","metricScore","runningScore",
                                  "coreEnrichment")

  # add core enrichment to enrich.meta
  sid <- unique(enrich.rank$id)
  # x = 1
  purrr::map_df(1:length(sid),function(x){
    tmp.meta <- enrich.meta %>%
      dplyr::filter(ID == sid[x])
    tmp.core <- enrich.rank %>%
      dplyr::filter(id == sid[x] & coreEnrichment == "Yes") %>%
      dplyr::arrange(rank)

    # add core enrichment
    tmp.meta <- tmp.meta %>%
      dplyr::mutate(core_enrichment = paste(tmp.core$symbol,collapse = "/"))
    return(tmp.meta)
  }) -> enrich.meta.new

  # ==================================================================================
  # 2.prepare ranked data
  # ==================================================================================
  # preprare geneList
  enrich.gene <- list.files(filePath,
                            # pattern = "^ranked_gene.*.tsv",
                            pattern = "^ranked_gene.*.(tsv|xls|xlsx)$",
                            full.names = T)

  genelist <- utils::read.table(enrich.gene,sep = "\t",header = TRUE)
  glist <- as.numeric(genelist$SCORE)
  names(glist) <- genelist$NAME

  # prepare geneSets
  gset <- clusterProfiler::read.gmt(paste(filePath,"edb/gene_sets.gmt",sep = ""))

  # return
  res <- list(meta = enrich.meta.new,
              glist = glist,
              gset = gset)

  class(res) <- "GSEAfromSoft"
  return(res)
}
