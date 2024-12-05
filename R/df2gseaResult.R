globalVariables(c("SCORE"))

#' Create gseaResult Object from Enrichment Analysis Output Dataframe
#'
#' This function creates a gseaResult object from the results of enrichment analysis
#' using gene set enrichment analysis (GSEA) for a given gene list and Gene Ontology (GO) data.
#'
#' @param enrich.df A data frame containing the results of enrichment analysis.
#' @param geneList A decreasing sorted vector of gene list.
#' @param OrgDb The organism-specific annotation database.
#' @param keytype The type of gene identifier used in the analysis (default is "ENTREZID").
#' @param setType The type of GO ontology to use ("BP" for Biological Process, "CC" for Cellular Component,
#'                "MF" for Molecular Function, or "ALL" for all ontologies).
#' @param pvalueCutoff The p-value cutoff for significance (default is 0.05).
#' @param eps A small value to avoid division by zero (default is 1e-10).
#' @param pAdjustMethod The p-value adjustment method (default is "BH" for Benjamini-Hochberg).
#' @param exponent The exponent for weighting p-values (default is 1).
#' @param minGSSize The minimum gene set size (default is 10).
#' @param maxGSSize The maximum gene set size (default is 500).
#'
#' @import yulab.utils gson AnnotationDbi  GO.db
#'
#' @return A gseaResult object containing the results of GSEA.
#'
#' @export
dfGO2gseaResult <- function(enrich.df = NULL,
                            geneList = NULL,
                            OrgDb = NULL,
                            keytype = "ENTREZID",
                            setType = c("BP","CC","MF","ALL"),
                            pvalueCutoff = 0.05,
                            eps = 1e-10,
                            pAdjustMethod = "BH",
                            exponent = 1,
                            minGSSize = 10,
                            maxGSSize = 500){

  setType <- match.arg(setType,c("BP","CC","MF","ALL"))

  # check geneList
  if (!is.sorted(geneList))
    stop("geneList should be a decreasing sorted vector...")

  # get genesets data
  GO_DATA <- get_GO_data(OrgDb = OrgDb, ont = setType, keytype = keytype)
  geneSets <- getGeneSet(GO_DATA)

  if(keytype == 'SYMBOL'){
    readable <- TRUE
  }else{
    readable <- FALSE
  }

  organism <- get_organism(OrgDb)

  # create gseaResult class
  gseaResult <- methods::new("gseaResult",
                             result = enrich.df,
                             organism = organism,
                             setType = setType,
                             geneSets = geneSets,
                             geneList = geneList,
                             keytype = keytype,
                             permScores = matrix(),
                             params = list(pvalueCutoff = pvalueCutoff,
                                           eps = eps,
                                           pAdjustMethod = pAdjustMethod,
                                           exponent = exponent,
                                           minGSSize = minGSSize,
                                           maxGSSize = maxGSSize),
                             gene2Symbol = character(),
                             readable = readable,
                             termsim = matrix(),
                             method = character(),
                             dr = list())
}



#' Create gseaResult Object from KEGG Enrichment Analysis Output Dataframe
#'
#' This function creates a gseaResult object from the results of KEGG enrichment analysis
#' for a given gene list and organism.
#'
#' @param enrich.df A data frame containing the results of KEGG enrichment analysis.
#' @param geneList A decreasing sorted vector of gene identifiers.
#' @param organism The organism for KEGG enrichment analysis (default is "hsa" for Homo sapiens).
#' @param keytype The type of gene identifier used in the analysis (default is "kegg").
#' @param setType The type of enrichment analysis ("KEGG" for KEGG pathways).
#' @param use_internal_data Logical value indicating whether to use internal KEGG data (default is FALSE).
#' @param pvalueCutoff The p-value cutoff for significance (default is 0.05).
#' @param eps A small value to avoid division by zero (default is 1e-10).
#' @param pAdjustMethod The p-value adjustment method (default is "BH" for Benjamini-Hochberg).
#' @param exponent The exponent for weighting p-values (default is 1).
#' @param minGSSize The minimum gene set size (default is 10).
#' @param maxGSSize The maximum gene set size (default is 500).
#'
#' @return A gseaResult object containing the results of KEGG enrichment analysis.
#'
#' @export
dfKEGG2gseaResult <- function(enrich.df = NULL,
                              geneList = NULL,
                              organism = "hsa",
                              keytype = "kegg",
                              setType = "KEGG",
                              use_internal_data = FALSE,
                              pvalueCutoff = 0.05,
                              eps = 1e-10,
                              pAdjustMethod = "BH",
                              exponent = 1,
                              minGSSize = 10,
                              maxGSSize = 500){

  # check geneList
  if (!is.sorted(geneList))
    stop("geneList should be a decreasing sorted vector...")

  # get genesets data
  if (inherits(organism, "character")) {
    if (organism == "cpd") {
      organism = gson_cpd()
    }
  }

  if (inherits(organism, "character")) {
    species <- organismMapper(organism)
    if (use_internal_data) {
      KEGG_DATA <- get_data_from_KEGG_db(species)
    } else {
      KEGG_DATA <- prepare_KEGG(species, "KEGG", keytype)
    }
  } else if (inherits(organism, "GSON")) {
    KEGG_DATA <- organism
    species <- KEGG_DATA@species
    keyType <- KEGG_DATA@keytype
  } else {
    stop("organism should be a species name or a GSON object")
  }

  geneSets <- getGeneSet(KEGG_DATA)


  # create gseaResult class
  gseaResult <- methods::new("gseaResult",
                             result = enrich.df,
                             organism = organism,
                             setType = setType,
                             geneSets = geneSets,
                             geneList = geneList,
                             keytype = keytype,
                             permScores = matrix(),
                             params = list(pvalueCutoff = pvalueCutoff,
                                           eps = eps,
                                           pAdjustMethod = pAdjustMethod,
                                           exponent = exponent,
                                           minGSSize = minGSSize,
                                           maxGSSize = maxGSSize),
                             gene2Symbol = character(),
                             readable = FALSE,
                             termsim = matrix(),
                             method = character(),
                             dr = list())
}
