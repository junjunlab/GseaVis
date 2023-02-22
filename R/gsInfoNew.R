#' @title gsInfoNew
#' @name gsInfoNew
#' @author Guangchuang Yu, modified by JunZhang
#' @param geneList geneList for GSEA software outputs which is saved in
#'  "gseaOutputs/ranked_gene_list_treat_versus_control*.tsv".
#' @param geneSetID gene set ID.
#' @param geneSet enrichment term sets for GSEA software outputs
#' which is saved in "*/edb/gene_sets".
#' @param exponent weight of each step, defalut is 1.
#' @return a data.frame
#' @export

# define function
gsInfoNew <- function(geneList = NULL,
                      geneSetID = NULL,
                      geneSet = NULL,
                      exponent = 1) {
  gseaScores <- utils::getFromNamespace("gseaScores", "DOSE")

  # geneList <- object@geneList
  #
  # if (is.numeric(geneSetID)) {
  #   geneSetID <- object@result[geneSetID, "ID"]
  # }

  # geneSet <- object@geneSets[[geneSetID]]
  geneSet.new <- geneSet %>%
    dplyr::filter(term == geneSetID) %>%
    dplyr::select(gene)
  geneSet.new <- geneSet.new$gene
  # names(geneSet.new) <- geneSetID
  # exponent <- object@params[["exponent"]]

  df <- gseaScores(geneList, geneSet.new, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList

  # df$Description <- object@result[geneSetID, "Description"]
  df$Description <- geneSetID
  return(df)
}
