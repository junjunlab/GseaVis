#' @title gsInfo
#' @name gsInfo
#' @author Guangchuang Yu
#' @param object gseaResult object
#' @param geneSetID gene set ID
#'
#' @return data.frame
#' @export

# define function
gsInfo <- function(object, geneSetID) {
  gseaScores <- utils::getFromNamespace("gseaScores", "DOSE")

  geneList <- object@geneList

  if (is.numeric(geneSetID)) {
    geneSetID <- object@result[geneSetID, "ID"]
  }

  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList

  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
