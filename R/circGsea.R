globalVariables(c("symbol"))

#' Circular Gene Set Enrichment Analysis Plot
#'
#' This function creates a circular plot for gene set enrichment analysis.
#'
#' @param object An object containing enrichment results from clusterProfiler.
#' @param geneSetID A vector of gene set IDs to be plotted.
#' @param htCol A vector of two colors for the heatmap gradient. Default is c("#08519C", "#A50F15").
#' @param segmentsCol A vector of two colors for the segment gradient. Default is c("#0099CC", "#FF0033").
#' @param bgCol Background color for the sectors. Default is "grey80".
#' @param curveCol Color for the runningScore curve. Default is "#CC3333".
#' @param pointSize Size of the points on the runningScore curve. Default is 0.5.
#' @param addPoints Logical indicating whether to add points to the runningScore curve. Default is FALSE.
#' @param type Type of plot to display: 'c' for heatmap only, 'h' for segments and heatmap, 'm' for segments only. Default is "c".
#' @param addGeneRank Logical indicating whether to add a track for gene ranks. Default is FALSE.
#' @param markGene A vector of gene symbols to be marked on the plot. Default is NULL.
#' @param markGeneSize Size of the marked gene symbols. Default is 0.5.
#' @param GoIDfacing Facing direction of the gene set ID labels. Default is "inside".
#' Additional choice are "inside", "outside", "reverse.clockwise", "clockwise","downward", "bending",
#' "bending.inside", "bending.outside".
#' @param GoIDsize Size of the gene set ID labels. Default is 0.75.
#' @param addDescription Logical indicating whether to add term descriptions. Default is FALSE.
#' @param descriptionFacing Facing direction of the term description text. Default is "bending.outside".
#' Additional choice are "inside", "outside", "reverse.clockwise", "clockwise","downward", "bending",
#' "bending.inside", "bending.outside".
#' @param descripShift Shift of the term description text. Default is 1.
#' @param descripGap Gap between lines of the term description text. Default is 0.3.
#' @param descripLength Maximum line length for the term description text. Default is 40.
#' @param descripSize Size of the term description text. Default is 0.75.
#'
#' @return A circular gene set enrichment analysis plot.
#'
#' @import circlize
#'
#' @export
circGsea <- function(object = NULL,
                     geneSetID = NULL,
                     htCol = c("#08519C", "#A50F15"),
                     segmentsCol = c("#0099CC", "#FF0033"),
                     bgCol = "grey80",
                     curveCol = "#CC3333",
                     pointSize = 0.5,
                     addPoints = FALSE,
                     type = c("c","h","m"),
                     addGeneRank = FALSE,
                     markGene = NULL,
                     markGeneSize = 0.5,
                     GoIDfacing = "inside",
                     GoIDsize = 0.75,
                     addDescription = FALSE,
                     descriptionFacing = "bending.outside",
                     descripShift = 1,
                     descripGap = 0.3,
                     descripLength = 40,
                     descripSize = 0.75){
  type <- match.arg(type,c("c","h","m"))
  # ============================================================================
  # 1_extract data
  # ============================================================================

  # all rank
  # setid <- "GO:0002705"
  gsdata <- purrr::map_df(geneSetID,function(setid){
    tmp <- gsInfo(object,geneSetID = setid) %>%
      dplyr::mutate("gene_name" = names(object@geneList)) %>%
      dplyr::mutate(id = setid)

    # chenge id name
    if(length(object@gene2Symbol) != 0){
      tmp$symbol <- object@gene2Symbol
    }else{
      tmp$symbol <- ""
    }

    return(tmp)
  })


  # target gene rank
  gsdata1 <- purrr::map_df(unique(gsdata$Description),function(setid){
    tmp <- gsdata %>%
      dplyr::filter(Description == setid) %>%
      # dplyr::mutate("gene_name" = names(object@geneList)) %>%
      dplyr::filter(position == 1)
  })

  # heatmap data
  ht <- purrr::map_df(unique(gsdata$id),function(setid){
    tmp <- gsdata %>%
      dplyr::filter(id == setid)

    v <- seq(1, sum(tmp$position), length.out = 9)
    inv <- findInterval(rev(cumsum(tmp$position)), v)
    if (min(inv) == 0) {
      inv <- inv + 1
    }

    # new color
    color <- grDevices::colorRampPalette(c(htCol[1], "white", htCol[2]))(10)

    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(xmin = xmin,
                    xmax = xmax,
                    col = color[unique(inv)],
                    id = setid)
  })

  # enrichment dataframe
  df_lb <- data.frame(object) %>%
    dplyr::filter(ID %in% geneSetID)

  # filter target mark genes
  if(!is.null(markGene)){
    mgene_df <- gsdata %>%
      dplyr::filter(symbol %in% markGene)
  }

  # ============================================================================
  # 2_plot
  # ============================================================================

  circos.clear()
  # circos.initialize(sectors = gsdata1$id,x = gsdata1$x)
  circos.par(circle.margin = rep(0.1,4))
  circos.initialize(sectors = gsdata1$id,xlim = c(1,max(gsdata$x)))

  # first GO id
  circos.track(sectors = geneSetID,ylim = c(0,1),
               bg.col = bgCol,
               panel.fun = function(x, y){
                 tmp <- subset(df_lb,ID == CELL_META$sector.index)

                 # assign pvalue
                 pval <- tmp$pvalue
                 if(pval < 0.05 & pval >= 0.01){
                   plabel <- "< 0.05"
                 }else if(pval < 0.01 & pval >= 0.001){
                   plabel <- "< 0.01"
                 }else if(pval < 0.001){
                   plabel <- "< 0.001"
                 }

                 goName <- paste(CELL_META$sector.index,"\n",
                                 "(",
                                 "NES:",round(tmp$NES,digits = 1)," | ",
                                 "Pvalue:",plabel,
                                 ")",sep = "")

                 circos.text(x = CELL_META$xcenter,y = CELL_META$ycenter,
                             labels = goName,
                             font = 2,
                             niceFacing = T,
                             facing = GoIDfacing,
                             cex = GoIDsize)

                 # add term Description
                 if(addDescription == TRUE){
                   str <- stringr::str_wrap(tmp$Description,width = descripLength)
                   str <- unlist(strsplit(str,split = "\n"))
                   gap <- rev(seq(1 - length(str)*descripGap,descripShift,descripGap))

                   # loop add description
                   lapply(seq_along(str),function(x){
                     circos.text(x = CELL_META$xcenter,y = CELL_META$ycenter + gap[x],
                                 labels = str[x],
                                 font = 2,
                                 niceFacing = T,
                                 facing = descriptionFacing,
                                 cex = descripSize)
                   })
                 }
               })

  # add lines and segments track
  circos.track(sectors = geneSetID,ylim = c(-1,1),
               panel.fun = function(x, y){
                 tmp <- subset(gsdata1,id == CELL_META$sector.index)

                 circos.segments(x0 = 0,x1 = max(gsdata1$x),y0 = 0,y1 = 0,
                                 lty = "dashed",col = "black")

                 # gene rank segment layer
                 if(type == "h"){
                   # assign colors for segments
                   colors <- grDevices::colorRampPalette(segmentsCol)(nrow(tmp))

                   # order colors according to value
                   sorted_colors <- colors[match(tmp$runningScore,sort(tmp$runningScore))]

                   circos.lines(x = c(1,tmp$x,max(gsdata$x)),y = c(0,tmp$runningScore,0),
                                type = "h",baseline = 0,col = sorted_colors)


                 }else if(type == "m"){
                   circos.segments(x0 = c(1,tmp$x,max(gsdata$x)),x1 = c(1,tmp$x,max(gsdata$x)),
                                   y0 = -0.25,y1 = 0.25,
                                   lty = "solid",col = "grey30")
                 }else{

                 }

                 # runningScore line plot
                 if(type == "h" & addPoints == TRUE){
                   circos.points(x = c(1,tmp$x,max(gsdata$x)),y = c(0,tmp$runningScore,0),
                                 col = sorted_colors,pch = 19,cex = pointSize)
                 }else{
                   circos.lines(x = c(1,tmp$x,max(gsdata$x)),y = c(0,tmp$runningScore,0),
                                lwd = 2,col = curveCol)
                 }


               })

  # add heatmap and segments track
  circos.track(sectors = geneSetID,ylim = c(0,1),
               panel.fun = function(x, y){

                 if(type == "c"){
                   # gene rank segment layer
                   tmp <- subset(gsdata1,id == CELL_META$sector.index)

                   circos.segments(x0 = tmp$x,x1 = tmp$x,
                                   y0 = 0,y1 = 1)

                   ytop <- 0.5
                 }else{
                   ytop <- 1
                 }

                 # gene rank heatmap layer
                 tmpht <- subset(ht,id == CELL_META$sector.index)

                 circos.rect(xleft = tmpht$xmin,xright = tmpht$xmax,
                             col = ggplot2::alpha(tmpht$col,alpha = 0.75),
                             border = ggplot2::alpha(tmpht$col,alpha = 0.75),
                             ybottom = rep(0,nrow(tmpht)),ytop = rep(ytop,nrow(tmpht)))

                 circos.rect(xleft = 0,xright = max(tmpht$xmax),
                             col = NA,border = "black",
                             ybottom = 0,ytop = ytop)

               })

  # add mark gene rank
  if(!is.null(markGene)){
    circos.labels(sectors = mgene_df$id,x = mgene_df$x,
                  labels = mgene_df$symbol,
                  cex = markGeneSize,
                  connection_height = mm_h(2.5))
  }

  # add gene rank track
  if(addGeneRank == TRUE){
    circos.track(sectors = geneSetID,ylim = range(gsdata$geneList),
                 panel.fun = function(x, y){
                   # gene rank col layer
                   tmpht <- subset(gsdata,id == CELL_META$sector.index)

                   circos.lines(x = tmpht$x,y = tmpht$geneList,
                                type = "h",baseline = 0,col = "grey")

                   circos.segments(x0 = 0,x1 = max(gsdata1$x),y0 = 0,y1 = 0,
                                   lty = "solid",col = "black")

                   # this will run much slowly
                   # circos.barplot(value = tmpht$geneList,pos = tmpht$x,
                   #                bar_width = 1,border = "grey",col = "grey")
                 })
  }


}
