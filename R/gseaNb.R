#' @title gseaNb
#' @name gseaNb
#' @author Jun Zhang
#' @param object GSEA enrich results.
#' @param subPlot which plot to show, 1/2/3, default is 3.
#' @param lineSize curve line size. default is 0.8.
#' @param geneSetID which pathway name to plot.
#' @param rmSegment whether to remove segment on the curve plot, default is FALSE.
#' @param termWidth the width or the term name, defalut is 40.
#' @param segCol segment color on the curves, defalut is "red".
#' @param addGene whether add gene name on the curve, defalut is FALSE.
#' @param geneCol gene name label color, defalut is NULL.
#' @param arrowAngle arrow angle, defalut is 20.
#' @param arrowLength arrow line length, defalut is 0.2.
#' @param arrowEnd arrow end, defalut is "last".
#' @param arrowType arrow type, defalut is "closed".
#' @param curveCol curve color, defalut is c("#76BA99", "#EB4747", "#996699").
#' @param htCol heatmap color, defalut is c("#08519C", "#A50F15").
#' @param rankCol gene rank fill color, defalut is c("#08519C", "white", "#A50F15").
#' @param rankSeq gene rank plot X axis breaks, defalt is 5000.
#' @param htHeight the relative height when "subplot = 2" to the vertical line plot, defalut is 0.3.
#' @param force the gene label force, refer to geom_text_repel function, defalut is 20.
#' @param max.overlaps refer to geom_text_repel function, defalut is 50.
#' @param geneSize gene label text size, defalut is 4.
#' @param newGsea whether show new style of plot, defalut is FALSE.
#' @param addPoint new style plot with point layer, defalut is TRUE.
#' @param newCurveCol new style plot curve color, defalut is c("#336699", "white", "#993399").
#' @param newHtCol new style plot heatmap color, defalut is c("#336699", "white", "#993399").
#' @param rmHt whether remove new style plot heatmap, defalut is FALSE.
#' @param addPval whether add pvalue and NES, defalut is FALSE.
#' @param pvalX set pvalue label x position, defalut is 0.9.
#' @param pvalY set pvalue label y position, defalut is 0.9.
#' @param pvalSize set pvalue label text size, defalut is 4.
#' @param pCol pvalue label color, defalut is "grey30".
#' @param pHjust pvalue label hjust, defalut is 1.
#'
#' @param rmPrefix whether remove GO term prefix like "GOBP/KEGG/CC/MF_*", defalut is TRUE.
#' @param nesDigit the NES score digits retained, defalut is 2.
#' @param pDigit the pvalue and pajust value digits retained, defalut is 2.
#' @param markTopgene whether add top n genes on plot, defalut is FALSE.
#' @param topGeneN the number of genes to be marked on plot, defalut is 5.
#' @param kegg whether input is gseKEGG object, defalut is FALSE.
#' @param legend.position the legend position, defalut is "right".
#'
#' @importFrom ggplot2 aes_
#' @import DOSE
#' @import RColorBrewer
#' @import ggpp
#' @return ggplot2 object
#' @export
#'
#' @examples
#'\donttest{# load data
#'test_data <- system.file("extdata", "gseaRes.RDS", package = "GseaVis")
#'gseaRes <- readRDS(test_data)
#'
#'# all plot
#'gseaNb(object = gseaRes,
#'       geneSetID = 'GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS')
#'}

globalVariables(c(".", "ID", "aes_", "gene_name","gseaRes", "position","x","y"))

# define function
gseaNb <- function(object = NULL,
                   subPlot = 3,
                   lineSize = 0.8,
                   geneSetID = NULL,
                   rmSegment = FALSE,
                   termWidth = 40,
                   segCol = "red",
                   addGene = NULL,
                   geneCol = NULL,
                   arrowAngle = 20,
                   arrowLength = 0.2,
                   arrowEnd = "first",
                   arrowType = "closed",
                   curveCol = c("#76BA99", "#EB4747","#996699"),
                   htCol = c("#08519C", "#A50F15"),
                   rankCol = c("#08519C", "white", "#A50F15"),
                   rankSeq = 5000,
                   htHeight = 0.3,
                   force = 20,
                   max.overlaps = 50,
                   geneSize = 4,
                   newGsea = FALSE,
                   addPoint = TRUE,
                   newCurveCol = c("#336699", "white", "#993399"),
                   newHtCol = c("#336699", "white", "#993399"),
                   rmHt = FALSE,
                   addPval = FALSE,
                   pvalX = 0.9,
                   pvalY = 0.9,
                   pvalSize = 4,
                   pCol = "grey30",
                   pHjust = 1,
                   rmPrefix = TRUE,
                   nesDigit = 2,
                   pDigit = 2,
                   markTopgene = FALSE,
                   topGeneN = 5,
                   kegg = FALSE,
                   legend.position = "right") {
  # get dat
  gsdata <- purrr::map_df(geneSetID,function(setid){
    gsInfo(object,geneSetID = setid)
  })

  if(kegg == FALSE){
    # filter in pathway gene
    # gsdata1 <- gsdata %>%
    #   dplyr::mutate("gene_name" = names(object@geneList)) %>%
    #   dplyr::filter(position == 1)

    gsdata1 <- purrr::map_df(unique(gsdata$Description),function(setid){
      tmp <- gsdata %>%
        dplyr::filter(Description == setid) %>%
        dplyr::mutate("gene_name" = names(object@geneList)) %>%
        dplyr::filter(position == 1)
    })
  }else{
    # filter in pathway gene
    gene2Symbol <- object@gene2Symbol %>% data.frame()

    # gsdata1 <- gsdata %>%
    #   dplyr::mutate("gene_name" = gene2Symbol$.) %>%
    #   dplyr::filter(position == 1)

    gsdata1 <- purrr::map_df(unique(gsdata$Description),function(setid){
      tmp <- gsdata %>%
        dplyr::filter(Description == setid) %>%
        dplyr::mutate("gene_name" = gene2Symbol$.) %>%
        dplyr::filter(position == 1)
    })
  }

  # to dataframe
  data_ga <- data.frame(object) %>%
    dplyr::filter(ID %in% geneSetID)
  data_ga <- data_ga[unique(gsdata$Description),]

  ################################################
  # nice title
  niceTit <- purrr::map_chr(unique(gsdata$Description),function(x){
    tit <- unlist(strsplit(x, split = "_"))

    if(length(tit) == 1){
      niceTit <-
        paste(stringr::str_to_title(tit[1:length(tit)]), collapse = " ") %>%
        stringr::str_wrap(., width = termWidth)
    }else{
      if(rmPrefix == TRUE){
        niceTit <-
          paste(stringr::str_to_title(tit[2:length(tit)]), collapse = " ") %>%
          stringr::str_wrap(., width = termWidth)
      }else{
        niceTit <-
          paste(stringr::str_to_title(tit[1:length(tit)]), collapse = " ") %>%
          stringr::str_wrap(., width = termWidth)
      }
    }
  })

  # rename term id
  if(length(geneSetID) != 1){
    ledend.t <- niceTit
    niceTit <- ''
  }

  ################################################
  if(length(geneSetID) == 1){
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~runningScore),
                               size = lineSize)
    line.col <- ggplot2::scale_color_gradient(low = curveCol[1], high = curveCol[2])
    legend.position = "none"
  }else{
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~Description),
                               size = lineSize)
    line.col <- ggplot2::scale_color_manual(values = curveCol,
                                            labels = ledend.t,
                                            name = 'Term Name')
    legend.position = legend.position
  }

  # plot
  pcurve <-
    ggplot2::ggplot(gsdata,ggplot2::aes_(x = ~x, y = ~runningScore)) +
    line +
    line.col +
    ggplot2::geom_hline( yintercept = 0,
                         size = lineSize,
                         color = "black",
                         lty = "dashed") +
    # ggplot2::geom_line(size = lineSize) +
    # geom_segment(data = gsdata1,aes(xend = x,yend = 0)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(legend.position = legend.position,
                   legend.box.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.margin = ggplot2::margin(t = .2,
                                                 r = .2,
                                                 b = 0,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::ylab("Running Enrichment Score") +
    ggplot2::ggtitle(niceTit)

  ###########################################
  # calculate midpoint
  midpoint <- sum(range(gsdata$runningScore)) / 2

  pnew <-
    ggplot2::ggplot(gsdata,
                    ggplot2::aes_(x = ~x, y = ~runningScore, color = ~runningScore)) +
    ggplot2::geom_hline(yintercept = 0,
                        size = lineSize,
                        color = "black",
                        lty = "dashed") +
    ggplot2::geom_line(size = lineSize) +
    ggplot2::geom_segment(data = gsdata1, ggplot2::aes_(xend = ~x, yend = 0)) +
    ggplot2::theme_bw(base_size = 14) +
    # scale_color_gradient(low = '#336699',high = '#993399') +
    ggplot2::scale_color_gradient2(low = newCurveCol[1],
                                   mid = newCurveCol[2],
                                   high = newCurveCol[3],
                                   midpoint = midpoint) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.margin = ggplot2::margin(t = .2,
                                                 r = .2,
                                                 b = 0,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::ylab("Running Enrichment Score") +
    ggplot2::ggtitle(niceTit)

  # whether add point
  if (addPoint == TRUE) {
    panother <- pnew +
      ggplot2::geom_point()
  } else {
    panother <- pnew
  }

  ###########################################
  if (newGsea == FALSE) {
    pcurveRes <- pcurve
  } else {
    pcurveRes <- panother
  }

  ###########################################
  # add gene name
  if (is.null(addGene)) {
    plabel <- pcurveRes
  } else {
    # whether mark top genes
    if(markTopgene == TRUE){
      # add gene name
      geneLabel <- gsdata1 %>% dplyr::arrange(x) %>% dplyr::slice_head(n = topGeneN)
    }else{
      # add gene name
      geneLabel <- gsdata1 %>% dplyr::filter(gene_name %in% addGene)
    }

    # add gene on plot
    if (nrow(geneLabel) == 0) {
      print("Your gene is not in this pathway! Please choose again!")
    } else {
      if (rmSegment == TRUE) {
        if (is.null(geneCol)) {
          plabel <- pcurveRes +
            ggrepel::geom_text_repel(
              data = geneLabel,
              ggplot2::aes_(label = ~gene_name),
              force = force,
              max.overlaps = max.overlaps,
              # nudge_y = 0.2,
              size = geneSize,
              fontface = "italic",
              arrow = ggplot2::arrow(
                angle = arrowAngle,
                length = ggplot2::unit(arrowLength, "cm"),
                ends = arrowEnd,
                type = arrowType
              )
            )
        } else {
          plabel <- pcurveRes +
            ggrepel::geom_text_repel(
              data = geneLabel,
              ggplot2::aes_(label = ~gene_name),
              force = force,
              max.overlaps = max.overlaps,
              # nudge_y = 0.2,
              size = geneSize,
              fontface = "italic",
              color = geneCol,
              arrow = ggplot2::arrow(
                angle = arrowAngle,
                length = ggplot2::unit(arrowLength, "cm"),
                ends = arrowEnd,
                type = arrowType
              )
            )
        }
      } else {
        if (is.null(geneCol)) {
          plabel <- pcurveRes +
            ggplot2::geom_segment(
              data = geneLabel,
              ggplot2::aes_(xend = ~x, yend = 0),
              color = segCol
            ) +
            ggrepel::geom_text_repel(
              data = geneLabel,
              ggplot2::aes_(label = ~gene_name),
              force = force,
              max.overlaps = max.overlaps,
              # nudge_y = 0.2,
              size = geneSize,
              fontface = "italic",
              arrow = ggplot2::arrow(
                angle = arrowAngle,
                length = ggplot2::unit(arrowLength, "cm"),
                ends = arrowEnd,
                type = arrowType
              )
            )
        } else {
          plabel <- pcurveRes +
            ggplot2::geom_segment(
              data = geneLabel,
              ggplot2::aes_(xend = ~x, yend = 0),
              color = segCol
            ) +
            ggrepel::geom_text_repel(
              data = geneLabel,
              ggplot2::aes_(label = ~gene_name),
              force = force,
              max.overlaps = max.overlaps,
              # nudge_y = 0.2,
              size = geneSize,
              fontface = "italic",
              color = geneCol,
              arrow = ggplot2::arrow(
                angle = arrowAngle,
                length = ggplot2::unit(arrowLength, "cm"),
                ends = arrowEnd,
                type = arrowType
              )
            )
        }
      }
    }
  }

  ###########################################
  # add NES Pvalue
  if (addPval == TRUE) {
    pLabel <- paste0(
      "NES: ",
      round(data_ga$NES, digits = nesDigit),
      "\n",
      "Pvalue: ",
      # round(data_ga$pvalue, digits = pDigit),
      ifelse(data_ga$pvalue < 0.001,"< 0.001",round(data_ga$pvalue, digits = pDigit)),
      "\n",
      "Ajusted Pvalue: ",
      # round(data_ga$p.adjust, digits = pDigit),
      ifelse(data_ga$p.adjust < 0.001,"< 0.001",round(data_ga$p.adjust, digits = pDigit)),
      "\n",
      sep = " "
    )

    # define pvalue position
    # if (data_ga$NES > 0) {
    #   px <- pvalX * nrow(gsdata)
    #   py <-
    #     pvalY * sum(abs(range(gsdata$runningScore))) + min(gsdata$runningScore)
    # } else {
    #   px <- pvalX * nrow(gsdata)
    #   py <-
    #     pvalY * sum(abs(range(gsdata$runningScore))) + min(gsdata$runningScore)
    # }

    px <- pvalX * nrow(gsdata[which(gsdata$Description == geneSetID[1]),])
    py <-
      pvalY * sum(abs(range(gsdata$runningScore))) + min(gsdata$runningScore)

    # add pvlaue label

    if(length(geneSetID) == 1){
      pLabelOut <-
        plabel +
        ggplot2::annotate(geom = "text",
                          x = px,
                          y = py,
                          label = pLabel,
                          size = pvalSize,
                          color = pCol,
                          fontface = "italic",
                          hjust = pHjust)
    }else{
      mytable <- tibble::tibble(x = px, y = py,
                                table = list(tibble::tibble('NES' = round(data_ga$NES, digits = nesDigit),
                                                            'Pvalue' = ifelse(data_ga$pvalue < 0.001,"< 0.001",round(data_ga$pvalue, digits = pDigit)),
                                                            'Ajusted Pvalue' = ifelse(data_ga$p.adjust < 0.001,"< 0.001",round(data_ga$p.adjust, digits = pDigit)))))


      pLabelOut <-plabel +
        ggpp::geom_table(data = mytable, ggplot2::aes(x, y, label = table))
    }

  } else {
    pLabelOut <- plabel
  }

  ################################################
  if(length(geneSetID) == 1){
    line.col <- ggplot2::scale_color_manual(values = 'black')
  }

  pseg <-
    ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore,color = ~Description)) +
    ggplot2::geom_segment(data = gsdata1,
                          ggplot2::aes_(x = ~x,
                                        xend = ~x,
                                        y = 0,
                                        yend = 1),
                          # color = "black",
                          show.legend = F) +
    line.col +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(0.1,'cm'),
                   plot.margin = ggplot2::margin(t = 0,
                                                 r = .2,
                                                 b = .2,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::xlab("Rank in Ordered Dataset") +
    ggplot2::facet_wrap(~Description,ncol = 1)

  ################################################
  # v <- seq(1, sum(gsdata$position), length.out = 9)
  # inv <- findInterval(rev(cumsum(gsdata$position)), v)
  # if (min(inv) == 0) {
  #   inv <- inv + 1
  # }
  #
  # # col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
  # # new color
  # color <- grDevices::colorRampPalette(c(htCol[1], "white", htCol[2]))(10)
  #
  # ymin <- 0
  # yy <- htHeight
  # xmin <- which(!duplicated(inv))
  # xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  # d <- data.frame(ymin = ymin,
  #                 ymax = yy,
  #                 xmin = xmin,
  #                 xmax = xmax,
  #                 col = color[unique(inv)])

  d <- purrr::map_df(unique(gsdata$Description),function(setid){
    tmp <- gsdata %>%
      dplyr::filter(Description == setid)

    v <- seq(1, sum(tmp$position), length.out = 9)
    inv <- findInterval(rev(cumsum(tmp$position)), v)
    if (min(inv) == 0) {
      inv <- inv + 1
    }

    # col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    # new color
    color <- grDevices::colorRampPalette(c(htCol[1], "white", htCol[2]))(10)

    ymin <- 0
    yy <- htHeight
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin,
                    ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = color[unique(inv)],
                    Description = setid)
  })

  pseg_ht <-
    pseg + ggplot2::geom_rect(
      ggplot2::aes_(xmin = ~xmin,
                    xmax = ~xmax,
                    ymin = ~ymin,
                    ymax = ~ymax,
                    fill = ~ I(col)),
      data = d,
      alpha = 0.8,
      inherit.aes = FALSE)

  ################################################
  # add gene rank
  pseg_ht1 <-
    pseg_ht + ggplot2::xlab("") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = -.1,
                                    r = .2,
                                    b = 0,
                                    l = .2,
                                    unit = "cm"))

  prank <-
    ggplot2::ggplot(gsdata[which(gsdata$Description == unique(gsdata$Description)[1]),],
                    ggplot2::aes_(x = ~x, y = ~geneList)) +
    # geom_col(width = 1,fill = 'grey80',color = NA) +
    ggplot2::geom_col(
      ggplot2::aes_(fill = ~geneList),
      width = 1,
      color = NA,
      show.legend = F) +
    ggplot2::scale_fill_gradient2(low = rankCol[1],
                                  mid = rankCol[2],
                                  high = rankCol[3],
                                  midpoint = 0) +
    ggplot2::geom_hline(yintercept = 0,
                        size = 0.8,
                        color = "black",
                        lty = "dashed") +
    ggplot2::scale_x_continuous(breaks = seq(0, nrow(gsdata), rankSeq)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = -.1,
                                                 r = .2,
                                                 b = .2,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::ylab("Ranked List") +
    ggplot2::xlab("Rank in Ordered Dataset")

  ###########################################
  # new color
  d <- purrr::map_df(unique(d$Description),function(x){
    tmp <- d %>%
      dplyr::filter(Description == x)

    # add color
    htcolor <- grDevices::colorRampPalette(newHtCol)(nrow(tmp))

    tmp <- tmp %>%
      dplyr::mutate(htcol = htcolor)
  })

  # heatmap
  ht <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore)) +
    ggplot2::geom_rect(ggplot2::aes_(xmin = ~xmin,
                                     xmax = ~xmax,
                                     ymin = ~ymin,
                                     ymax = ~ymax,
                                     fill = ~ I(htcol)),
                       data = d,
                       alpha = 0.8,
                       inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 0,
                                                 r = .2,
                                                 b = .2,
                                                 l = .2,
                                                 unit = "cm"))

  ###########################################
  if (newGsea == FALSE) {
    # subplots to show
    if (subPlot == 1) {
      pres <- pLabelOut
      # return(pres)
    } else if (subPlot == 2) {
      # combine
      if (rmHt == FALSE){
        pres <- aplot::plot_list(
          gglist = list(pLabelOut, pseg_ht),
          ncol = 1,
          heights = c(0.8, 0.2)
        )
      }else{
        pres <- aplot::plot_list(
          gglist = list(pLabelOut, pseg),
          ncol = 1,
          heights = c(0.8, 0.2)
        )
      }

      # return(pres)
    } else if (subPlot == 3) {
      # combine
      if (rmHt == FALSE){
        pres <-
          aplot::plot_list(
            gglist = list(pLabelOut, pseg_ht1, prank),
            ncol = 1,
            heights = c(0.5, 0.2, 0.3)
          )
      }else{
        pres <-
          aplot::plot_list(
            gglist = list(pLabelOut, pseg, prank),
            ncol = 1,
            heights = c(0.5, 0.2, 0.3)
          )
      }

      # return(pres)
    } else {
      print("Please give 1/2/3 parameters!")
    }
  } else {
    # combine
    if (rmHt == FALSE) {
      pres <- aplot::plot_list(
        gglist = list(pLabelOut, ht),
        ncol = 1,
        heights = c(0.9, 0.1)
      )
    } else {
      pres <- pLabelOut
    }
  }

  # output
  return(pres)
}
