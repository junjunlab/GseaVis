#' @title gseaNb
#' @name gseaNb
#' @author Jun Zhang
#' @param object GSEA enrich results.
#' @param filePath filePath the path of the GSEA software enrichment
#'  outputs or "readGseaFile" object, defalut is NULL.
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
#' @param whether add target gene expression heatmap, defalut is FALSE.
#' @param exp the expression matrix,tpm/fpkm/rpkm format, defalut is NULL.
#' @param scale.exp whether scale the expression matrix, defalut is TRUE.
#' @param sample.order the expression matrix sample orders, defalut is NULL.
#' @param exp.col the expression colors, defalut is c('blue','white','red').
#' @param ht.legend whether show the heatmap legend, defalut is TRUE.
#' @param ght.relHight the relative height to the main plot, defalut is 0.4.
#' @param ght.geneText.size the gene lable text size, defalut is 6.
#' @param ght.facet whether facet expression heatmap, defalut is FALSE.
#' @param ght.facet.scale the facet plot scale argumrnt, defalut is "free".
#' @param termID.order the facet term ID orders, defalut is NULL.
#'
#' @param rank.gene add your gene label on rank plot, defalut is NULL.
#' @param rank.gene.nudgey the gene label nudge y on rank plot, defalut is 2.
#' @param rm.newGsea.ticks whether remove right axis when you plot multiple terms with newGsea plot, defalut is TRUE.
#' @param pFill the pvalue table fill color when you plot multiple terms with newGsea plot, defalut is transparent.
#' @param base_size the plot theme font size, defalut is 12.
#' @param ncol the columns for newGSEA plot with multiple terms, defalut is 1
#'
#' @importFrom ggplot2 aes_
#' @import DOSE
#' @import RColorBrewer
#' @import ggpp
#' @return ggplot2 object
#' @export
#'
#' @examples
#'# load data
#'test_data <- system.file("extdata", "gseaRes.RDS", package = "GseaVis")
#'gseaRes <- readRDS(test_data)
#'
#'# all plot
#'gseaNb(object = gseaRes,
#'       geneSetID = 'GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS',
#'       subPlot = 2)

globalVariables(c(".", "ID", "aes_", "gene_name","gseaRes",
                  "position","x","y","value","variable","logfc","nudge_y","vjust",
                  "pLabel","px","py","id"))

# define function
gseaNb <- function(object = NULL,
                   filePath = NULL,
                   subPlot = 3,
                   lineSize = 1,
                   base_size = 12,
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
                   addPoint = FALSE,
                   ncol = 1,
                   newCurveCol = c("#336699", "grey80", "#993399"),
                   newHtCol = c("#993399", "white", "#336699"),
                   rm.newGsea.ticks = TRUE,
                   rmHt = FALSE,
                   addPval = FALSE,
                   pvalX = 0.9,
                   pvalY = 0.9,
                   pvalSize = 4,
                   pCol = "grey0",
                   pFill = "transparent",
                   pHjust = 1,
                   rmPrefix = TRUE,
                   nesDigit = 2,
                   pDigit = 2,
                   markTopgene = FALSE,
                   topGeneN = 5,
                   kegg = FALSE,
                   legend.position = "right",
                   add.geneExpHt = FALSE,
                   exp = NULL,
                   scale.exp = TRUE,
                   sample.order = NULL,
                   exp.col = c('blue','white','red'),
                   ht.legend = TRUE,
                   ght.relHight = 0.4,
                   ght.geneText.size = 6,
                   ght.facet = FALSE,
                   ght.facet.scale = "free",
                   termID.order = NULL,
                   rank.gene = NULL,
                   rank.gene.nudgey = 2) {
  ##################################################################################
  # prepare data for plot
  ##################################################################################
  if(is.null(filePath)){
    gsdata <- purrr::map_df(geneSetID,function(setid){
      gsInfo(object,geneSetID = setid) %>%
        dplyr::mutate(id = setid)
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
    data_ga <- data_ga[unique(gsdata$id),]
  }else{
    # load GSEA software outputs
    if(is.character(filePath)){
      gsea.res <- readGseaFile(filePath = filePath)
    }else{
      gsea.res <- filePath
    }

    # get rank data
    gsdata <- purrr::map_df(geneSetID,function(setid){
      gsInfoNew(geneList = gsea.res$glist,
                geneSet = gsea.res$gset,
                geneSetID = setid) %>%
        dplyr::mutate(id = setid)
    })

    # filter terms gene
    gsdata1 <- purrr::map_df(unique(gsdata$Description),function(setid){
      tmp <- gsdata %>%
        dplyr::filter(Description == setid) %>%
        dplyr::mutate("gene_name" = names(gsea.res$glist)) %>%
        dplyr::filter(position == 1)
    })

    # to dataframe
    data_ga <- gsea.res$meta %>%
      dplyr::filter(ID %in% geneSetID) %>%
      tibble::column_to_rownames("ID") %>%
      dplyr::mutate(ID = rownames(.))
    data_ga <- data_ga[unique(gsdata$id),]
  }


  ##################################################################################
  # plot
  ##################################################################################
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
  }else{
    ledend.t <- niceTit
  }

  ################################################
  if(length(geneSetID) == 1){
    if(length(curveCol) == 1){
      line <- ggplot2::geom_line(color = curveCol,size = lineSize)
      line.col <- NULL
    }else{
      line <- ggplot2::geom_line(ggplot2::aes_(color = ~runningScore),
                                 size = lineSize)
      line.col <- ggplot2::scale_color_gradient(low = curveCol[1], high = curveCol[2])
    }

    legend.position = "none"
  }else{
    # assign colors
    if(newGsea == FALSE){
      mulcol <- curveCol

      # if(gsdata$id[1] == gsdata$Description[1]){
      #   names(mulcol) <- geneSetID
      # }else{
      #   names(mulcol) <- unique(gsdata$Description)
      # }

      names(mulcol) <- geneSetID

      line.col <- ggplot2::scale_color_manual(values = mulcol,
                                              labels = ledend.t,
                                              name = 'Term Name')
    }else{
      line.col <- ggplot2::scale_color_brewer()
    }

    # layers
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~id),
                               size = lineSize)

    legend.position = legend.position
  }

  # order
  # gsdata$Description <- factor(gsdata$Description,levels = geneSetID)
  if(gsdata$id[1] == gsdata$Description[1]){
    gsdata$id <- factor(gsdata$id,levels = geneSetID)
  }else{
    gsdata$id <- factor(gsdata$id,levels = data_ga$ID)
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
    ggplot2::theme_bw(base_size = base_size) +
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

  # order
  # gsdata1$Description <- factor(gsdata1$Description,levels = geneSetID)
  if(gsdata$id[1] == gsdata$Description[1]){
    gsdata1$id <- factor(gsdata1$id,levels = geneSetID)
  }else{
    gsdata1$id <- factor(gsdata1$id,levels = data_ga$ID)
  }

  # facet labels
  facet.label <- ledend.t
  names(facet.label) <- geneSetID

  # plot
  if(length(newCurveCol) > 3){
    gseaNewCol <- ggplot2::scale_color_gradientn(colors = newCurveCol)
  }else{
    gseaNewCol <- ggplot2::scale_color_gradient2(low = newCurveCol[1],
                                                 mid = newCurveCol[2],
                                                 high = newCurveCol[3],
                                                 midpoint = midpoint)
  }

  pnew <-
    ggplot2::ggplot(gsdata,
                    ggplot2::aes_(x = ~x, y = ~runningScore, color = ~runningScore)) +
    ggplot2::geom_line(size = lineSize) +
    ggplot2::geom_hline(yintercept = 0,
                        size = lineSize,
                        color = "black",
                        lty = "dashed") +
    ggplot2::geom_segment(data = gsdata1, ggplot2::aes_(xend = ~x, yend = 0)) +
    ggplot2::theme_bw(base_size = base_size) +
    gseaNewCol +
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

  # facet plot if geneSetID more than 1
  if(length(geneSetID) > 1){
    pnew <- pnew +
      ggplot2::facet_wrap(~id,
                          ncol = ncol,
                          scales = "free_y",
                          strip.position = "left",
                          labeller = ggplot2::labeller(id = facet.label)) +
      ggplot2::theme(panel.spacing = ggplot2::unit(0,"mm"),
                     panel.grid = ggplot2::element_blank(),
                     strip.background = ggplot2::element_rect(fill = "grey95"),
                     strip.text.y.left = ggplot2::element_text(angle = 0,size = 12)) +
      ggplot2::scale_y_continuous(position = "right")
  }

  # whether remove strip aixs
  if(rm.newGsea.ticks == TRUE){
    pnew <- pnew +
      ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank())
  }

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
      # geneLabel <- gsdata1 %>% dplyr::arrange(x) %>% dplyr::slice_head(n = topGeneN)
      geneLabel <- purrr::map_df(unique(gsdata1$Description),function(dc){
        topg <- gsdata1 %>%
          dplyr::filter(Description == dc) %>%
          dplyr::arrange(x)

        # check NES score
        nes <- data_ga %>% dplyr::filter(Description == dc)

        if(nes$NES > 0){
          topg <- topg %>% dplyr::slice_head(n = topGeneN)
        }else{
          topg <- topg %>% dplyr::slice_tail(n = topGeneN)
        }

      })
    }else{
      # add gene name
      geneLabel <- gsdata1 %>% dplyr::filter(gene_name %in% addGene)
    }

    # add gene on plot
    if (nrow(geneLabel) == 0) {
      message("Your gene is not in this pathway! Please choose again!")
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
    if(length(geneSetID) == 1 | newGsea == FALSE){
      pLabel <- paste0(
        "NES: ",
        round(data_ga$NES, digits = nesDigit),
        "\n",
        "Pvalue: ",
        # round(data_ga$pvalue, digits = pDigit),
        ifelse(data_ga$pvalue < 0.001,"< 0.001",round(data_ga$pvalue, digits = pDigit)),
        "\n",
        "Adjusted Pvalue: ",
        # round(data_ga$p.adjust, digits = pDigit),
        ifelse(data_ga$p.adjust < 0.001,"< 0.001",round(data_ga$p.adjust, digits = pDigit)),
        sep = " "
      )

      # define pvalue position
      px <- pvalX * nrow(gsdata[which(gsdata$id == geneSetID[1]),])
      py <-
        pvalY * sum(abs(range(gsdata$runningScore))) + min(gsdata$runningScore)

    }else{
      purrr::map_df(unique(data_ga$ID),function(j){
        tmp <- data_ga %>% dplyr::filter(ID == j)
        tmpgsdata <- gsdata %>%
          dplyr::filter(id == j)

        # pvalue label
        pLabel <- paste0(
          "NES: ",
          round(tmp$NES, digits = nesDigit),
          "\n",
          "Pvalue: ",
          # round(data_ga$pvalue, digits = pDigit),
          ifelse(tmp$pvalue < 0.001,"< 0.001",round(tmp$pvalue, digits = pDigit)),
          "\n",
          "Adjusted Pvalue: ",
          # round(data_ga$p.adjust, digits = pDigit),
          ifelse(tmp$p.adjust < 0.001,"< 0.001",round(tmp$p.adjust, digits = pDigit)),
          sep = " "
        )

        # xy coordinate
        px <- pvalX * nrow(tmpgsdata)
        py <- pvalY * sum(abs(range(tmpgsdata$runningScore))) + min(tmpgsdata$runningScore)

        return(data.frame(pLabel = pLabel,px = px,py = py,id = j,
                          x = 0,runningScore = 0))
      }) -> ptable

      # order
      if(gsdata$id[1] == gsdata$Description[1]){
        ptable$id <- factor(ptable$id,levels = geneSetID)
      }else{
        ptable$id <- factor(ptable$id,levels = data_ga$ID)
      }
    }

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
      if(newGsea == FALSE){
        mytable <- tibble::tibble(x = px, y = py,
                                  table = list(tibble::tibble('NES' = round(data_ga$NES, digits = nesDigit),
                                                              'Pvalue' = ifelse(data_ga$pvalue < 0.001,"< 0.001",round(data_ga$pvalue, digits = pDigit)),
                                                              'Adjusted Pvalue' = ifelse(data_ga$p.adjust < 0.001,"< 0.001",round(data_ga$p.adjust, digits = pDigit)))))


        pLabelOut <- plabel +
          ggpp::geom_table(data = mytable, ggplot2::aes(px, py, label = table))
      }else{
        pLabelOut <- plabel +
          ggplot2::geom_label(data = ptable,ggplot2::aes(x = px,y = py,
                                                         label = pLabel),
                              color = pCol,
                              fill = pFill,
                              size = pvalSize,
                              fontface = "italic",
                              hjust = pHjust)
      }

    }

  } else {
    pLabelOut <- plabel
  }

  ################################################
  if(length(geneSetID) == 1){
    line.col <- ggplot2::scale_color_manual(values = 'black')
  }

  # bottomn space
  if(add.geneExpHt == TRUE){
    pseg.b = 0
  }else{
    pseg.b = 0.2
  }

  # order
  # gsdata1$Description <- factor(gsdata1$Description,levels = geneSetID)
  if(gsdata$id[1] == gsdata$Description[1]){
    gsdata1$id <- factor(gsdata1$id,levels = geneSetID)
  }else{
    gsdata1$id <- factor(gsdata1$id,levels = data_ga$ID)
  }

  # plot
  pseg <-
    ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore,color = ~id)) +
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
                                                 b = pseg.b,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::xlab("Rank in Ordered Dataset") +
    ggplot2::facet_wrap(~id,ncol = 1)

  if(subPlot > 2){
    pseg <-
      pseg +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }else{
    pseg <- pseg
  }

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

  # bottomn space
  if(add.geneExpHt == TRUE){
    prank.b = 0
  }else{
    prank.b = 0.2
  }

  # plot
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
                                                 b = prank.b,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::ylab("Ranked List") +
    ggplot2::xlab("Rank in Ordered Dataset")

  if(add.geneExpHt == TRUE){
    prank <- prank +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }else{
    prank <- prank
  }
  # =================================================
  # prank add gene label
  if(is.null(filePath)){
    if(kegg == FALSE){
      rank.g <- data.frame(logfc = object@geneList,
                           gene_name = names(object@geneList)) %>%
        dplyr::mutate(x = 1:length(object@geneList))
    }else{
      # check gene2Symbol in object
      if(length(object@gene2Symbol) > 0){
        rank.g <- data.frame(logfc = object@geneList,
                             gene_name = object@gene2Symbol) %>%
          dplyr::mutate(x = 1:length(object@geneList))
      }else{
        message("Please use readble gene symbols for your
                KEGG results('clusterProfiler::setReadable()')!")
      }
    }
  }else{
    rank.g <- data.frame(logfc = gsea.res$glist,
                         gene_name = names(gsea.res$glist)) %>%
      dplyr::mutate(x = 1:length(gsea.res$glist))
  }

  # filter rank genes
  if(!is.null(rank.gene)){
    target.rank.g <- rank.g %>%
      dplyr::filter(gene_name %in% rank.gene) %>%
      dplyr::mutate(vjust = ifelse(logfc > 0,"bottom","top"),
                    nudge_y = ifelse(logfc > 0,-rank.gene.nudgey,rank.gene.nudgey))

    # add rank gene label
    prank <-
      prank +
      ggrepel::geom_text_repel(data = target.rank.g,
                               ggplot2::aes(x = as.numeric(x),y = 0,
                                            label = gene_name,
                                            vjust = vjust,
                                            nudge_y = nudge_y),
                               max.overlaps = 200,direction = "x",angle = 90,
                               fontface = "italic",
                               size = geneSize)
  }

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


  # =========================================================
  # whether add gene expression heatmap
  if(add.geneExpHt == TRUE){
    # gene
    target.g <- purrr::map_df(data_ga$ID,function(x){
      tmp <- data_ga %>%
        dplyr::filter(ID == x)
      coregene <- unique(unlist(strsplit(tmp$core_enrichment,split = '\\/')))
      output <- data.frame(gene_name = coregene,
                           ID = x,
                           Description = tmp$Description) %>%
        dplyr::distinct(.,gene_name,.keep_all = TRUE)
    })

    # get gene position
    if(is.null(filePath)){
      gpos <- if(kegg == TRUE){
        match(target.g$gene_name,object@gene2Symbol)
      }else{
        match(target.g$gene_name,names(object@geneList))
      }
    }else{
      gpos <- match(target.g$gene_name,names(gsea.res$glist))
    }

    ginfo <- target.g %>%
      dplyr::mutate(gpos = gpos) %>%
      # data.frame(gpos = gpos,gene = target.g) %>%
      dplyr::arrange(gpos)

    # get expression data
    if(scale.exp == TRUE){
      gexp <- t(scale(t(exp[,2:ncol(exp)]),scale = TRUE,center = TRUE)) %>%
        data.frame()
      gexp$gene_name <- exp[,1]
    }else{
      gexp <- exp
      colnames(gexp)[1] <- "gene_name"
    }

    # wide2long
    exp.long <- gexp %>%
      dplyr::filter(gene_name %in% unique(ginfo$gene_name)) %>%
      dplyr::left_join(.,ginfo[,1:2],by = "gene_name") %>%
      reshape2::melt(.,id.vars = c("gene_name","ID"))

    # order
    exp.long$gene_name <- factor(exp.long$gene_name,levels = unique(ginfo$gene_name))
    if(!is.null(sample.order)){exp.long$variable <- factor(exp.long$variable,levels = sample.order)}
    if(!is.null(termID.order)){exp.long$ID <- factor(exp.long$ID,levels = termID.order)}

    # plot
    ght <-
      ggplot2::ggplot(exp.long) +
      ggplot2::geom_tile(ggplot2::aes(x = gene_name,y = variable,fill = value),
                         color = NA,
                         show.legend = ht.legend) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::coord_cartesian(expand = 0) +
      ggplot2::scale_fill_gradient2(low = exp.col[1],mid = exp.col[2],high = exp.col[3],
                                    midpoint = 0,name = 'Z-Score') +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5,hjust = 1,size = ght.geneText.size),
                     axis.text = ggplot2::element_text(color = "black"),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(t = -0.1,
                                                   r = .2,
                                                   b = .2,
                                                   l = .2,
                                                   unit = "cm")) +
      ggplot2::scale_y_discrete(position = "right") +
      ggplot2::xlab('') + ggplot2::ylab('')

    # whether facet heatmap
    if(ght.facet == TRUE){
      fght <- ght +
        ggplot2::facet_wrap(~ID,ncol = 1,scales = ght.facet.scale,
                            strip.position = "left") +
        ggplot2::theme(strip.background = ggplot2::element_rect(color = NA,fill = "grey90"),
                       strip.placement = "outside")
    }else{
      fght <- ght
    }
  }

  ###########################################
  if (newGsea == FALSE) {
    # subplots to show
    if (subPlot == 1) {
      pres <- pLabelOut
      # return(pres)
    } else if (subPlot == 2) {
      # combine
      if(length(geneSetID) == 1){
        seg.ht = 0.2
      }else{
        seg.ht = 0.5
      }

      if (rmHt == FALSE){
        pres <- aplot::plot_list(
          gglist = list(pLabelOut, pseg_ht),
          ncol = 1,
          heights = c(0.8, seg.ht)
        )
      }else{
        pres <- aplot::plot_list(
          gglist = list(pLabelOut, pseg),
          ncol = 1,
          heights = c(0.8, seg.ht)
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
      message("Please give 1/2/3 parameters!")
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

  # whether add gene expression heatmap
  if(add.geneExpHt == TRUE){
    pfinal <- aplot::plot_list(
      gglist = list(pres +
                      ggplot2::xlab(''), fght),
      ncol = 1,
      heights = c(1-ght.relHight, ght.relHight)
    )
  }else{
    pfinal <- pres
  }

  # output
  return(pfinal)
}


#' This is a test data for this package
#' test data describtion
#'
#' @name intergrated
#' @docType data
#' @author JunZhang
"intergrated"
