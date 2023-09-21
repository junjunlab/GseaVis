globalVariables(c("val","label","hjust"))
#' Perform GSEA Multi-group Plotting
#'
#' This function generates a multi-group plot for Gene Set Enrichment Analysis
#' (GSEA) results.
#'
#' @author Jun Zhang
#'
#' @param gsea_list A list of GSEA results for multiple experiments.
#' @param geneSetID The ID of the gene set to be visualized.
#' @param exp_name Names of the experiments corresponding to the GSEA results.
#' @param addPval Logical, indicating whether to add NES (Normalized Enrichment Score)
#' and p-value labels to the plot. Default is FALSE.
#' @param curve.col A vector of colors for the curves representing different
#' experiments. If NULL, random colors are assigned.
#' @param kegg Logical, indicating whether the gene set is a KEGG pathway.
#' Default is FALSE.
#' @param lineSize The size of the lines in the enrichment score curve plot.
#' Default is 1.
#' @param base_size The base font size for the plot. Default is 12.
#' @param nesDigit The number of digits to round NES values to. Default is 2.
#' @param pDigit The number of digits to round p-values to. Default is 2.
#' @param pvalX The X-coordinate for placing p-value labels on the plot. Default is 0.9.
#' @param pvalY The Y-coordinate for placing p-value labels on the plot. Default is 0.9.
#' @param rect.bm.col A vector of colors for the bottom rectangle representing
#' up-regulated and down-regulated genes. Default colors c("#CC3333","white","#003366")
#' are provided.
#' @param subplot.heights Heights of subplots in the multi-group plot. Default
#' values c(0.4,0.2,0.08) are provided.
#' @param legend.position The position of the legend in the plot. Default is (0.85, 0.85).
#' @param rect.bm.label Labels for the bottom rectangle, specifying "Up regulated"
#' and "Down regulated". Default labels are provided.
#' @param breaks.n The number of X axis breaks. Default is 6.
#'
#' @return Returns a multi-panel plot for GSEA results.
#'
#' @examples
#' \dontrun{
#' # Example Usage
#' result <- GSEAmultiGP(gsea_list = gsea_results,
#'                       geneSetID = "gene_set_1",
#'                       exp_name = c("Exp1", "Exp2"),
#'                       addPval = TRUE,
#'                       curve.col = c("red", "blue"),
#'                       pvalX = 0.9,
#'                       pvalY = 0.9,
#'                       rect.bm.col = c("#CC3333", "white", "#003366"),
#'                       subplot.heights = c(0.4, 0.2, 0.08),
#'                       legend.position = c(0.85, 0.85),
#'                       rect.bm.label = c("Up regulated", "Down regulated"))
#'                       }
#'
#' @import ggplot2
#' @import purrr
#' @import dplyr
#' @importFrom ggpp geom_table
#' @import circlize
#'
#' @seealso
#' \code{\link{gsInfo}}, \code{\link{ggplot2}}
#'
#' @export
GSEAmultiGP <- function(gsea_list = NULL,
                        geneSetID = NULL,
                        exp_name = NULL,
                        addPval = FALSE,
                        curve.col = NULL,
                        kegg = FALSE,
                        lineSize = 1,
                        base_size = 12,
                        nesDigit = 2,pDigit = 2,
                        pvalX = 0.9,pvalY = 0.9,
                        rect.bm.col = c("#CC3333","white","#003366"),
                        subplot.heights = c(0.4,0.2,0.08),
                        legend.position = c(0.85,0.85),
                        rect.bm.label = c("Up regulated","Down regulated"),
                        breaks.n = 6){
  # ============================================================================
  # data process
  # ============================================================================
  if(length(geneSetID) > 1){
    message("Please give one geneSetID when multi.group set to TRUE!")
    break
  }

  # process data
  gsdata <- purrr::map_df(1:length(gsea_list),function(x){
    tmp <- gsInfo(gsea_list[[x]],geneSetID = geneSetID)
    # tmp$Description <- paste(exp_name[x],tmp$Description,sep = "_")
    # tmp$id <- paste(exp_name[x],geneSetID,sep = "_")
    tmp$id <- exp_name[x]

    return(tmp)
  })

  gsdata1 <- purrr::map_df(1:length(gsea_list),function(x){
    tmp <- gsInfo(gsea_list[[x]],geneSetID = geneSetID)
    # tmp$Description <- paste(exp_name[x],tmp$Description,sep = "_")
    # tmp$id <- paste(exp_name[x],geneSetID,sep = "_")
    tmp$id <- exp_name[x]

    if(kegg == FALSE){
      # filter in pathway gene
      gn <- names(gsea_list[[x]]@geneList)
      tmp2 <- tmp %>%
        dplyr::mutate(gene_name = gn) %>%
        dplyr::filter(position == 1)
    }else{
      # filter in pathway gene
      gene2Symbol <- gsea_list[[x]]@gene2Symbol %>% data.frame()

      tmp2 <- tmp %>%
        dplyr::filter(Description == geneSetID) %>%
        dplyr::mutate("gene_name" = gene2Symbol$.) %>%
        dplyr::filter(position == 1)
    }

    return(tmp2)
  })


  data_ga <- purrr::map_df(1:length(gsea_list),function(x){
    # to dataframe
    data_ga <- data.frame(gsea_list[[x]]) %>%
      dplyr::filter(ID %in% geneSetID)
    # data_ga$Description <- paste(exp_name[x],data_ga$Description,sep = "_")
    # data_ga$ID <- paste(exp_name[x],geneSetID,sep = "_")
    data_ga$ID <- exp_name[x]
    return(data_ga)
  })

  # ============================================================================
  # plot
  # ============================================================================
  # colors
  if(is.null(curve.col)){
    line.col <- circlize::rand_color(n = length(exp_name))
    names(line.col) <- exp_name
  }else{
    line.col <- curve.col
    names(line.col) <- exp_name
  }

  # curve plot
  pcurve <-
    ggplot2::ggplot(gsdata,ggplot2::aes_(x = ~x, y = ~runningScore)) +
    ggplot2::geom_line(ggplot2::aes(color = id)) +
    ggplot2::scale_color_manual(values = line.col,name = "") +
    ggplot2::geom_hline(yintercept = 0,
                        size = lineSize,
                        color = "black",
                        lty = "dashed") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::scale_x_continuous(expand = c(0, 0),limits = c(0,max(gsdata$x)),
                                breaks = seq(0,max(gsdata$x),length.out = breaks.n),
                                labels = seq(0,max(gsdata$x),length.out = breaks.n)) +
    ggplot2::theme(legend.position = legend.position,
                   legend.box.background = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(hjust = 0.5,face = "bold.italic"),
                   plot.title = ggplot2::element_text(hjust = 0.5,face = "bold.italic"),
                   panel.grid = ggplot2::element_blank(),
                   axis.text = element_text(colour = "black"),
                   axis.line.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.margin = ggplot2::margin(t = .2,
                                                 r = .2,
                                                 b = 0,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::ylab("Enrichment Score") +
    ggplot2::ggtitle(unique(data_ga$Description))


  # add NES Pvalue
  if (addPval == TRUE) {
    # add pvlaue label
    mytable <- tibble::tibble(x = pvalX * max(gsdata$x),
                              y = pvalY * sum(abs(range(gsdata$runningScore))) + min(gsdata$runningScore),
                              table = list(tibble::tibble("Exp" = data_ga$ID,
                                                          'NES' = round(data_ga$NES, digits = nesDigit),
                                                          'Pvalue' = ifelse(data_ga$pvalue < 0.001,"< 0.001",round(data_ga$pvalue, digits = pDigit)),
                                                          'Adjusted Pvalue' = ifelse(data_ga$p.adjust < 0.001,"< 0.001",round(data_ga$p.adjust, digits = pDigit)))))

    pLabelOut <- pcurve + ggpp::geom_table(data = mytable, ggplot2::aes(x, y, label = table))
  }else{
    pLabelOut <- pcurve
  }

  # segment plot
  pseg <-
    ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore,color = ~id)) +
    ggplot2::geom_segment(data = gsdata1,
                          ggplot2::aes_(x = ~x,xend = ~x,y = 0,yend = 1),
                          show.legend = F) +
    ggplot2::scale_color_manual(values = line.col,name = "") +
    ggplot2::scale_x_continuous(expand = c(0, 0),name = "") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::facet_wrap(~id,ncol = 1) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_blank(),
                   panel.spacing.y = ggplot2::unit(1,"mm"),
                   plot.margin = ggplot2::margin(t = 0,
                                                 r = .2,
                                                 b = 0,
                                                 l = .2,
                                                 unit = "cm"))


  # prepare bottom rect
  rect_b <- data.frame(x = 1:max(gsdata$x),
                       val = seq(-2,2,length.out = length(1:max(gsdata$x))))
  label_b <- data.frame(x = c(0,max(gsdata$x)),
                        y = 1,
                        label = rect.bm.label,
                        hjust = c(0,1))

  # rect plot
  rect_bp <-
    ggplot2::ggplot(rect_b) +
    ggplot2::geom_tile(ggplot2::aes(x = x,y = 1,fill = val),show.legend = F) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0,xmax = max(gsdata$x),ymin = 0.5,ymax = 1.5),
                       fill = NA,color = "black") +
    ggplot2::scale_fill_gradient2(low = rect.bm.col[1],mid = rect.bm.col[2],high = rect.bm.col[3],
                                  midpoint = 0) +
    ggplot2::theme(plot.background = ggplot2::element_rect(colour = "black")) +
    ggplot2::geom_label(data = label_b,
                        ggplot2::aes(x = x,y = y,label = label,hjust = hjust),
                        fontface = "bold.italic") +
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
                                                 b = 0.5,
                                                 l = .2,
                                                 unit = "cm")) +
    ggplot2::xlab("Rank in Ordered Dataset")



  # combine plot
  aplot::plot_list(gglist = list(pLabelOut,pseg,rect_bp),
                   ncol = 1,
                   heights = subplot.heights)

}
