globalVariables(c("gene_ratio", 'next_node', "next_x", "node", "term_ratio"))

#' sankeyGoPlot
#'
#' @param goData GO data frame from Clusterprofiler.
#' @param topGenes the top genes to be shown for each term, default 5.
#' @param keep_all_gene Whether keep all genes, default "FALSE", `topGenes` will be useless if TRUE.
#' @param sankeyExpand the sankey plot expand for left and side, default c(0.5,1).
#' @param flow_fill The flow fill color, default `grey`.
#' @param nodeSize node text size for sankey plot, default 2.5.
#' @param nodeColor node fill color for sankey plot, default NULL.
#' @param goCol go color for go plot, default NULL.
#' @param xShift the shift on horizontal for go plot, default 0.05.
#' @param downShift the shift on vertical for go plot, default 4.5.
#' @param upShift the shift on vertical for go plot, default 0.25.
#' @param geom_layer Custom geom_layers, default NULL.
#' @param enrich_type The enrichment type, "go" or "gsea".
#'
#' @return ggplot obeject
#'
#' @import ggsankey
#' @importFrom cols4all c4a
#' @importFrom tidyr separate_longer_delim
#'
#' @export
sankeyGoPlot <- function(goData = NULL,
                         topGenes = 5,
                         keep_all_gene = FALSE,
                         sankeyExpand = c(0.5,1),
                         flow_fill = "grey",
                         nodeSize = 2.5,
                         nodeColor = NULL,
                         goCol = NULL,
                         xShift = 0.05,
                         downShift = 4.5,
                         upShift = 0.25,
                         geom_layer = NULL,
                         enrich_type = c("go","gsea")){
  enrich_type <- match.arg(enrich_type,c("go","gsea"))
  # ============================================================================
  # 1.prepare data
  # ============================================================================
  if(enrich_type == "gsea"){
    ego_df <- goData %>%
      dplyr::rename(geneID = core_enrichment) %>%
      dplyr::group_by(Description) %>%
      dplyr::arrange(pvalue)
  }else{
    ego_df <- goData %>%
      dplyr::group_by(Description) %>%
      dplyr::mutate(gene_ratio = eval(parse(text = GeneRatio))) %>%
      dplyr::arrange(pvalue)

    # calculate term_ratio(diffGenes / termGenes in pathway)
    ego_df$term_gene <- sapply(strsplit(ego_df$BgRatio,split = "\\/"),"[",1) %>%
      as.numeric()
    ego_df$term_ratio <- ego_df$Count/ego_df$term_gene
  }


  # term order
  ego_df$Description <- factor(ego_df$Description,levels = rev(ego_df$Description))


  # select genes from every term
  sankey_df <- ego_df %>%
    dplyr::select(Description,geneID) %>%
    tidyr::separate_longer_delim(geneID,delim = "/") %>%
    as.data.frame() %>%
    dplyr::mutate(Description = as.character(Description))

  if(keep_all_gene == FALSE){
    sankey_df <- sankey_df %>%
      dplyr::group_by(Description) %>%
      dplyr::slice_head(n = topGenes) %>%
      dplyr::ungroup()
  }

  sankey_df_long <- sankey_df %>% ggsankey::make_long(geneID,Description)

  # order
  sankey_df_long$node <- factor(sankey_df_long$node,
                                levels = c(unique(sankey_df$Description),
                                           unique(sankey_df$geneID)))

  # ============================================================================
  # 2.draw plot
  # ============================================================================
  if(is.null(nodeColor)){
    mycol <- cols4all::c4a('rainbow_wh_rd',length(unique(sankey_df_long$node)))
  }else{
    mycol <- grDevices::colorRampPalette(colors = nodeColor)(length(unique(sankey_df_long$node)))
  }

  # sankey plot
  ps <-
    ggplot(data = sankey_df_long,
           mapping = aes(x = x,
                         next_x = next_x,
                         node = node,
                         next_node = next_node,
                         fill = factor(node),
                         label = node)) +
    geom_sankey(flow.alpha = 0.5,
                flow.fill = flow_fill,
                # flow.color = 'grey',
                # width = 0.1,
                node.fill = mycol) +
    theme_void() +
    scale_x_discrete(expand = expansion(mult = sankeyExpand)) +
    theme(legend.position = 'none')

  # extract data from sankey plot
  ps_data <- ggplot_build(ps)
  ps_data_info <- ps_data$data[[2]]

  # re-adding text labels
  ps2 <-
    ps +
    geom_text(data = ps_data_info,
              mapping = aes(next_x = 0,next_node = 0,
                            x = xmin,y = (ymin + ymax)/2,
                            label = label,hjust = 1),
              fontface = "bold",size = nodeSize)

  # go plot
  if(is.null(goCol)){
    pcol <- scale_fill_viridis_c(option = "plasma",direction = -1)
  }else{
    pcol <- scale_fill_gradient(low = goCol[1],high = goCol[2])
  }

  if(!is.null(geom_layer)){
    plyer <- geom_layer
  }else{
    if(enrich_type == "gsea"){
      plyer <- geom_point(aes(x = -log10(pvalue),y = Description,
                              size = -log10(p.adjust),fill = -log10(p.adjust)),
                          color = "black",shape = 21)
    }else{
      plyer <- geom_point(aes(x = -log10(pvalue),y = Description,
                              size = term_ratio,fill = gene_ratio),
                          color = "black",shape = 21)
    }
  }

  pp <-
    ggplot(ego_df) +
    plyer +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(colour = "black"),
          plot.background = element_blank()) +
    ylab("") + pcol

  # ============================================================================
  # calculate relative position
  # ============================================================================
  xmin <- max(ps_data_info$xmax)
  y_range <- subset(ps_data_info,x == "2")
  ymin <- min(c(y_range$ymin,y_range$ymax))
  ymax <- max(c(y_range$ymin,y_range$ymax))


  # insert plot
  cb <-
    ps2 +
    annotation_custom(grob = ggplotGrob(pp),
                      xmin = xmin - xShift,xmax = 2 + sankeyExpand[2],
                      ymin = ymin - downShift,ymax = ymax + upShift)

  return(cb)
}
