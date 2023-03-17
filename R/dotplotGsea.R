#' @title dotplotGsea
#' @name dotplotGsea
#' @author Jun Zhang
#' @param data GSEA enrich object from clusterProfiler, defalut is NULL.
#' @param pval pvalue cutoff to select significant terms, defalut is NULL.
#' @param pajust adjusted pvalue cutoff to select significant terms, defalut is 0.05.
#' @param order.by the X axis, defalut is "GeneRatio".
#' @param str.width the width of term name, defalut is 50.
#' @param base_size theme base size, defalut is 12.
#' @param topn show the top terms, defalut is NULL.
#' @param scales facet scales, defalut is "free_x".
#' @param add.seg whether add segment line to point, defalut is "FALSE".
#' @param line.col segment line color, defalut is "grey80".
#' @param line.size segment line size, defalut is 1.5.
#' @param line.type segment line type, defalut is "solid".
#'
#' @return a ggplot object.
#' @export

globalVariables(c('.data', 'Description', 'GeneRatio', 'NES', 'p.adjust', 'pvalue', 'type'))
dotplotGsea <- function(data = NULL,
                        pval = NULL,
                        pajust = 0.05,
                        order.by = "GeneRatio",
                        str.width = 50,
                        base_size = 12,
                        topn = NULL,
                        scales = "free_x",
                        add.seg = FALSE,
                        line.col = 'grey80',
                        line.size = 1.5,
                        line.type = 'solid'){
  # select sig terms
  if(is.null(pval)){
    df <- data.frame(data) %>%
      dplyr::filter(p.adjust <= pajust)

    geomPoint <- ggplot2::geom_point(ggplot2::aes(color = p.adjust,size = Count))
  }else{
    df <- data.frame(data) %>%
      dplyr::filter(pvalue <= pval)

    geomPoint <- ggplot2::geom_point(ggplot2::aes(color = pvalue,size = Count))
  }

  # add type
  df$type <- ifelse(df$NES > 0, "activated","suppressed")

  # add gene ratio
  Count <- strsplit(df$core_enrichment,split = "\\/") %>%
    lapply(., length) %>% unlist()

  df$Count <- Count
  df$GeneRatio <- Count/df$setSize

  # filter top n
  if(order.by == "NES"){
    if(!is.null(topn)){
      df.up <- df %>%
        dplyr::filter(type == "activated") %>%
        dplyr::arrange(dplyr::desc(.data[[order.by]])) %>%
        dplyr::slice_head(n = topn)

      df.down <- df %>%
        dplyr::filter(type == "suppressed") %>%
        dplyr::arrange(.data[[order.by]]) %>%
        dplyr::slice_head(n = topn)

      df <- rbind(df.up,df.down)
    }
  }else{
    if(!is.null(topn)){
      df <- df %>%
        dplyr::group_by(type) %>%
        dplyr::arrange(dplyr::desc(.data[[order.by]])) %>%
        dplyr::slice_head(n = topn)
    }
  }

  # wrap term name
  df <- df %>% dplyr::arrange(.data[[order.by]])
  df$Description <- stringr::str_wrap(df$Description,width = str.width)
  df$Description <- factor(df$Description,levels = df$Description)
  df$type <- factor(df$type,levels = c('suppressed','activated'))

  # plot
  p <-
    ggplot2::ggplot(df,ggplot2::aes_string(x = order.by,
                                           y = "Description"))

  # whether add aegment
  if(add.seg == TRUE){
    if(order.by == "GeneRatio"){
      aes.seg <- ggplot2::aes(x = 0,xend = GeneRatio,
                              y = Description,yend = Description)
    }else{
      aes.seg <- ggplot2::aes(x = 0,xend = NES,
                              y = Description,yend = Description)
    }

    p1 <- p +
      ggplot2::geom_segment(aes.seg,
                            color = line.col,
                            size = line.size,
                            lty = line.type)
  }else{
    p1 <- p
  }

  # add main
  p2 <- p1 +
    geomPoint +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(color = 'black',size = base_size),
                   strip.background = ggplot2::element_rect(fill = 'grey90'),
                   strip.text = ggplot2::element_text(size = base_size + 1)) +
    ggplot2::facet_wrap(~type,nrow = 1,scales = scales) +
    ggplot2::scale_size(range = c(3,8)) +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(color = 'grey60')),
                    color = ggplot2::guide_colorbar(frame.colour = 'black')) +
    ggsci::scale_color_gsea(reverse = T) +
    ggplot2::ylab('')

  # if(order.by == "NES"){
  #   vl.data <- data.frame(type = c('suppressed','activated'),
  #                         x = c(-1,1))
  #   p3 <- p2 +
  #     # facetted_pos_scales(x = list(scale_x_continuous(),scale_x_reverse()))
  #     # geom_vline(data = vl.data,aes(xintercept = x),lty = 'dashed',color = 'grey50',size = 1)
  # }else{
  #   p3 <- p2
  # }

  # return
  res <- list(df = df,plot = p2)
  return(res)
}
