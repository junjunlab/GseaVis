#' @title volcanoGsea
#' @name volcanoGsea
#' @author Jun Zhang
#' @param data GSEA enrich object from clusterProfiler, defalut is NULL.
#' @param NES.cutoff NES cutoff to select significant terms, defalut is 1.
#' @param pvalue.cutoff pvalue cutoff to select significant terms, defalut is NULL.
#' @param p.adjust.CUTOFF adjusted pvalue cutoff to select significant terms, defalut is 0.05.
#' @param nudge.y y shift to ajust label, defalut is c(0,0).
#' @param topN top term to show, defalut is 5.
#' @param point.size point size, defalut is 3.
#' @param point.color point color, defalut is c('#CC3333','#CCCCCC','#0099CC').
#' @param ... other arguments passed by geom_text_repel.
#'
#' @return a ggplot object.
#' @export
volcanoGsea <- function(data = NULL,
                        NES.cutoff = 1,
                        pvalue.cutoff = NULL,
                        p.adjust.CUTOFF = 0.05,
                        nudge.y = c(0,0),
                        topN = 5,
                        point.size = 3,
                        point.color = c('#CC3333','#CCCCCC','#0099CC'),
                        ...){
  # to data.frame
  df <- data.frame(data)

  # assign type
  if(is.null(pvalue.cutoff)){
    pcol <- "p.adjust"
    df$type <- dplyr::case_when(df$NES >= NES.cutoff & df$p.adjust < p.adjust.CUTOFF ~ 'sig-Activated',
                                df$NES <= -NES.cutoff & df$p.adjust < p.adjust.CUTOFF ~ 'sig-Repressed',
                                TRUE ~ 'none sig')

    # add -log10P
    df <- df %>% dplyr::mutate(logp = -log10(p.adjust))
  }else{
    pcol <- "pvalue"
    df$type <- dplyr::case_when(df$NES >= NES.cutoff & df$pvalue < pvalue.cutoff ~ 'sig-Activated',
                                df$NES <= -NES.cutoff & df$pvalue < pvalue.cutoff ~ 'sig-Repressed',
                                TRUE ~ 'none sig')

    # add -log10P
    df <- df %>% dplyr::mutate(logp = -log10(pvalue))
  }


  # top term
  topterm <- purrr::map_df(c('sig-Activated','sig-Repressed'),function(x){
    tmp <- df %>%
      dplyr::filter(type == x) %>%
      dplyr::arrange(.data[[pcol]]) %>%
      dplyr::slice_head(n = topN)

    # add nudge_y
    if(nrow(tmp) == 0){
      return(NULL)
    }else{
      if(x == 'sig-Activated'){
        tmp$nudgey <- nudge.y[1]
      }else{
        tmp$nudgey <- nudge.y[2]
      }
      return(tmp)
    }
  })

  # plot
  p <-
    ggplot2::ggplot(df,ggplot2::aes_string(x = "logp",y = "NES")) +
    ggplot2::geom_point(ggplot2::aes(color = type),alpha = 0.5,size = point.size) +
    ggplot2::geom_vline(xintercept = -log10(ifelse(p.adjust.CUTOFF,p.adjust.CUTOFF,pvalue.cutoff)),
                        size = 1,lty = 'solid',color = 'grey75') +
    ggplot2::geom_hline(yintercept = c(-NES.cutoff,NES.cutoff),size = 1,lty = 'dashed',color = 'grey75') +
    ggplot2::scale_colour_manual(name = '',
                                 values = c('sig-Activated' = point.color[1],
                                            'none sig' = point.color[2],
                                            'sig-Repressed' = point.color[3])) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black", size = ggplot2::rel(1)),
                   axis.title = ggplot2::element_text(colour = "black", size = ggplot2::rel(1.1)),
                   legend.position = 'top') +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggplot2::ylab('Normalized enriched score') +
    ggplot2::xlab(paste("-log10",pcol,sep = ' ')) +
    ggrepel::geom_text_repel(data = topterm,
                             ggplot2::aes_string(x = "logp",y = "NES",label = "Description"),
                             fontface = 'italic',
                             max.overlaps = 80,
                             force = 50,
                             nudge_y = topterm$nudgey,
                             min.segment.length = ggplot2::unit(0.1, "cm"),
                             ...)

  return(p)
}
