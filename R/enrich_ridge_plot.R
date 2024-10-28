globalVariables(c("core_enrichment"))

#' Generate Enrichment Ridge Plot
#'
#' This function creates a ridge plot for visualizing gene enrichment analysis results using either `enrichResult` or `gseaResult` objects.
#' Ridge plots allow for the visualization of the distribution of genes across multiple terms and enriched pathways.
#'
#' @param object An `enrichResult` or `gseaResult` object obtained from enrichment analysis.
#' @param terms_ID A vector of term IDs (corresponding to enriched pathways or terms) to be included in the plot.
#' @param gene_list A data frame with two columns("id","logfc") providing the gene fold changes (log fold changes).
#' @param term_name_width (optional) Maximum width of the term name for wrapping text in the y-axis. Default is 60.
#' @param geom_density_ridges_params (optional) A list of additional parameters passed to `geom_density_ridges` to customize the ridge plot appearance.
#'
#' @details This function accepts two types of enrichment result objects: `enrichResult` and `gseaResult`. For both types, the function handles gene identifiers
#' (either raw gene IDs or human-readable names) and generates the appropriate data frames with corresponding log fold change values.
#' It then filters and plots the specified terms.
#'
#' The ridge plot helps visualize how genes for each term distribute across log fold change values, and it also highlights the statistical
#' significance of the terms using a color gradient based on the p-value.
#'
#' @return A `ggplot` object which contains the ridge plot visualizing the log fold change distribution of enriched terms.
#'
#' @importFrom tidyr separate_longer_delim
#' @importFrom ggridges position_points_jitter
#' @importFrom utils modifyList
#'
#' @export
enrich_ridge_plot <- function(object = NULL,
                              terms_ID = NULL,
                              gene_list = NULL,
                              term_name_width = 60,
                              geom_density_ridges_params = list()){
  # ============================================================================
  # prepare data
  # ============================================================================
  # check type of enrichment object
  if(inherits(object,"enrichResult")){
    # check gene_list
    if(is.null(gene_list)) stop("Please supply gene list with fold changes!")

    if(object@readable == TRUE){
      idm <- data.frame(id = names(object@gene2Symbol),gname = object@gene2Symbol)
      gids <- idm %>%
        dplyr::left_join(y = gene_list,by = "id")

      ids <- "gname"
    }else{
      idm <- data.frame(id = object@gene)
      gids <- idm %>%
        dplyr::left_join(y = gene_list,by = "id")

      ids <- "id"
    }

    # filter data
    df <- data.frame(object) %>%
      tidyr::separate_longer_delim(cols = geneID,delim = "/") %>%
      dplyr::left_join(y = gids,by = c("geneID" = ids))

    df_plot <- subset(df,ID %in% terms_ID)
  }else if(inherits(object,"gseaResult")){
    if(object@readable == TRUE){
      idm <- data.frame(id = names(object@gene2Symbol),gname = object@gene2Symbol)
      gids <- data.frame(id = names(object@geneList),logfc = object@geneList) %>%
        dplyr::left_join(y = idm,by = "id")

      ids <- "gname"
    }else{
      gids <- data.frame(id = names(object@geneList),logfc = object@geneList)

      ids <- "id"
    }

    # filter data
    df <- data.frame(object) %>%
      tidyr::separate_longer_delim(cols = core_enrichment,delim = "/") %>%
      dplyr::left_join(y = gids,by = c("core_enrichment" = ids))

    df_plot <- subset(df,ID %in% terms_ID)
  }

  # ============================================================================
  # plot
  # ============================================================================

  ggplot(df_plot, aes(x = logfc, y = Description)) +
    do.call(ggridges::geom_density_ridges,
            modifyList(list(mapping = aes(fill = -log10(pvalue)),
                            jittered_points = TRUE,
                            position = position_points_jitter(width = 0.05, height = 0),
                            scale = 1,
                            point_shape = '|', point_color = "#CC0033"),
                       geom_density_ridges_params)) +
    # geom_density_ridges(
    #   aes(fill = -log10(pvalue)),
    #   jittered_points = TRUE,
    #   position = position_points_jitter(width = 0.05, height = 0),
    #   scale = 1,
    #   point_shape = '|', point_color = "grey50"
    # ) +
    scale_y_discrete(labels = function(x){stringr::str_wrap(x,width = term_name_width)}) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank()) +
    xlab("log2FoldChange") + ylab("") +
    scale_fill_viridis_c(option = "magma")
}
