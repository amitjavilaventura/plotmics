# expressionHeatmap2() -------------------------

#' @title expressionHeatmap2
#' @author amitjavilaventura
#'
#' @description
#' Function that takes a list of dataframes with expression data and the columns 'Geneid' and 'log2FoldChange'.
#' and draws a heatmap using ggplot.
#'
#' @seealso `plotDendogram`
#' @seealso `stats::hclust`
#' @seealso `stats::dist`
#'
#' @usage `expressionHeatmap2((expr_list,expr_names = names(expr_list),genes = c("Lef1", "Tcf7l1", "Tcf7l2", "Tcf7"),clust_rows = T, clust_cols = F, show_dend_rows = F, show_dend_cols = F, dist_method = "euclidean", hclust_method = "ward.D", write_label = T, label_size = 3, label_color = "black", label_digits = 2, hm_height = length(genes)*10, hm_width = length(expr_list)*10, hm_colors = c("cornflowerblue", "white", "gold3"), legend_scale = c(-1.5, 1.5), legend_breaks_num = 5, legend_breaks_by = .5, legend_midpoint = 0, legend_height = hm_height, legend_title = NULL, dend_cols_prop = .1, dend_rows_prop = .2, title = "", subtitle = "", caption = NULL, xlab = "", ylab = NULL, axis_text_size = 10, x_axis_angle = 90))`
#'
#'
#' @param expr_list List of dataframes with, at least, the columns 'Geneid' and 'log2FoldChange'. Better if it's a named list.
#' @param expr_names Character of namesof equal length to 'expr_list'. Names of the elements in 'expr_list'. Default: names(expr_list)
#' @param genes Character. Names of the genes to be plotted. They must present in the column 'Geneid' of 'df'.
#' @param clust_rows Logical of length 1. Whether to cluster the rows (TRUE) or not (FALSE) using `hclust`. Default: T.
#' @param clust_cols Logical of length 1. Whether to cluster the columns (TRUE) or not (FALSE) using `hclust`. Default: T.
#' @param show_dend_rows Logical of length 1. Whether to draw the dendogram of the rows clustering (TRUE) or not (FALSE). It works only with 'clust_rows = T'. Default: T.
#' @param show_dend_cols Logical of length 1. Whether to draw the dendogram of the columns clustering (TRUE) or not (FALSE). It works only with 'clust_cols = T'. Default: T.
#' @param hclust_method Character of length 1. Method for hierarchical clustering. One of `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"` (= UPGMA), `"mcquitty"` (= WPGMA), `"median"` (= WPGMC) or `"centroid"` (= UPGMC). Default: "ward.D"
#' @param dist_method Character of length 1. Distance calculation. One of `"euclidean"`, `"maximum"`, `"manhattan"`, "`canberra"`, `"binary"` or `"minkowski"`. Default: "euclidean".
#' @param write_label Logical of length 1. Whether to write the expression values in the heatmap (TRUE) or not (FALSE). Default: T.
#' @param label_size Numerical of length 1. Size of the expression values written in each cell of the heatmap. Default: 3.
#' @param label_color Character of length 1. Color of the expression values written in each cell of the heatmap. Default: "black".
#' @param label_digits Numerical of length 1. Number of digits the expression values are rounded to. Defalut: 2.
#' @param hm_height Numerical of length 1. Height of the heatmap in mm. Default: length(genes)*10.
#' @param hm_width Numerical of length 1. Width of the heatmap in mm. Default: (ncol(df)-1)*10.
#' @param hm_colors Character of length 3. Colors in the lower limit, midpoint (defined by 'legend_midpoint') and higher limit, respectively. Default: c("cornflowerblue", "white", "gold3").
#' @param legend_scale Numerical of length 2 or NULL. If NULL, the color scale of the heatmap will take the minimum and the maximum values as limits. If numerical, the color scale will take the first as the lower limit and the second element as the higher limit. Default: c(-1.5, 1.5).
#' @param legend_breaks_num Numerical of length 1. Only if 'legend_scale = NULL'. The number of breaks you want in the legend. Default: 5.
#' @param legend_breaks_by Numerical of length 1. Only 'legend_scale' is numerical (e.g. `c(-1,1)`). The distance between the breaks of the legend. Default: .5.
#' @param legend_midpoint Numerical of length 1. Only if scale is not NULL. Point where the central color of the legend will placed. Default: 0
#' @param legend_height Numerical of length 1. Height of the legend which, by default, is the height of the heatmap. Default: hm_height.
#' @param legend_title Character of length 1 or NULL. Title of the legend, placed in the right part of it. Default: NULL.
#' @param dend_cols_prop Numerical of length 1. Proportion (from 0 to 1) of the columns dendogram compared to the heatmap height. Default: 0.1.
#' @param dend_rows_prop Numerical of length 1. Proportion (from 0 to 1) of the rows dendogram compared to the heatmap width Default: 0.2.
#' @param title Character of length 1 or NULL. Title of the plot. Default: "".
#' @param subtitle Character of length 1 or NULL. Subtitle of the plot. Default: "".
#' @param caption Character of length 1 or NULL. Caption of the plot; placed at the bottom. Default: NULL.
#' @param xlab Character of length 1 or NULL. Title of the X axis. Default: "
#' @param ylab Character of length 1 or NULL. Title of the Y axis. If row dendogram is plotted, the Y axis is placed to the right. Default: NULL.
#' @param legend_title_size Numerical of length 1. Size of the legend title. Default: 8
#' @param title_size Numerical of length 1. Size of the plot title. Default: 13.
#' @param subtitle_size Numerical of length 1. Size of the plot subtitle. Default: 13.
#' @param caption_size Numerical of length 1. Size of the plot caption. Default: 7.
#' @param axis_title_size Numerical of length 1. Size of the axis titles. Default: 10
#' @param axis_text_size Numerical of length 1. Size of the text in the axes. Default: 10
#' @param title_hjust Numerical of length 1. Justification of the plot title and subtitle. Default: 0.
#' @param scale Character of length 1 or NULL. One of c("rows", "cols"). If "rows", it scales data by row; if "cols", it scales data by columns; if NULL (the default), it does not scale the data. Default: NULL.
#'
#' @export

expressionHeatmap2 <- function(expr_list,
                               expr_names = names(expr_list),
                               genes = c("Lef1", "Tcf7l1", "Tcf7l2", "Tcf7"),
                               clust_rows = F,
                               clust_cols = F,
                               show_dend_rows = F,
                               show_dend_cols = F,
                               dist_method = "euclidean",
                               hclust_method = "ward.D",
                               write_label = T,
                               label_size = 3,
                               label_color = "black",
                               label_digits = 2,
                               hm_height = length(genes)*10,
                               hm_width = length(expr_list)*10,
                               hm_colors = c("cornflowerblue", "white", "gold3"),
                               legend_scale = c(-1.5, 1.5),
                               legend_breaks_num = 5,
                               legend_breaks_by = .5,
                               legend_midpoint = 0,
                               legend_height = hm_height,
                               legend_title = NULL, legend_title_size = 8,
                               dend_cols_prop = .1,
                               dend_rows_prop = .2,
                               title = "", title_size = 13, title_hjust = 0,
                               subtitle = "", subtitle_size = 12,
                               caption = NULL, caption_size = 7,
                               xlab = "",
                               ylab = NULL,
                               axis_title_size = 11, axis_text_size = 10,
                               scale = NULL){

  # Load required packages
  require(plyr)
  require(dplyr)
  require(purrr)
  require(reshape2)
  require(ggplot2)
  require(ggpubr)
  require(ggh4x)
  require(ggdendro)
  require(patchwork)

  # Check that inputs are OK
  if(!is.list(expr_list)){ stop("'expr_list' must be a (named) list of dataframes with, at least, the columns 'Geneid' and 'log2FoldChange'.") }
  else if(!(sum(sapply(expr_list, class) == "data.frame") == length(expr_list))){ stop("'expr_list' must be a (named) list of dataframes with, at least, the columns 'Geneid' and 'log2FoldChange'.") }
  else if(!(sum(sapply(expr_list, colnames) == "Geneid") == length(expr_list))){ stop("Each data frame in the list must have the column 'Geneid' with the names of the genes.") }
  else if(!(sum(sapply(expr_list, colnames) == "log2FoldChange") == length(expr_list))){ stop("Each data frame in the list must have the column 'log2FoldChange' with the expression values of each gene.") }
  if(show_dend_rows & !clust_rows){ stop("To plot the rows dendogram (`show_dend_rows = T`), `clust_rows` must be TRUE") }
  if(show_dend_cols & !clust_cols){ stop("To plot the columns dendogram (`show_dend_cols = T`), `clust_cols` must be TRUE") }


  # Bind Log2FCs from all genes in all contrasts
  expr <- expr_list %>%
    purrr::set_names(expr_names) %>%
    purrr::map(~dplyr::select(.x, Geneid, log2FoldChange)) %>%
    purrr::map(~dplyr::filter(.x, Geneid %in% genes)) %>%
    purrr::imap(~magrittr::set_colnames(.x, value = c("Geneid", .y))) %>%
    plyr::join_all(by = "Geneid") %>%
    dplyr::relocate(Geneid)

  ## Scale
  if(!is.null(scale)){
    if("rows" %in% scale) { expr[,2:length(expr)] <- expr[,2:length(expr)] %>% t() %>% scale() %>% t() }
    if("cols" %in% scale) { expr[,2:length(expr)] <- expr[,2:length(expr)] %>% scale() }
  }

  # Melt df with all Log2FCs
  expr.m <- reshape2::melt(expr, variable.name = "Condition", id.vars = "Geneid", value.name = "Expr")

  # Initialize heatmap -----
  # Get the legend breaks.
  if(is.null(legend_scale)){
    lower    <- min(expr.m$Expr)
    higher   <- max(expr.m$Expr)
    midpoint <- (higher+(lower))/2
    breaks   <- round(seq(from = lower, to = higher, by = (higher-lower)/legend_breaks_num))
  } else {
    lower    <- legend_scale[1]
    higher   <- legend_scale[2]
    midpoint <- legend_midpoint
    breaks   <- (seq(from = lower, to = higher, by = legend_breaks_by))
  }

  ## Plot heatmap
  hm <- ggplot(expr.m, aes(Condition, Geneid) ) +
    geom_tile(aes(fill = Expr)) +
    force_panelsizes(rows = unit(hm_height, "mm"),
                     cols = unit(hm_width, "mm"))  +
    coord_fixed() +
    scale_fill_gradient2(limits = c(legend_scale[1], legend_scale[2]), midpoint = midpoint,
                         breaks = breaks,
                         oob = scales::squish,
                         low = hm_colors[1], mid = hm_colors[2], high = hm_colors[3],
                         guide = guide_colorbar(title = legend_title,
                                                title.position = "right",
                                                title.theme = element_text(angle = 270,
                                                                           hjust = .5, vjust = .5,
                                                                           size = legend_title_size),
                                                frame.colour = "black", frame.linewidth = 1.5,
                                                ticks.colour = NA,
                                                barheight = unit(legend_height, "mm"))) +

    theme_pubr(border = T, legend = "right", margin = T) + # format the plot
    theme(plot.title = element_text(face = "bold", size = title_size, hjust = title_hjust),
          plot.subtitle = element_text(face = "italic", size = subtitle_size, hjust = title_hjust),
          axis.title = element_text(face = "bold", size = axis_title_size),
          axis.text  = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = 90, hjust = .5, vjust = .5),
          panel.border = element_rect(size = 1.1)) +
    xlab(xlab) + ylab(ylab)

  ##Write lfc values
  if(write_label){ hm <- hm + geom_text(aes(label = round(Expr, label_digits)), size = label_size, color = label_color) }

  # Hierarchical clustering -----
  if(clust_rows){
    # Hierarchical clustering of rows alone
    hclust_rows <- hclust( dist(expr[-1], method = dist_method), method = hclust_method )

    # Position of Y axis
    if(show_dend_rows) { ytext_pos = "right" }
    else { ytext_pos = "left" }

    # Order rows
    hm <- hm +
      scale_y_discrete(limits=expr.m$Geneid[hclust_rows$order],
                       labels = expr.m$Geneid, position = ytext_pos,
                       expand = c(0,0))

  } else { hm <- hm + scale_y_discrete(expand = c(0,0)) }

  if(clust_cols){
    # Hierarchical clustering of columns alone
    hclust_cols <- hclust( dist(t(expr[-1]), method = dist_method), method = hclust_method )

    # Order columns
    hm <- hm +
      scale_x_discrete(limits = colnames(expr[,-1])[hclust_cols$order],
                       labels = colnames(expr[,-1])[hclust_cols$order],
                       expand = c(0,0))
  } else { hm <- hm + scale_x_discrete(expand = c(0,0)) }

  ### Plot dendograms
  if(clust_rows & clust_cols & show_dend_rows & show_dend_cols){
    dend_rows <- plotDendogram(hclust_rows, axis = "rows")
    dend_cols <- plotDendogram(hclust_cols, axis = "cols")

    # Plot heatmap and dendogram(s) together with patchwork and the created layout
    heatmap <- hm +
      inset_element(dend_cols, left = 0, bottom = 1, right = 1, top = 1+dend_cols_prop) +
      inset_element(dend_rows, left = -dend_rows_prop, bottom = 0, right = 0, top = 1) +
      plot_annotation(title = title, subtitle = subtitle, caption = caption) &
      theme(plot.margin = unit(x = c(1,1,1,1), units = "mm"))

  } else if(clust_cols & show_dend_cols) {
    dend_cols <- plotDendogram(hclust_cols, axis = "cols")

    # Plot heatmap and dendogram(s) together with patchwork and the created layout
    heatmap <- hm +
      inset_element(dend_cols, left = -0.1, bottom = 1, right = 1.1, top = 1+dend_cols_prop) +
      plot_annotation(title = title, subtitle = subtitle, caption = caption) &
      theme(plot.margin = unit(x = c(1,1,1,1), units = "mm"))

  } else if(clust_rows & show_dend_rows) {
    dend_rows <- plotDendogram(hclust_rows, axis = "rows")

    # Plot heatmap and dendogram(s) together with patchwork and the created layout
    heatmap <- hm +
      inset_element(dend_rows, left = -dend_rows_prop, bottom = 0, right = 0, top = 1) +
      plot_annotation(title = title, subtitle = subtitle, caption = caption) &
      theme(plot.margin = unit(x = c(1,1,1,1), units = "mm"))

  } else {
    heatmap <- hm +
      labs(title = title, subtitle = subtitle, caption = caption) &
      theme(plot.margin = unit(x = c(1,1,1,1), units = "mm"))
  }



  return(heatmap)

}
