# expressionHeatmap() -------------------------

#' @title expressionHeatmap
#' @author amitjavilaventura
#'
#' @description
#' Function that takes a data frame with a Geneid column and several columns for expression data (e.g. TPMs, Log2FC...)
#' and draws a heatmap using ggplot.
#'
#' @seealso `plotDendogram`
#' @seealso `stats::hclust`
#' @seealso `stats::dist`
#'
#' @usage expressionHeatmap(expr_df, genes, clust_rows = T, clust_cols = F, show_dend_rows = F, show_dend_cols = F, dist_method = "euclidean", hclust_method = "ward.D", write_label = T, label_size = 4, label_color = "black", label_digits = 2, hm_height = length(genes)*10, hm_width = (ncol(df)-1)*10, hm_colors = c("cornflowerblue", "white", "gold3"), legend_scale = NULL, legend_breaks_num = 5, legend_midpoint = 0, legend_height = hm_height, legend_title = NULL, dend_cols_prop = .1, dend_rows_prop = .2, title = "", subtitle = "", caption = NULL, xlab = "", ylab = NULL, axis_text_size = 10, x_axis_angle = 90)
#'
#' @param expr_df Dataframe with a 'Geneid' column and several columns with numerical expression data for different samples, such as TPMs or Log2FC.
#' @param genes Character. Names of the genes to be plotted. They must present in the column 'Geneid' of 'df'.
#' @param clust_rows Logical of length 1. Whether to cluster the rows (TRUE) or not (FALSE) using `hclust`. Default: T.
#' @param clust_cols Logical of length 1. Whether to cluster the columns (TRUE) or not (FALSE) using `hclust`. Default: F.
#' @param show_dend_rows Logical of length 1. Whether to draw the dendogram of the rows clustering (TRUE) or not (FALSE). It only works with 'clust_rows = T'. Default: F.
#' @param show_dend_cols Logical of length 1. Whether to draw the dendogram of the columns clustering (TRUE) or not (FALSE). It only works with 'clust_cols = T'. Default: F.
#' @param hclust_method Character of length 1. Method for hierarchical clustering. One of `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"` (= UPGMA), `"mcquitty"` (= WPGMA), `"median"` (= WPGMC) or `"centroid"` (= UPGMC). Default: "ward.D"
#' @param dist_method Character of length 1. Distance calculation. One of `"euclidean"`, `"maximum"`, `"manhattan"`, "`canberra"`, `"binary"` or `"minkowski"`. Default: "euclidean".
#' @param write_label Logical of length 1. Whether to write the expression values in the heatmap (TRUE) or not (FALSE). Default: T.
#' @param label_size Numerical of length 1. Size of the expression values written in each cell of the heatmap. Default: 3.
#' @param label_color Character of length 1. Color of the expression values written in each cell of the heatmap. Default: "black".
#' @param label_digits Numerical of length 1. Number of digits the expression values are rounded to. Defalut: 2.
#' @param hm_height Numerical of length 1. Height of the heatmap in mm. Default: length(genes)*10.
#' @param hm_width Numerical of length 1. Width of the heatmap in mm. Default: (ncol(df)-1)*10.
#' @param hm_colors Character of length 3. Colors in the lower limit, midpoint (defined by 'legend_midpoint') and higher limit, respectively. Default: c("cornflowerblue", "white", "gold3").
#' @param na_color Charachter of length 1. Color for the NA values. Default: "gray".
#' @param border_color Character of length 1 or NA. Color of the border of each cell in the heatmap. Default: NA.
#' @param legend_scale Numerical of length 2 or NULL. If NULL, the color scale of the heatmap will take the minimum and the maximum values as limits. If numerical, the color scale will take the first as the lower limit and the second element as the higher limit. Default: NULL.
#' @param legend_breaks_num Numerical of length 1. Only 'legend_scale' is NULL. The number of breaks you want in the legend. Default: 5.
#' @param legend_breaks_by Numerical of length 1. Only 'legend_scale' is numerical (e.g. `c(-1,1)`). The distance between the breaks of the legend. Default: .5.
#' @param legend_midpoint Numerical of length 1. Only if 'legend_scale' is not NULL. Point where the central color of the legend will placed. Default: 0
#' @param legend_height Numerical of length 1. Height of the legend which, by default, is the height of the heatmap. Default: hm_height.
#' @param legend_title Character of length 1 or NULL. Title of the legend, placed in the right part of it. Default: NULL.
#' @param dend_cols_prop Numerical of length 1. Proportion (from 0 to 1) of the columns dendogram compared to the heatmap height. Default: 0.1.
#' @param dend_rows_prop Numerical of length 1. Proportion (from 0 to 1) of the rows dendogram compared to the heatmap width Default: 0.2.
#' @param title Character of length 1 or NULL. Title of the plot. Default: "".
#' @param subtitle Character of length 1 or NULL. Subtitle of the plot. Default: "".
#' @param caption Character of length 1 or NULL. Caption of the plot; placed at the bottom. Default: NULL.
#' @param xlab Character of length 1 or NULL. Title of the X axis. Default: "
#' @param ylab Character of length 1 or NULL. Title of the Y axis. If row dendogram is plotted, the Y axis is placed to the right. Default: NULL.
#' @param title_size Numerical of length 1. Size of the plot title. Default: 13.
#' @param subtitle_size Numerical of length 1. Size of the plot subtitle. Default: 12.
#' @param caption_size Numerical of length 1. Size of the plot caption. Default: 7.
#' @param axis_title_size Numerical of length 1. Size of the axes titles. Default: 11.
#' @param axis_text_size Numerical of length 1. Size of the text in the axes. Default: 10.
#' @param title_hjust Numerical of length 1. Horizontal justification for the title and the subtitle. Default: 0.
#' @param scale Character of length 1 or NULL. One of c("rows", "cols"). If "rows", it scales data by row; if "cols", it scales data by columns; if NULL (the default), it does not scale the data. Default: NULL.
#'
#' @export

expressionHeatmap <- function(expr_df,
                              genes,
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
                              hm_width = (ncol(expr_df)-1)*10,
                              hm_colors = c("cornflowerblue", "white", "gold3"),
                              na_color = "gray",
                              border_color = NA,
                              legend_scale = NULL,
                              legend_breaks_num = 5,
                              legend_breaks_by = .5,
                              legend_midpoint = 0,
                              legend_height = hm_height,
                              legend_title = NULL,
                              legend_title_size = 9,
                              dend_cols_prop = .1,
                              dend_rows_prop = .2,
                              title = "",
                              subtitle = "",
                              caption = NULL,
                              xlab = "",
                              ylab = NULL,
                              title_size = 13,
                              subtitle_size = 12,
                              caption_size = 7,
                              axis_text_size = 10,
                              axis_title_size = 11,
                              title_hjust = 0,
                              scale = NULL,
                              show_gene_names = T){

  # Load requireed packages
  require(dplyr)
  require(reshape2)
  require(ggplot2)
  require(ggh4x)
  require(ggpubr)
  require(patchwork)

  # Check that inputs are OK -----
  if(!is.data.frame(expr_df)){ stop("'df' must be a data frame with 'Geneid' in the first column and the other columns with the expression values to plot.") }
  if(!("Geneid" %in% colnames(expr_df))){ stop("'df' must have a character first column named 'Geneid'") }
  if(show_dend_rows & !clust_rows){ stop("To plot the rows dendogram (`show_dend_rows = T`), `clust_rows` must be TRUE") }
  if(show_dend_cols & !clust_cols){ stop("To plot the columns dendogram (`show_dend_cols = T`), `clust_cols` must be TRUE") }

  # Format data frame -----
  ## Filter dataframe to get only desired genes
  expr <- expr_df %>%
    dplyr::relocate(Geneid) %>%
    dplyr::filter(Geneid %in% genes)



  ## Scale
  if(!is.null(scale)){
    if("rows" %in% scale) { expr[,2:length(expr)] <- expr[,2:length(expr)] %>% t() %>% scale() %>% t() }
    if("cols" %in% scale) { expr[,2:length(expr)] <- expr[,2:length(expr)] %>% scale() }
  }

  ## Melt data frame
  expr.m <- expr %>%
    reshape2::melt(value.name = "Expr", variable.name = "Condition")

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
    breaks   <- round(seq(from = lower, to = higher, by = legend_breaks_by))
  }

  # Initialize the heatmap with ggplot
  hm <- ggplot(expr.m, aes(Condition, Geneid) ) +
    geom_tile(aes(fill = Expr), color = border_color) +
    ggh4x::force_panelsizes(rows = unit(hm_height, "mm"),  # force height of the heatmap
                            cols = unit(hm_width, "mm")) + # force width of the heatmap
    coord_fixed() +  # fix the coordinates
    scale_fill_gradient2(limits = c(lower, higher), midpoint = midpoint, # setup the legend
                         breaks = breaks, na.value = na_color,
                         oob = scales::squish, ## to put out of bound values into scale
                         low = hm_colors[1], high = hm_colors[3], mid = hm_colors[2],
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
          plot.caption = element_text(size = caption_size),
          axis.title = element_text(face = "bold", size = axis_title_size),
          axis.ticks.y = element_blank(),
          axis.text  = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = 90, hjust = .5, vjust = .5),
          panel.border = element_rect(size = 1.1)) +
    xlab(xlab) + ylab(ylab)

  # Write thee label of the expression values
  if(write_label){ hm <- hm + geom_text(aes(label = round(Expr, label_digits)), size = label_size, color = label_color) }

  # Remove row names
  if(!show_gene_names){ hm <- hm + theme(axis.text.y = element_blank()) }

  # Hierarchical clustering -----
  if(clust_rows){
    # Hierarchical clustering of rows alone
    hclust_rows <- hclust( dist(expr[-1], method = dist_method), method = hclust_method )

    # Position of the Y axis
    if(show_dend_rows) { ytext_pos = "right" }
    else { ytext_pos = "left" }

    # Order rows
    hm <- hm +
      scale_y_discrete(limits=expr.m$Geneid[hclust_rows$order],
                       labels = expr.m$Geneid[hclust_rows$order], position = ytext_pos,
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


  # Return heatmap
  return(heatmap)

}
