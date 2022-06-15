# expressionCor() -------------------------

#' @title expressionCor
#' @author amitjavilaventura
#'
#' @description
#' Function that takes a data frame with a Geneid column and several columns for
#' expression data of different samples (e.g. TPMs, Log2FC...) and draws a correlation plot
#'
#'
#' @param df Dataframe with a 'Geneid' column and several columns with numerical expression data for different samples, such as TPMs or Log2FC.
#' @param genes Character or NULL. Names of the genes to be used for the correlation analysis. They must present in the column 'Geneid' of 'df'. If NULL, all genes in the data frame are used. Default: NULL
#' @param samples Character or NULL. Names of the samples to be used for the correlation analysis. They must be column names of 'df'. If NULL, all samples are used. Default: NULL.
#' @param corr_method Character of length 1. Correlation method to be passed through `cor()`. One of `"pearson"`, `"spearman"` or `"kendall"`. Default: "pearson".
#' @param plot_type Character of length 1. Type of the correlation plot. One of `"full"`, `"upper"` or `"lower"`. Default: "lower".
#' @param plot_diagonal Logical of length 1. Whether to plot the diagonal or not if the plot type is upper or lower. Default: TRUE.
#' @param plot_size Numerical of length 1. Size in milimeters of the plot. The same number is used for plot height and plot width Default: 85.
#' @param plot_border Logical of length 1. Whether to draw a border in the panel or not. Default: TRUE.
#' @param plot_colors Character of length 3. Colors of the lower, midpoint and higher limits of the scale. Default: c("Gold3", "White", "Cornflowerblue").
#' @param plot_title Character of length 1 or NULL. Title of the plot. Defalut: NULL.
#' @param plot_subtitle Character of length 1 or NULL. Subtitle of the plot. Defalut: NULL.
#' @param plot_caption Character of length 1 or NULL. Caption of the plot. Defalut: NULL.
#' @param cell_border Charachter of length 1, NULL or NA. Color of the cell border. Default: "Black".
#' @param legend_pos Charachter of length 1. Position of the legend, to be passed through `ggpubr::theme_pubr()`. Default: "right".
#' @param legend_size Numerical of length 2. Size of the legend. The first element is the width and the second is the height of the legend bar. Default: c(8,plot_size)
#' @param legend_limits Numerical of length 2. Lower and upper limits of the correlation scale. Default: c(-1,1).
#' @param legend_breaks_by Numerical of length 1. Size of the breaks in the legend. Default: 0.5.
#' @param legend_title Character of length 1 or NULL. Title of the legend. Default: `paste(stringr::str_to_sentence(corr_method), "correlation", sep = " ")`.
#' @param legend_title_size Numerical of length 1. Size of the legend title. Default: 10.
#' @param coeffs_color Character of length 1 or NULL. Color of the correlation coeficients to be plotted. If NULL, no correlation coefficients are plotted. Default: "Black".
#' @param coeffs_size Numerical of length 1. Size of the correlation coefficients. Default: 4.
#' @param title_hjust Numerical of length 1. Horizontal justification of the title and the subtitle. Default: 0.5.
#' @param title_face Character of length 1. Face of the title text. One of "plain", "italic", "bold". Default: "plain".
#' @param title_size Numerical of length 1. Size of the title text. Default: 12.
#' @param subtitle_face Character of length 1. Face of the subtitle text. One of "plain", "italic", "bold". Default: "italic".
#' @param subtitle_size Numerical of length 1. Size of the subtitle text. Default: 11.
#' @param caption_size Numerical of length 1. Size of the plot caption. Default: 6.
#' @param axis_text_size Numerical of length 1. Size of the text in the axes. Default: 7.
#' @param axis_text_color Character of length 1. Color of the text in the axes. Default: "Black".
#' @param by_groups Logical of length 1. If TRUE, the function searches an extra column named 'group', whose values will be used to separate the plot into facets using facet_wrap(). Default: FALSE.

#'
#' @export
expressionCor <- function(df,
                          genes              = NULL,
                          samples            = NULL,
                          samples_order      = samples,
                          corr_method        = "pearson",
                          plot_type          = "lower",
                          plot_diagonal      = TRUE,
                          plot_size          = 85,
                          plot_border        = TRUE,
                          plot_colors        = c("Gold3", "White", "Cornflowerblue"),
                          plot_title         = NULL,
                          plot_subtitle      = NULL,
                          plot_caption       = NULL,
                          cell_border        = "Black",
                          legend_pos         = "right",
                          legend_size        = c(8,plot_size),
                          legend_limits      = c(-1,1),
                          legend_breaks_by   = 0.5,
                          legend_title       = paste(stringr::str_to_sentence(corr_method), "correlation", sep = " "),
                          legend_title_size  = 10,
                          coeffs_color       = "Black",
                          coeffs_size        = 4,
                          title_hjust        = .5,
                          title_face         = "plain",
                          title_size         = 12,
                          subtitle_face      = "italic",
                          subtitle_size      = 11,
                          caption_size       = 6,
                          axis_text_size     = 8,
                          axis_text_color    = "black",
                          by_groups          = F) {


  # Load packages -----
  require(dplyr)
  require(tibble)
  require(reshape2)
  require(ggplot2)
  require(ggh4x)
  require(ggpubr)
  require(scales)

  # Check that inputs are OK
  if(by_groups & !"group" %in% colnames(df)) { stop("If 'by_groups' is TRUE, a column named 'group' must be in the input 'df'.") }

  # Format and filter data frame -----
  df_filt <- df %>% na.omit()

  # Filter genes and samples if those are indicated
  if(!is.null(genes)) { df_filt <- df_filt %>% dplyr::filter(Geneid %in% genes) }
  if(!is.null(samples)) { df_filt <- df_filt %>% dplyr::select(c("Geneid", samples)) }

  # If argument genes is NULL, don't use Geneid column
  if(is.null(genes)){
    if("Geneid" %in% colnames(df_filt)) { df_filt <- df_filt %>% dplyr::select(-Geneid)}
    else { df_filt <- df_filt }
  }

  # If the by_groups option is set to FALSE (default)
  if(!by_groups){

    # Remove group column is present
    if("group" %in% colnames(df_filt)){ df_filt <- df_filt %>% dplyr::select(-group) }

    # Do correlation and format the correlation matrix -----
    corr <- df_filt %>% tibble::column_to_rownames("Geneid") %>% cor(method = corr_method)

    # Set plot type as full, upper or lower, including whether or not to plot the diagonal
    if(plot_type == "full"){ corr <- corr }
    else if(plot_type == "upper") { corr[lower.tri(corr, diag = !plot_diagonal)] <- NA }
    else if(plot_type == "lower") { corr[upper.tri(corr, diag = !plot_diagonal)] <- NA }

    # Melt correlation matrix
    corr.m <- reshape2::melt(corr) %>% dplyr::as_tibble()
    # Add order for samples (or columns used to compute the correlation)
    if(!is.null(samples_order)) { corr.m <- corr.m %>% dplyr::mutate(Var1 = factor(Var1, levels = samples_order), Var2 = factor(Var2, levels = samples_order)) }

    # Draw the plot ----
    # Initialize the plot and the squares
    g <- ggplot(corr.m, aes(Var1, Var2, fill = value)) + geom_tile(color = cell_border, na.rm = T)

  } else {

    # Initialize list for correlation values
    corr_list <- list()
    # For each group, do:
    for(i in unique(df_filt$group)){
      # Filter the data frame by group and remove group variable
      df_filt_group <- df_filt %>% dplyr::filter(group == i) %>% dplyr::select(-group)

      # Compute correlation
      if(!is.null(genes)) { corr <- df_filt_group %>% tibble::column_to_rownames("Geneid") %>% cor(method = corr_method) }
      else { corr <- df_filt_group %>% cor(method = corr_method) }

      # Format the correlation matrix to draw a full, upper or lower plot, including whether to plot the diagonal or not.
      if(plot_type == "full"){ corr <- corr }
      else if(plot_type == "upper") { corr[lower.tri(corr, diag = !plot_diagonal)] <- NA }
      else if(plot_type == "lower") { corr[upper.tri(corr, diag = !plot_diagonal)] <- NA }

      # Format the correlation matrix into a dataframe and add it to the correlation list
      corr_list[[i]] <- corr %>% as.data.frame()  %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(group = i)
    }

    # Bind all the elements in the correlation list into a unique data frame
    corr_list <- dplyr::bind_rows(corr_list)

    # Reshape the correlation dataframe
    corr_list_melt <- corr_list %>% reshape2::melt() %>% dplyr::as_tibble()
    if(!is.null(samples_order)) { corr_list_melt <- corr_list_melt %>% dplyr::mutate(sample = factor(sample, levels = samples_order), variable = factor(variable, levels = samples_order)) }

    # Initialize the plot with facets
    g <- ggplot(corr_list_melt, aes(sample, variable, fill = value)) +
      geom_tile(color = cell_border, na.rm = T) +
      facet_wrap(~group)
  }

  # Force size of the heatmap
  g <- g + ggh4x::force_panelsizes(rows = unit(plot_size, "mm"), cols = unit(plot_size, "mm")) +
    scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0))

  # Theme pubr
  g <- g + theme_pubr(border = plot_border, legend = legend_pos)

  # Setup the legend
  if(legend_pos %in% c("right")) { legend_height <- legend_size[2]; legend_width <- legend_size[1]; legend_title_angle = 270; legend_title_pos = "right" }
  else if(legend_pos %in% c("left")) { legend_height <- legend_size[2]; legend_width <- legend_size[1]; legend_title_angle = 90; legend_title_pos = "left" }
  else if(legend_pos %in% c("bottom")) { legend_size <- rev(legend_size); legend_height <- legend_size[2]; legend_width <- legend_size[1]; legend_title_angle = 0; legend_title_pos = "bottom" }
  else if(legend_pos %in% c("top")) { legend_size <- rev(legend_size); legend_height <- legend_size[2]; legend_width <- legend_size[1]; legend_title_angle = 0; legend_title_pos = "top" }

  # Format the legend
  g <- g + scale_fill_gradient2(limits = legend_limits, midpoint = median(legend_limits),
                                breaks = seq(legend_limits[1],legend_limits[2],legend_breaks_by), na.value = NA,
                                oob = scales::squish, ## to put out of bound values into scale
                                low = plot_colors[1], high = plot_colors[3], mid = plot_colors[2],
                                guide = guide_colorbar(title = legend_title,
                                                       title.position = legend_title_pos,
                                                       title.theme = element_text(angle = legend_title_angle,
                                                                                  hjust = .5, vjust = .5,
                                                                                  size = legend_title_size),
                                                       frame.colour = "black", frame.linewidth = 1.5,
                                                       ticks.colour = NA,
                                                       barheight = unit(legend_height, "mm"),
                                                       barwidth  = unit(legend_width, "mm")))

  # Draw the coefficients
  if(!is.null(coeffs_color)) { g <- g + geom_text(aes(label = round(value, 2)), color = coeffs_color, size = coeffs_size, na.rm = T) }

  # Add title, subtitle and caption
  if(!is.null(plot_title)) { g <- g + ggtitle( label = plot_title, subtitle = plot_subtitle ) }
  if(!is.null(plot_caption)) { g <- g + labs( caption  = plot_caption) }

  # Format theme
  g <- g + theme(plot.title    = element_text(face = title_face, size = title_size, hjust = title_hjust),
                 plot.subtitle = element_text(face = subtitle_face, size = subtitle_size, hjust = title_hjust),
                 plot.caption  = element_text(size = caption_size),
                 axis.title    = element_blank(),
                 axis.ticks    = element_blank(),
                 axis.text     = element_text(size = axis_text_size, colour = axis_text_color),
                 axis.text.x   = element_text(angle = 90, hjust = 1, vjust = .5),
                 panel.border  = element_rect(size = 1.1))

  return(g)
}
