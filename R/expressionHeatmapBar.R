# expressionHeatmapBar()
# ===============================================

#' @title expressionHeatmapBar
#' @author amitjavilaventura
#'
#' @description
#' It takes 3 data.frames with expression values regardless the strand, sense or antisense (i.e., featureCounts -s 0|1|2)
#' and draws a heatmap with the expression values and an hotizontal barplot with the expression sense and antisense.
#'
#' @returns A patchwork combination of a ggplot-based heatmap and a ggplot-based barplot (or a list with the individual plots).
#'
#' @param expr Dataframe with a Geneid column and expression values or NULL. Column names must be equal to 'expr_sense' and 'expr_antisense'.
#' @param expr_sense Dataframe with a Geneid column and expression values . Column names must be equal to 'expr' and 'expr_antisense'.
#' @param expr_antisense Dataframe with a Geneid column and expression values . Column names must be equal to 'expr' and 'expr_sense'.
#' @param cols Character vector with the names of the columns to take into account (Geneid does not count) or NULL. If NULL, all columns are taken into account. Default: NULL
#' @param genes Character vector with the names of the genes to take into account or NULL. If NULL, all genes are taken into account. Default: NULL
#' @param hm_values Logical of length 1. If TRUE (the default), write the corresponding expression values in each cell of the heatmap
#' @param hm_value_round Numerical of length 1. The number of decimals the expression values will be rounded to. Default: 0.
#' @param hm_value_size Numerical of length 1. Size of the expression values written in the heatmap. Default: 3.
#' @param hm_value_color Charachter of length 1. Color of the expression values written in the heatmap. Default: "black".
#' @param hm_legend_colors Character of length 2. Colors for the lower and upper limits of the heatmap, respectively. Passed through scale_fill_gradient. Default: c("white", "cornflowerblue"),
#' @param hm_legend_limits Numerical of length 2. Lower and upper limits of the heatmap, respectively. Passed through scale_fill_gradient. Default: c(0,1e6),
#' @param hm_legend_breaks Numerical vector. Breaks shown in the legend of the heatmap. Passed through scale_fill_gradient. Default: seq(0,1e6, 1e5),
#' @param hm_legend_title Character of length 1 or NULL. Title of the legend of the heatmap. Default: "Counts"
#' @param hm_legend_size Numerical of length 2. Height and width (in mm) of the legend of the heatmap, respectively. Passed through scale_fill_gradient. Default: c(40, 7),
#' @param bar_pos Character of length 1. Either "fill" or "stack", Position of the columns, passed through the argument position in `ggplot2::geom_col()`. Default: "fill".
#' @param bar_xlab Character of length 1 or NULL. Title of the X axis of the barplot. Default: NULL
#' @param bar_colors Character of length 2. Colors for the antisense and sense bars in the barplot. Default: c("green", "tomato")
#' @param combine_plots Logical of length 1. If TRUE (the default), it combines the heatmap and the barplot using `patchwork`. If FALSE, it returns a list with the individual plots.
#' @param plot_title Character of length 1 or NULL. Title of the plot. Only if combine_plots = T. Default: NULL.
#'
#' @export

expressionHeatmapBar <- function(expr,
                                 expr_sense,
                                 expr_antisense,
                                 cols             = NULL,
                                 genes            = NULL,
                                 hm_values        = T,
                                 hm_value_round   = 0,
                                 hm_value_size    = 3,
                                 hm_value_color   = "black",
                                 hm_legend_colors = c("white", "cornflowerblue"),
                                 hm_legend_limits = c(0,1e6),
                                 hm_legend_breaks = seq(0,1e6, 1e5),
                                 hm_legend_title  = "Counts",
                                 hm_legend_size   = c(40, 7),
                                 bar_pos          = "fill",
                                 bar_xlab         = NULL,
                                 bar_colors       = c("green", "tomato"),
                                 combine_plots    = T,
                                 plot_title       = NULL){

  # Load packages
  require(dplyr)
  require(reshape2)
  require(purrr)
  require(ggh4x)
  require(ggplot2)
  require(ggmitji)

  # Check that inputs are OK
  if(!is.null(expr) & !is.data.frame(expr)) { stop("'expr' should be NULL or a data.frame with Geneid and expression values.") }
  if(!is.data.frame(expr_sense)) { stop("'expr_sense' should be a data.frame with Geneid and expression values.") }
  if(!is.data.frame(expr_antisense)) { stop("'expr_antisense' should be a data.frame with Geneid and expression values.") }
  if(!"Geneid" %in% colnames(expr_sense)){ stop("The column names of the expression data.frames should be equal between them and have a Geneid column") }
  if(is.null(cols)) { cols <- colnames(expr)[which(colnames(expr) != "Geneid")] }
  if(is.null(expr)) {
    require(tibble)
    expr <- (expr_sense %>%  dplyr::arrange(Geneid) %>% tibble::column_to_rownames("Geneid")) + (expr_antisense %>%  dplyr::arrange(Geneid) %>% tibble::column_to_rownames("Geneid"))
  }
  if(any(colnames(expr) != colnames(expr_sense))) { stop("The column names of the expression data.frames should be equal between them and have a Geneid column") }
  if(any(colnames(expr) != colnames(expr_antisense))) { stop("The column names of the expression data.frames should be equal between them and have a Geneid column") }

  # Heatmap - ALL reads regardless the strand -----
  ## Format expression dataframe, select desired cols and reshape it
  expr_df <- expr %>%
    dplyr::select(Geneid, cols) %>%
    reshape2::melt()

  ## Retain only desired genes
  if(!is.null(genes)) { expr_df <- expr_df %>% dplyr::filter(Geneid %in% genes)}

  ## Draw heatmap
  expr_heatmap <- expr_df %>%
    ggplot(aes(variable, Geneid, fill = value)) +
    geom_tile(color = "gray") +
    # Split rows by ID, set scales and space to free
    ggh4x::facet_grid2(rows = vars(Geneid), scales = "free_y", space = "free") +
    # Format the legend (colorbar)
    scale_fill_gradient(low = hm_legend_colors[1], high = hm_legend_colors[2], limits = hm_legend_limits, breaks = hm_legend_breaks, oob = scales::squish,
                        guide = guide_colorbar(title = hm_legend_title, title.position = "left", frame.colour = "black", barheight = unit(hm_legend_size[1], "mm"),
                                               barwidth = unit(hm_legend_size[2], "mm"), frame.linewidth = 1.5, ticks.colour = NA,
                                               title.theme = element_text(angle = 90, hjust = 0.5,  vjust = 0.5, size = 8))) +
    # Format theme, remove facet strips, remove whitespace between bars and panel border, and remove axis labels
    ggmitji::theme_custom(legend = "left", x.text.angle = 90, x.text.hjust = 1, x.text.vjust = .5, title.face = "italic", title.size = 9) +
    scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
    labs(x = NULL, y = NULL)

  ## Write values in each cell
  if(hm_values){  expr_heatmap <- expr_heatmap + geom_text(aes(label = round(value, hm_value_round)), size = hm_value_size, color = hm_value_color) }

  # Barplot - Proportion of sense and antisense counts -----
  ## Format the counts data frames with sense and antisense
  ## Make list and set names to the list
  expr_strand <- list(expr_sense, expr_antisense) %>%
    purrr::set_names(c("Sense", "Antisense")) %>%
    # Retain the desired columns and create a column for sense antisense
    purrr::imap(~dplyr::select(.x, Geneid, cols) %>% dplyr::mutate(sense = .y)) %>%
    # Bind lists by rows
    # Reshape the resulting dataframe
    dplyr::bind_rows() %>%
    unique() %>%
    reshape2::melt()

  ## Retain only desired genes
  if(!is.null(genes)) { expr_strand <- expr_strand %>% dplyr::filter(Geneid %in% genes)}

  ## Draw horizontal barplot
  expr_strand_barplot <-  expr_strand %>%
    ggplot(aes(y = variable, x = value, fill = sense)) +
    geom_col(position = bar_pos, color = "black") +
    # Split rows by ID, then put the strips at the left side
    ggh4x::facet_grid2(rows = vars(Geneid), switch = "y") +
    # Customize labs, theme, position of the Y axis and remove strips
    labs(x = bar_xlab, y = NULL, title = NULL) +
    ggmitji::theme_custom(legend = "right", x.text.angle = 90, x.text.hjust = 1, x.text.vjust = .5, title.face = "italic", axis.title.face = "plain", vgrid.major = .3, vgrid.minor = .3) +
    scale_y_discrete(guide = guide_axis(position = "right")) +
    scale_fill_manual(bar_colors)

  # Combine plots or not -----
  if(combine_plots){
    # Load patchwork
    require(patchwork)

    # Combine plots
    final_plot <-
    expr_heatmap+ggmitji::rm_strips() +
      expr_strand_barplot+ggmitji::rm_strips() +
      patchwork::plot_layout(design = "AAAB") +
      patchwork::plot_annotation(title = plot_title, theme = theme(plot.title = element_text(face = "italic",  size = 9, hjust = .5)))

    # Return plot
    return(final_plot)

  } else {
    # Return list of individual plots
    return(list(expr_heatmap, expr_strand_barplot))
  }

}
