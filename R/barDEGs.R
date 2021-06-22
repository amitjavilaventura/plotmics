# barDEGs()
# ===============================================

#' @title barDEGs
#' @author amitjavilaventura
#'
#' It takes a named list with the DE data of different contrasts and draws an horizontal barplot with the up (right) and downregulated (left) genes.
#'
#' @returns A ggplot2-based barplot with the upregulated and downregulated genes in each contrast.
#'
#' @param deg_list Named list of dataframes with, at least, a column named DEG with the values 'Downregulated', "Upregulated" and "NS"
#' @param deg_names Character vector of equal length to 'deg_list' names of the elements in 'deg_list'. Default: names(deg_list)
#' @param name_pos Character of length 1. Where to write the name of each contrast. One of c('min', 'left', 'right', 'none'). If 'min', the name of the contrast will be written at the side where there are less DE genes; if 'right', the name of the contrast will be written to the right side of the plot; if 'left', the name of the contrast will be written to the left side of the plot; if 'none' the name of the contrast won't be written. Default: 'min'.
#' @param xlim Numeric of length 2 or NULL. The limits of the X axis. If NULL, the limits of the X axis are the default. Useful when 'name_pos' is 'left' or 'right'. Default: NULL
#' @param position_num Numeric of length 1. The position of the number of DEGs in each bar and the position of the contrast label. Default = 10.
#' @param xaxis Logical of length 1. Whether to draw the text of the X axis (number of DE) or not. Default: FALSE.
#' @param yaxis Logical of length 1. Whether to draw the text of the Y axis (contrast names) or not. Default: FALSE.
#' @param colors Character of length 2. Colors for the downregulated and upregulated genes. Default: c("green", "red")
#' @param alpha Numeric of length 1. Alpha (transparency) value for the vars. Deefault: 0.5
#'
#' @export

barDEGs <- function(deg_list, deg_names = names(deg_list),
                       name_pos = "min", xlim = NULL, position_num = 10,
                       xaxis = F, yaxis = F,
                       colors = c("green", "red"), alpha = 0.5){

  # Load packages
  require(dplyr)
  require(purrr)
  require(ggplot2)

  # Check that inputs are OK
  if(!is.list(deg_list)){ stop("'deg_list' must be a named list of data frames.") }
  else if(!is.character(deg_names)){ stop("'deg_names' must be a character vector.") }
  else if(length(deg_list) != length(deg_names)){ stop("'deg_list' and 'deg_names' must have the same lenght.") }
  else if(!name_pos %in% c('min', 'left', 'right', 'none')){ stop("'name_pos' must be one of c('min', 'left', 'right', 'none').") }
  else if(!is.null(xlim) & !is.numeric(xlim)){ stop("'xlim' must be NULL or a numeric vector of length 2.") }
  else if(!is.logical(xaxis) | !is.logical(yaxis)){ stop("Both 'xaxis' and 'yaxis' must be a logical vector of length 1.") }
  else if(!is.character(colors) | length(colors) != 2){ stop("'colors' must be a character vector of length 2 with valid color names/codes.") }
  else if(!is.numeric(alpha) | alpha < 0 | alpha > 1){ stop("'alpha' must be a numeric vector of length 1 with a value between 0 and 1")}

  # Named list of DEGs
  deg_numbers <- deg_list %>%
    # Add names to each element in the list
    purrr::set_names(deg_names) %>%
    # Count number of upregulated and downregulated genes in each contrast
    purrr::map(~dplyr::count(.x, DEG)) %>%
    purrr::map(~dplyr::filter(.x, DEG != "NS")) %>%
    # Add a contrast variable using the name of each element in the list
    purrr::imap(~dplyr::mutate(.x, contrast = .y)) %>%
    # Bind all dataframes in one
    bind_rows() %>%
    # Change the number of downregulated genes to negative
    mutate(number = if_else(DEG == "Downregulated", -n, n)) %>%
    # Add a variable for the position of the number of DEGs and the hjust
    mutate(pos_num   = if_else(DEG == "Downregulated", -position_num, position_num)) %>%
    mutate(hjust_num = if_else(DEG == "Downregulated", 1, 0)) %>%
    # Add a variable for the position of the name of the contrast
    mutate(contrast_pos = if_else(DEG == "Downregulated", number-position_num, number+position_num))

  # Set at which site will the contrast name be written
  if(name_pos == "min"){
    deg_numbers <- deg_numbers %>%
      group_by(contrast) %>%
      mutate(contrast_name = if_else(n == min(n), contrast, NULL)) %>%
      ungroup()
  } else if(name_pos == "left"){
    deg_numbers <- deg_numbers %>%
      mutate(contrast_name = if_else(n == -number, contrast, NULL))
  } else if(name_pos == "right"){
    deg_numbers <- deg_numbers %>%
      mutate(contrast_name = if_else(n == number, contrast, NULL))
  } else if(name_pos == "none"){
    deg_numbers <- deg_numbers %>%
      mutate(contrast_name = "")
  }

  # Do the plot
  updown_bar <- deg_numbers %>%
    ggplot(aes(number, contrast, fill = DEG)) +
    geom_col(color = "black", width = 0.5, alpha = alpha) +
    # Add the number of DEGs to each bar
    geom_text(mapping = aes(label = n, x = pos_num, hjust = hjust_num), size = 3) +
    # Add the name of the contrast
    geom_text(mapping = aes(label = contrast_name, x = contrast_pos, hjust = hjust_num), size = 3)  +
    # Change default colors
    scale_fill_manual(values = c(colors[1], colors[2])) +
    # Add vertical line at 0
    geom_vline(aes(xintercept = 0), size = 0.5) +
    # Custom theme
    theme_pubr(border = F, margin = T, legend = "bottom") +
    theme(axis.title = element_blank(),
          legend.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())

  if(!xaxis){ updown_bar <- updown_bar + theme(axis.text.x = element_blank()) }
  if(!yaxis){ updown_bar <- updown_bar + theme(axis.text.y = element_blank()) }

  # Change limits of xaxis
  if(!is.null(xlim)){updown_bar <- updown_bar + coord_cartesian(xlim = c(xlim[1], xlim[2]))}

  # Return plot
  return(updown_bar)
}
