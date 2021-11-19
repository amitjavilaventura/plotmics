# barDEGs()
# ===============================================

#' @title barDEGs
#' @author amitjavilaventura
#'
#' @description
#' It takes a named list with the DE data of different contrasts and
#' draws an horizontal barplot with the up (right) and downregulated (left) genes.
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
#' @param alpha Numeric of length 1. Alpha (transparency) value for the vars. Deefault: 0.5.
#' @param num_size Numeric of length 1. Size of the number of DEGs written in the bars. Default: 3.
#' @param name_size Numeric of length 1. Size of the names of the contrasts  written next to the bars. Default: 3.
#' @param log2FC Numeric of length 1. Threshold for the log2FoldChange to define a gene as differentially expressed. Default: 1.
#' @param pval Numeric of length 1. Threshold for the adjusted p-value to define a gene as differentially expressed. Default: 0.05.
#' @param total Logical of length 1. If TRUE, it counts the total number of features in each elements of deg_list and writes them in the caption. Default: FALSE.
#' @param title Charachter of length 1 or NULL. Title of the plot. Default: NULL.
#' @param subtitle Charachter of lenght 1 or NULL. Subtitle of the plot. Default: NULL.
#' @param prop_test Logical of length 1. Whether to do a `prop.test()` to check whether the proportions of up and downregulated genes are significantly different from the expected proportions of 0.5. Default: FALSE.
#'
#' @export

barDEGs <- function(deg_list,
                    deg_names = names(deg_list),
                    name_pos = "min",
                    xlim = NULL,
                    position_num = 10,
                    xaxis = F,
                    yaxis = F,
                    colors = c("green", "red"),
                    alpha = 0.5,
                    num_size = 3,
                    name_size = 3,
                    log2FC = 1,
                    pval = 0.05,
                    total = FALSE,
                    title = NULL,
                    subtitle = NULL,
                    prop_test = FALSE){

  # Load packages
  require(dplyr)
  require(purrr)
  require(ggplot2)
  require(ggpubr)

  # Check that inputs are OK
  if(!is.list(deg_list)){ stop("'deg_list' must be a named list of data frames.") }
  else if(!is.character(deg_names)){ stop("'deg_names' must be a character vector.") }
  else if(length(deg_list) != length(deg_names)){ stop("'deg_list' and 'deg_names' must have the same lenght.") }
  else if(!name_pos %in% c('min', 'left', 'right', 'none')){ stop("'name_pos' must be one of c('min', 'left', 'right', 'none').") }
  else if(!is.null(xlim) & !is.numeric(xlim)){ stop("'xlim' must be NULL or a numeric vector of length 2.") }
  else if(!is.logical(xaxis) | !is.logical(yaxis)){ stop("Both 'xaxis' and 'yaxis' must be a logical vector of length 1.") }
  else if(!is.character(colors) | length(colors) != 2){ stop("'colors' must be a character vector of length 2 with valid color names/codes.") }
  else if(!is.numeric(alpha) | alpha < 0 | alpha > 1){ stop("'alpha' must be a numeric vector of length 1 with a value between 0 and 1")}

  # Mutate deg_list to be able to change the log2FC and pvalue thresholds for defining the DEGs
  deg_list <- deg_list %>% purrr::map(~dplyr::mutate(.x, DEG = "NS"))
  deg_list <- deg_list %>% purrr::map(~dplyr::mutate(.x, DEG = if_else(log2FoldChange > log2FC & padj < pval, "Upregulated", DEG)))
  deg_list <- deg_list %>% purrr::map(~dplyr::mutate(.x, DEG = if_else(log2FoldChange < -log2FC & padj < pval, "Downregulated", DEG)))

  # Count total elements in the list
  total_counts <- deg_list %>%
    # Add names to each element in the list
    purrr::set_names(deg_names) %>%
    # Add a contrast variable using the name of each element in the list
    purrr::imap(~dplyr::mutate(.x, contrast = .y)) %>%
    # Concatenate the dataframes
    dplyr::bind_rows() %>%
    # Count the number of genes in each df
    dplyr::count(contrast, name = "total_counts")

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
    dplyr::bind_rows() %>%
    # Change the number of downregulated genes to negative
    dplyr::mutate(number = if_else(DEG == "Downregulated", -n, n)) %>%
    # Add a variable for the position of the number of DEGs and the hjust
    dplyr::mutate(pos_num   = if_else(DEG == "Downregulated", -position_num, position_num)) %>%
    dplyr::mutate(hjust_num = if_else(DEG == "Downregulated", 1, 0)) %>%
    # Add a variable for the position of the name of the contrast
    dplyr::mutate(contrast_pos = if_else(DEG == "Downregulated", number-position_num, number+position_num)) %>%
    # Join with total number of elements
    dplyr::left_join(total_counts, by = "contrast")

  # Do prop.test()
  if(prop_test){

    # Set name position to left
    name_pos <- "left"

    # Start the prop.test
    deg_numbers <- deg_numbers %>%
      # Group by contrast
      dplyr::group_by(contrast) %>%
      # Compute the total number of DEGs
      # Compute the proportion of up and downregulated genes (not really necessary)
      # Add the expected proportions (0.5)
      # Do prop.test() for each contrast and add the p-value to the dataframe
      # Since the p-value is the same between up and downregulated genes, add only the p-value for the up.
      dplyr::mutate(sum_de = sum(n),
                    prop_de = n/sum_de,
                    expected_prop = 0.5,
                    p = prop.test(n,sum_de,p=expected_prop)$p.value,
                    p = ifelse(number < 0, NA, p),
                    p = ifelse(!is.na(p), paste("p:", signif(p,2)), NA)) %>%
      dplyr::ungroup()
  }

  # Set at which site will the contrast name be written (if prop_test = T, this is always set to left)
  if(name_pos == "min"){
    deg_numbers <- deg_numbers %>%
      dplyr::group_by(contrast) %>%
      dplyr::mutate(contrast_name = if_else(n == min(n), contrast, NULL)) %>%
      dplyr::ungroup()
  } else if(name_pos == "left"){
    deg_numbers <- deg_numbers %>%
      dplyr::mutate(contrast_name = if_else(n == -number, contrast, NULL))
  } else if(name_pos == "right"){
    deg_numbers <- deg_numbers %>%
      dplyr::mutate(contrast_name = if_else(n == number, contrast, NULL))
  } else if(name_pos == "none"){
    deg_numbers <- deg_numbers %>%
      dplyr::mutate(contrast_name = "")
  }

  # Do the plot
  updown_bar <- deg_numbers %>%
    ggplot(aes(number, contrast, fill = DEG)) +
    geom_col(color = "black", width = 0.5, alpha = alpha) +
    # Add the number of DEGs to each bar
    geom_text(mapping = aes(label = n, x = pos_num, hjust = hjust_num), size = num_size, na.rm = TRUE) +
    # Add the name of the contrast
    geom_text(mapping = aes(label = contrast_name, x = contrast_pos, hjust = hjust_num), size = name_size, na.rm = TRUE, fontface = "bold")  +
    # Change default colors
    scale_fill_manual(values = c(colors[1], colors[2])) +
    # Add vertical line at 0
    geom_vline(aes(xintercept = 0), size = 0.5) +
    # Custom theme
    theme_pubr(border = F, margin = T, legend = "bottom") +
    theme(axis.title = element_blank(),
          legend.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5, face = "italic"),
          plot.caption = element_text(hjust = 1))

  # Add total number of elements
  if(total){
    # If all elements in degs_list have the same number of regions/genes, write the number only one time.
    # If else, write the number of genes in each element
    if(length(unique(total_counts$total_counts)) == 1) { updown_bar <- updown_bar + labs(caption = paste("N = ", unique(total_counts$total_counts))) }
    else { updown_bar + labs(caption = paste(total_counts$contrast, total_counts$total_counts, sep = " = ", collapse = "; ")) }
  }

  if(!xaxis){ updown_bar <- updown_bar + theme(axis.text.x = element_blank()) }
  if(!yaxis){ updown_bar <- updown_bar + theme(axis.text.y = element_blank()) }

  # Change limits of xaxis
  if(!is.null(xlim)){ updown_bar <- updown_bar + coord_cartesian(xlim = c(xlim[1], xlim[2])) }

  # If title == NULL, turn subtitle to NULL, else write title and subtitle
  if(is.null(title)){ subtitle <- NULL }
  else {  updown_bar <- updown_bar + labs(title = title, subtitle = subtitle) }

  # Add p-val from prop.test
  if(prop_test){ updown_bar <- updown_bar + geom_text(aes(label=p, x = number+pos_num), hjust = 0, size = name_size, fontface = "italic", na.rm = T) }

  # Return plot
  return(updown_bar)
}
