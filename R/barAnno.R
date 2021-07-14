### barAnno()
### ===============================================================================================

#' @title barAnno
#' @author amitjavilaventura
#'
#' Function for ChIP-seq and ATAC-seq.
#' It must be used after the function annotatePeak() from the R package ChIPseeker. @seealso \code{\link{annotatePeak}}
#'
#' It takes a list of annotation objects that come as output of annotatePeak() and changes the features to "Promoter", "Distal" and "Gene body" (or to "Promoter" and "Distal"). Finally it plots a bargraph with the distribution of all the proportions
#' As a ggplot2-based function, it allows to add more layers to format the plot.
#'
#' @param anno_list Named list of annotation objects that come from annotatePeak().
#' @param anno_names Charachter vector of the same length as 'anno_list'. Names that will be given to each of the objects in anno_list. Not that will be the names plotted in the bargraph
#' @param names_order Character vector with the same entries as 'names' with the order wanted to plot the data. Defaults to: unique(names).
#' @param protein Character vector of the same length as 'anno_list' with the protein chipped in each dataframe of anno_list. This names will be passed through `facet_wrap()`. Only works if facet is set to TRUE. Default: NULL.
#' @param protein_order Character vector with the same entries as 'protein' with the order wanted to plot the data. Defaults to: unique(protein).
#' @param main Character of lenght 1. Title of the plot. Default: NULL.
#' @param subtitle Character of lenght 1. Subtitle of the plot. Default: NULL.
#' @param ylab Character of lenght 1. Title of the Y axis. Default: NULL.
#' @param xlab Character of lenght 1. Title of the X axis. Default: NULL
#' @param palette Character of lenght 1. Color palette used to color the bars through the function `scale_fill_brewer()`. Default: "Set1".
#' @param legend_position Character of lenght 1. Position of the legend. One of c("none", "bottom", "right", "left," "top"). Default: "right"
#' @param anno_num Numerical or character of length 1. Number of annotations to plot, either 2 (Promoter/Distal), 3 (Promoter/Gene body/Distal) or 'all' (the annotatePeak() default). Default: 2.
#' @param fill_position Logical. If TRUE (default), it the plotted bars will represent proportion of peaks in each feature. If FALSE, the bars will have the height of the total number of peaks with the correspondent feature color.
#' @param xangle Numerical of length 1. Angle of the text in the X axis. Default: 20.
#' @param width Numerical of length 1. Width of the bar in relative units. Default: 0.6
#'
#' @export

barAnno  <- function(anno_list,
                     anno_names = names(anno_list), names_order = unique(anno_names),
                     protein = NULL, protein_order = unique(protein),
                     main = NULL, subtitle = NULL, ylab = NULL, xlab = NULL,
                     color_palette = "Set2", legend_position = "right", anno_num = 2,
                     fill_position = T, xangle = 20, width = 0.6){

  # Load packages
  library(dplyr)
  library(purrr)
  library(magrittr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)



  # Check that inputs are OK.
  if(!is.list(anno_list)) { stop("'anno_list' must be a (named) list of data frames with the genomic coordinates of peak and its genomic annotation. Alternatively, elements in 'anno_list' can be of the class 'csAnno' (output of annotatePeak()).") }
  else if(!is.null(anno_names) & (length(anno_names) != length(anno_list))){ stop("'anno_names' must be a not NULL character vector with the same length of 'anno_list'.") }
  else if(!all(names_order %in% anno_names)) { stop("All the elements of 'names_order' must be inside 'names'") }
  else if(!is.null(protein) & (length(protein) != length(anno_list))){ stop("'protein' must be a character vector with the same length of 'anno_list' or NULL.") }

  # Set the names of each annotatePeak object in the list with the vector anno_names
  anno <- purrr::set_names(x = anno_list, nm = anno_names) %>%

    # Convert to tibble
    purrr::map(~as_tibble(.x)) %>%

    # Reformat annotation. Just change the "Promoter ..." annotation to "Promoter".
    purrr::map(~dplyr::mutate(.x, annotation = annotation %>% gsub(" \\(.*$", "", .))) %>%

    # Write an extra column to each dataframe with the name of the dataframe (provided in anno_names)
    purrr::imap(~dplyr::mutate(.data = .x, condition = .y)) %>%

    # Select only the annotation and condition columns
    purrr::map(~dplyr::select(.data = .x, annotation, condition))

  if(!is.null(protein)){
    # Set the names to protein
    anno <- purrr::set_names(x = anno, nm = protein) %>%

      # Write an extra column to each dataframe with the name of the dataframe (provided in protein)
      purrr::imap(~dplyr::mutate(.data = .x, prot = .y))
  }

  # Bind dataframes by rows
  anno_df <- anno %>% bind_rows() %>%
    # Convert the condition to factor and give it an order
    dplyr::mutate(condition = factor(condition, levels = names_order))

  # Convert the protein to factor and give it an order
  if(!is.null(protein)){ anno_df <- anno_df %>% dplyr::mutate(prot = factor(prot, levels = protein_order)) }


  # Rewrite annotation as promoter/distal,  promoter/distal/gene body or don't change it
  if(anno_num == "all"){ anno_df <- anno_df }
  else if(anno_num == 2){ anno_df <- anno_df %>% dplyr::mutate(annotation = if_else(annotation == "Promoter", "Promoter", "Distal")) }
  else if(anno_num == 3){
    anno_df <- anno_df %>%
      dplyr::mutate(annotation = dplyr::recode(annotation,
                                               "Distal Intergenic" = "Distal", "Downstream" = "Distal",
                                               "Intron" = "Gene body", "Exon" = "Gene body",
                                               "5' UTR" = "Gene body", "3' UTR" = "Gene body"))
  }


  # Start ggplot -----
  g <- ggplot(anno_df, aes(condition, fill = annotation))

  # Create gglpot2-based barplot
  if(fill_position & !is.null(protein)){ g <- g + geom_bar(position = "fill", color = "Black", width = width) + facet_wrap(~prot) }
  else if(fill_position & is.null(protein)){ g <- g + geom_bar(position = "fill", color = "Black", width = width) }
  else if(!is.null(protein)){ g <- g + geom_bar(color = "Black", width = width) + facet_wrap(~prot) }
  else{ g <- g + geom_bar(color = "Black", width = width) }


  if(length(color_palette)==1) { scale_colors <- scale_fill_brewer(palette = color_palette) }
  else if(length(color_palette)>1) { scale_colors <- scale_fill_manual(values = color_palette) }

  # Write plot title, subtitle and axis labels
  g <- g + ggtitle(main, subtitle) + ylab(ylab) + xlab(xlab) +

    # Format colors with ggplot2 function scale_fill_brewer
    scale_colors +

    # Basic formatting
    theme_pubr(legend = legend_position, x.text.angle = xangle, border = T) +
    theme(legend.title = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          axis.title = element_text(face = "bold"))

  # Return plot
  return(g)
}

