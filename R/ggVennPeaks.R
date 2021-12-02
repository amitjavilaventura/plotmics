
# ggVennPeaks() -------------------------

#' @title ggVennPeaks
#' @author amitjavilaventura
#'
#' @description
#' Function that calls 'getVennCounts()' and draws a VennDiagram of
#' peak intersections using the `ggvenn` package.
#'
#' @seealso `ggvenn::ggvenn`
#' @seealso `ChIPseeker::makeVennDiagram`
#'
#' @param peak_list List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param peak_names Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param percent Logical of length 1. Whether to show the percentage corresponding to each part of the intersection. Default = T.
#' @param stranded Logical of length 1. If TRUE, it considers the strand information to do the overlaps. Default = FALSE.
#' @param true_overlaps Logical of length 1. Only if peak_list has length 2. If TRUE, it counts the overlapping regions from each set using plyranges. If TRUE, stranded and percent are set to FALSE. Default: FALSE
#' @param in_fill Character of length equal to 'peak_list'. Colors of each set in 'peak_list'. Default: c("blue", "gold3").
#' @param alpha Numeric of length 1. Transparency of the the Venn circles. Default: 0.4.
#' @param out_color Character of length 1 or equal to 'peak_list'. Color of the circles outer lines. Default: "black".
#' @param name_color Character of length 1 or equal to 'peak_list'. Color of the titles of each set of peaks. Default: "black".
#' @param text_color Character of length 1. Color of the text inside the Venn circles (i.e. number of peaks and percentage in each intersection). Default: "black".
#' @param name_size Numeric of length 1. Size of the name of each peak set. Default: 5.
#' @param label_size Numeric of length 1. Size of the text inside the Venn circles. Default: 3.
#' @param title Character of length 1. Title of the plot. Can be NULL. Default: "".
#' @param subtitle Character of length 1. Title of the plot. Can be NULL. Default: "".
#' @param return_peaks Logical of length 1. Whether to return the peaks of each set overlapping with the other. Only if peak_list has length 2. If TRUE, it won't draw the plot. Default: FALSE.
#'
#' @export

ggVennPeaks <- function(peak_list,
                        peak_names = names(peak_list),
                        percent = TRUE,
                        stranded = FALSE,
                        true_overlaps = FALSE,
                        in_fill = c("blue", "gold3"),
                        alpha = .4,
                        out_color = "black",
                        name_color = "black",
                        text_color = "black",
                        name_size = 5,
                        label_size = 3,
                        title = "",
                        subtitle = "",
                        return_peaks = FALSE){

  # Load requireed packages
  suppressPackageStartupMessages(require(ggvenn))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(reshape2))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(purrr))
  suppressPackageStartupMessages(require(tidyr))
  suppressPackageStartupMessages(require(plyranges))

  # Check that inputs are OK
  if(!is.list(peak_list)){ stop("'peak_list' must be a (named) list of dataframes with, at least, the columns 'seqnames', 'start' and 'end'.") }
  if(length(peak_list) != length(peak_names)){ stop("'peak_names' must be a character vector with the same length as 'peak_list'.") }
  if(stranded){
    if(!"strand" %in% colnames(bind_rows(peak_list))){ stop("If 'stranded' is TRUE, a 'strand' column must be present in all the elements from peak list.") }
    else { peak_list <- peak_list %>% purrr::map(~dplyr::mutate(.x, strand = if_else(strand %in% c("\\.", "\\*", ".", "*"), "*", strand))) }
  }

  # Check true number of overlaps for each region; convert stranded and percent to FALSE
  if(true_overlaps & length(peak_list)==2) { stranded <- FALSE; percent  <- FALSE }

  # Change strand values to accepted values (e.g. . to *)
  peak_list <- peak_list %>% purrr::map(~dplyr::mutate(.x, strand = if_else(strand %in% c(".", "*", "\\.", "\\*"), "*", strand)))

  # Compute true number of overlaps and get overlapping peaks
  peaks1 <- peak_list[[1]] %>% plyranges::as_granges();  peaks2 <- peak_list[[2]] %>% plyranges::as_granges()
  overlaps1 <- plyranges::filter_by_overlaps(peaks1, peaks2); overlaps2 <- plyranges::filter_by_overlaps(peaks2, peaks1)

    # If desired, return overlapping peaks
  if(return_peaks & length(peak_list)==2){
    # List of overlaps
    overlaps <- list(overlaps1, overlaps2) %>% purrr::set_names(peak_names) %>% purrr::map(~as.data.frame(.x))
    # Return overlapping peaks
    return(overlaps)
  }

  # Get Venn Counts and the peaks in each set
  x <- getVennCounts(peaks = peak_list, conds = peak_names, stranded = stranded)


  # Set default colors in case the number of specified colors
  #   does not match the number elements in 'peak_list'
  if(length(in_fill) != length(peak_list)) {
    if(length(peak_list) == 2){ in_fill <- c("blue", "gold3") }
    else if(length(peak_list) == 3){ in_fill <- c("blue", "gold3", "pink") }
    else if(length(peak_list) == 4){ in_fill <- c("blue", "gold3", "pink", "green") }
    else if(length(peak_list) == 5){ in_fill <- c("blue", "gold3", "pink", "green", "orange") }
  }

  # Transform the matrix of the peaks in each set in order to plot the Venn
  y <- x$matrix %>%
    as_tibble() %>%
    magrittr::set_colnames(c("Peak", paste("cond", 1:length(peak_list), sep = ""))) %>%
    reshape2::melt() %>% mutate(value = if_else(value == 1, Peak, NULL)) %>%
    tidyr::pivot_wider(names_from = "variable", values_from = "value") %>%
    dplyr::select(-Peak) %>%
    dplyr::as_tibble() %>%
    as.list() %>%
    purrr::set_names(peak_names) %>%
    purrr::map(~na.omit(.x))

  # Draw the Venn diagram with ggvenn
  venn <- ggvenn::ggvenn(data = y, show_percentage = percent,
                         fill_color = in_fill, fill_alpha = alpha,
                         stroke_color =  out_color,
                         set_name_color = name_color, set_name_size = name_size,
                         text_color = text_color, text_size = label_size) +
    scale_y_discrete(expand = c(0,0.5))

  if(length(peak_list) %in% 4:5) { venn <- venn + scale_x_discrete(expand = c(0,0.5))}


  # Add titles
  if(!is.null(title)){
    venn <- venn + labs(title = title, subtitle = subtitle) +
      theme(plot.title = element_text(face = "bold", hjust = .5),
            plot.subtitle = element_text(face = "italic", hjust = .5))
  }

  # Add number or true overlaps for each set of regions
  if(true_overlaps & length(peak_list)==2) {

    # Load required packages
    suppressPackageStartupMessages(require(colorspace))

    # Count the numbers
    num1 <- length(overlaps1)
    num2 <- length(overlaps2)

    # Annotate true numbers
    venn <- venn +
      annotate("text", 0, 0.3, label = as.character(num1), color = colorspace::darken(in_fill[1], amount = .5), fontface = "bold", size = label_size*0.8) +
      annotate("text", 0, -0.3, label = as.character(num2), color = colorspace::darken(in_fill[2], amount = .5), fontface = "bold", size = label_size*0.8)
  }

  # Return the Venn diagram.
  return(venn)

}
