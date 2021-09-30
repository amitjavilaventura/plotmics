
# ggVennPeaks() -------------------------

#' @title ggVennPeaks
#' @author amitjavilaventura
#'
#' @description
#' Function that calls 'getVennCounts()' and draws a VennDiagram of
#' peak intersections using the `ggvenn` package.
#'
#' @usage ggVennPeaks(peak_list, peak_names = names(peak_list), percent = T, fill = c("blue", "gold3"), alpha = .4, color = "black", text_color = "black", name_size = 5, label_size = 3, title = "", subtitle = "")
#'
#' @seealso `ggvenn::ggvenn`
#' @seealso `ChIPseeker::makeVennDiagram`
#'
#' @param peak_list List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param peak_names Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param percent Logical of length 1. Whether to show the percentage corresponding to each part of the intersection. Default = T.
#'  Character of the same length as 'conds'. Same values as in 'conds' but with the desired order of priority. Default: 'conds'.
#' @param in_fill Character of length equal to 'peak_list'. Colors of each set in 'peak_list'. Default: c("blue", "gold3").
#' @param alpha Numeric of length 1. Transparency of the the Venn circles. Default: 0.4
#' @param out_color Character of length 1 or equal to 'peak_list'. Color of the circles outer lines. Default: "black"
#' @param name_color Character of length 1 or equal to 'peak_list'. Color of the titles of each set of peaks. Default: "black"
#' @param text_color Character of length 1. Color of the text inside the Venn circles (i.e. number of peaks and percentage in each intersection). Default: "black"
#' @param name_size Numeric of length 1. Size of the name of each peak set. Default: 5
#' @param label_size Numeric of length 1. Size of the text inside the Venn circles. Default: 3
#' @param title Character of length 1. Title of the plot. Can be NULL. Default: ""
#' @param subtitle Character of length 1. Title of the plot. Can be NULL. Default: ""
#'
#'
#' @export

ggVennPeaks <- function(peak_list, peak_names = names(peak_list), percent = T,
                        in_fill = c("blue", "gold3"), alpha = .4,
                        out_color = "black", name_color = "black", text_color = "black",
                        name_size = 5, label_size = 3, title = "", subtitle = ""){

  # Load requireed packages
  require(ggvenn)
  require(dplyr)
  require(reshape2)
  require(magrittr)
  require(purrr)
  require(tidyr)

  # Check that inputs are OK
  if(!is.list(peak_list)){ stop("'peak_list' must be a (named) list of dataframes with the columns 'seqnames', 'start' and 'end'.") }
  if(length(peak_list) != length(peak_names)){ stop("'peak_names' must be a character vector with the same length as 'peak_list'.") }

  # Get Venn Counts and the peaks in each set
  x <- getVennCounts(peaks = peak_list, conds = peak_names)

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

  if(!is.null(title)){
    venn <- venn + labs(title = title, subtitle = subtitle) +
      theme(plot.title = element_text(face = "bold", hjust = .5),
            plot.subtitle = element_text(face = "italic", hjust = .5))
  }

  # Return the Venn diagram.
  return(venn)

}
