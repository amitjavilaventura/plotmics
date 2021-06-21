getVennCounts <- function(peaks, conds = names(peaks), conds_order = conds, plot = F){

  require(pkgcond)
  require(purrr) %>% suppress_messages() %>% suppress_warnings()
  require(dplyr) %>% suppress_messages() %>% suppress_warnings()
  require(plyranges) %>% suppress_messages() %>% suppress_warnings()
  require(magrittr) %>% suppress_messages() %>% suppress_warnings()
  require(ChIPpeakAnno) %>% suppress_messages() %>% suppress_warnings()

  if(!is.list(peaks)){ stop("'peaks' must be a (named) list of dataframes with the columns 'seqnames', 'start' and 'end'.") }
  else if(is.null(conds)){ stop("'conds' must a not-NULL character vector with the conditions of the data frames in 'peaks.") }
  else if(length(peaks) != length(conds)){ stop("'peaks' and 'conds' must have the same length.") }

  len <- length(peaks)

  peaks <- peaks %>% set_names(nm = conds)

  overlaps <- peaks[conds_order] %>%
    purrr::map(~as_granges(.x)) %>%
    makeVennDiagram(plot = plot) %>%
    suppress_messages() %>% suppress_warnings()

  overlaps <- overlaps$vennCounts

  matrix <- matrix(data = rep(0, len), ncol = len, byrow=T) %>% set_colnames(conds_order)
  for(row in 1:nrow(overlaps)){

    counts <- overlaps[row, len+1]

    m <- matrix(rep(overlaps[row, 1:len], counts), ncol = len, byrow = T)

    matrix <- rbind(matrix, m)

  }

  x <- matrix %>%
    na.omit %>%
    as.data.frame() %>%
    dplyr::mutate(rowSum = rowSums(.)) %>%
    dplyr::filter(rowSum != 0) %>%
    dplyr::mutate(peak = paste("peak", 1:nrow(.), sep = "")) %>%
    dplyr::select(peak, everything(), -rowSum)

  return(list("matrix" = x, "vennCounts" = overlaps))

}

# upsetPeaks() -------------------------

#' @title upsetPeaks
#' @author amitjavilaventura
#'
#' @seealso plyranges
#' @seealso UpSetR
#'
#' @usage upsetPeaks(peaks, names = names(peaks), names_order = names)
#'
#' @param peaks List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param conds Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param conds_order Character of the same length as 'conds'. Same values as in 'conds' but with the desired order of priority. Default: 'conds'.
#' @param title (Not available yet) Character of length 1 or NULL. Title of the plot shown at the top. Passed through grid.arrange(). If NULL, title is not shown. Default: NULL
#' @param order.by Same as in UpSetR::upset(). One of "freq", "degree" or both.
#' @param mainbar.y.label Same as in UpSetR::upset().
#' @param sets.x.label Same as in UpSetR::upset().
#' @param ... Further arguments to be passed through 'UpSetR::upset()'.
#'
#' @export

upsetPeaks <- function(peaks, conds = names(peaks), conds_order = conds, order.by = "freq",
                       mainbar.y.label = "Intersect size", sets.x.label = "Set size",
                       ...){

  # Load required packages
  require(UpSetR)
  require(grid)
  require(gridExtra)

  #
  # call makeVenn4upSet to retrieve the peak matrix
  x <- getVennCounts(peaks, conds, conds_order, plot = F)

  # draw upset plot with options
  upset <- UpSetR::upset(data = x$matrix, order.by = order.by, sets.x.label = sets.x.label, mainbar.y.label = mainbar.y.label)
  return(upset)
}

# ggUpsetPeaks() -------------------------

#' @title ggUpsetPeaks
#' @author amitjavilaventura
#'
#' @seealso plyranges
#'
#' @usage ggUpsetPeaks(peaks, conds = names(peaks), conds_order = cond)
#'
#' @param peaks List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param conds Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param conds_order Character of the same length as 'conds'. Same values as in 'conds' but with the desired order of priority. Default: 'conds'.
#' @param num_size Numeric of length 1. Size of the labels indicating the size of each intersection. Default: 4.
#' @param order_by_freq Logical of length 1. If TRUE (the default), the intersections are arranged by descending order of size. Default: TRUE.
#' @param title Character of length 1. Title of the plot. Default: NULL.
#' @param subtitle Character of length 1. Subtitle of the plot. Default: NULL.
#' @param caption Character of length 1. Caption written at the bottom of the plot. Default: NULL.
#'
#' @export

ggUpsetPeaks <- function(peaks, conds = names(peaks), conds_order = conds,
                         order_by_freq = T, num_size = 4,
                         title = NULL, subtitle = NULL, caption = NULL){

  # Load required packages
  require(dplyr)
  require(purrr)
  require(reshape2)
  require(magrittr)
  require(patchwork)
  require(ggplot2)
  require(ggpubr)
  require(ggforestplot)

  # Calculate number of peaks for each condition
  peakNum <- peaks %>% purrr::set_names(conds) %>%
    purrr::imap(~dplyr::mutate(.x, set = .y)) %>%
    bind_rows() %>%
    group_by(set) %>% summarize(n = n()) %>% ungroup()

  # Call makeVenn4upSet to retrieve the count matrix
  x <- getVennCounts(peaks, conds, conds_order, plot = F)

  # Transform count matrix to a data frame
  countMatrix <- x$vennCounts %>%
    matrix(ncol = ncol(x$vennCounts), nrow = nrow(x$vennCounts)) %>%
    as.data.frame() %>%
    set_colnames(colnames(x$vennCounts)) %>%
    set_rownames(rownames(x$vennCounts)) %>%
    filter_all(any_vars(. != 0)) %>%
    tibble::rownames_to_column("Intersection")

  # Arrange by size of intersection
  if(order_by_freq){
    orderInt <- countMatrix %>% dplyr::arrange(desc(Counts)) %>% pull(Intersection)
    countMatrix <- countMatrix %>% dplyr::mutate(Intersection = factor(Intersection, levels = orderInt))
  }

  # Draw the barplot of the intersections
  barPlot <- ggplot(countMatrix, aes(Intersection, Counts)) +
    geom_col(fill = "black", width = 0.5) +
    geom_text(aes(label = Counts, y = Counts + max(Counts)*0.04), size = num_size) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(countMatrix$Counts)*1.1)) +
    ylab("Intersection size") +
    theme_pubclean(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          plot.title = element_text(face = "bold"))

  # Draw barplot with size of each set
  setSize <- peakNum %>%
    ggplot(aes(n, set)) +
    geom_col(fill = "black", width = 0.6) +
    ggforestplot::geom_stripes(odd = "#22222222", even = "#00000000") +
    scale_x_reverse(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    theme_pubclean(base_size = 10, flip = T) +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x = element_line(color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          plot.title = element_text(face = "bold")) +
    xlab("Set size")

  # Draw the points for the combination of sets in each intersection
  setPoints <- countMatrix %>%
    dplyr::select(-Counts) %>%
    reshape2::melt() %>%
    filter(value != 0) %>%
    ggplot(aes(y = variable, x = Intersection)) +
    ggforestplot::geom_stripes(odd = "#22222222", even = "#00000000") +
    geom_point(size = 4) +
    geom_line(aes(group = Intersection), size = 1.5) +
    xlab("Intersections") +
    theme_pubclean(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          plot.title = element_text(face = "bold"),
          axis.line = element_line(color = "black"))


  # Create layout to use with patchwork
  upsetLayout <- "
  ##AAAAA
  ##AAAAA
  ##AAAAA
  ##AAAAA
  BBCCCCC
  "

  # Build upset plot
  upsetPlot <- barPlot + setSize + setPoints +
    plot_layout(design = upsetLayout) +
    plot_annotation(title = title, subtitle = subtitle, caption = caption) &
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"))

  # Return plot
  return(upsetPlot)

}


# ggVennPeaks() -------------------------

#' @title ggVennPeaks
#' @author amitjavilaventura
#'
#' @seealso ggvenn
#' @seealso makeVennDiagram
#'
#' @usage ggVennPeaks(peak_list, peak_names = names(peak_list), percent = T, fill = c("blue", "gold3"), alpha = .4, color = "black", text_color = "black", name_size = 5, label_size = 3
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
#'
#' @export

ggVennPeaks <- function(peak_list, peak_names = names(peak_list), percent = T,
                        in_fill = c("blue", "gold3"), alpha = .4,
                        out_color = "black", name_color = "black", text_color = "black",
                        name_size = 5, label_size = 3){

  # Load requireed packages
  require(ggvenn)
  require(dplyr)
  require(reshape2)
  require(magrittr)
  require(purrr)
  require(tidyr)

  # Check that inputs are OK
  if(!is.list(peak_list)){ stop("'peak_list' must be a (named) list of dataframes with the columns 'seqnames', 'start' and 'end'.") }
  else if(length(peak_list) != length(peak_names)){ stop("'peak_names' must be a character vector with the same length as 'peak_list'.") }

  # Get Venn Counts and the peaks in each set
  x <- getVennCounts(peak_list)

  # Transform the matrix of thee peaks in each set in ordeer to plot the Venn
  y <- x$matrix %>%
    as_tibble() %>%
    magrittr::set_colnames(c("Peak", paste("cond", 1:length(peak_list), sep = ""))) %>%
    reshape2::melt() %>% mutate(value = if_else(value == 1, Peak, NULL)) %>%
    tidyr::pivot_wider(names_from = "variable", values_from = "value") %>%
    dplyr::select(-Peak) %>%
    dplyr::as_tibble() %>%
    as.list() %>%
    purrr::set_names(peak_names)

  # Draw the Venn diagram with ggvenn
  venn <- ggvenn::ggvenn(data = y, show_percentage = percent,
                         fill_color = in_fill, fill_alpha = alpha,
                         stroke_color =  out_color,
                         set_name_color = name_color, set_name_size = name_size,
                         text_color = text_color, text_size = label_size)

  # Return the Venn diagram.
  return(venn)

}
