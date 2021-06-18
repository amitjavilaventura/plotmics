makeVenn4upSet <- function(peaks, conds = names(peaks), conds_order = conds, plot = F){

  require(pkgcond)
  require(purrr) %>% suppress_messages() %>% suppress_warnings()
  require(dplyr) %>% suppress_messages() %>% suppress_warnings()
  require(plyranges) %>% suppress_messages() %>% suppress_warnings()
  require(magrittr) %>% suppress_messages() %>% suppress_warnings()
  require(ChIPpeakAnno) %>% suppress_messages() %>% suppress_warnings()

  if(is.null(conds)) { stop("'conds' must not be NULL")}

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

upsetPeaks <- function(peaks, conds = names(peaks), conds_order = conds,
                       title = NULL, order.by = "freq",
                       mainbar.y.label = "Intersect size", sets.x.label = "Set size",
                       ...){

  # Load required packages
  require(UpSetR)
  require(grid)
  require(gridExtra)

  #
  # call makeVenn4upSet to retrieve the peak matrix
  x <- makeVenn4upSet(peaks, conds, conds_order, plot = F)

  # draw upset plot with options
  if(is.null(title)) {
    upset <- UpSetR::upset(data = x$matrix, order.by = order.by, sets.x.label = sets.x.label, mainbar.y.label = mainbar.y.label)
    return(upset)

  } else {
    plot.new()
    UpSetR::upset(data = x$matrix, order.by = order.by, sets.x.label = sets.x.label, mainbar.y.label = mainbar.y.label)
    #grid.edit(gPath = 'arrange', name="upset")
    #vp <- grid.grab()
    #grid.arrange( grobs = list( vp ), top=title, cols=1 )
  }
}

#upsetPeaks(peaks)

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
#' @param title
#' @param subtitle
#' @param caption
#'
#' @export

ggUpsetPeaks <- function(peaks, conds = names(peaks), conds_order = conds,
                         order_by_freq = T,
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
  x <- makeVenn4upSet(peaks, conds, conds_order, plot = F)

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
    geom_col(fill = "black", width = 0.6) +
    #ggtitle("Title", "Subtitile") +
    ylab("Intersection size") +
    theme_pubclean(base_size = 9) +
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
    theme_pubclean(base_size = 9, flip = T) +
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
    theme_pubr(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          plot.title = element_text(face = "bold"))


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

#ggUpsetPeaks(peaks, title = "This is an upset plots", subtitle = "With two sets of peaks", caption = "the intersections are made with makeVennDiagram", order_by_freq = T)


