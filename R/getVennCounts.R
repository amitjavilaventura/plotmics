# getVennCounts() -------------------------

#' @title getVennCounts
#' @author amitjavilaventura
#'
#' @description
#' Function that calls 'ChIPpeakAnno::makeVennDiagram()' to obtain the venn counts.
#' It also does a list of peaks and with several variables corresponding to the presence of these peaks in each set.
#' It's a helper for other functions to draw Venn diagrams or UpSet plots for peak intersections.
#'
#' @seealso `ChIPpeakAnno::makeVennDiagram`
#'
#' @param peaks List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param conds Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param conds_order Character of the same length as 'conds'. Same values as in 'conds' but with the desired order of priority. Default: 'conds'.
#' @param stranded Logical of length 1. If TRUE, it considers the strand information to do the overlaps. Default = FALSE.
#' @param plot Logical of lenght 1. Whether to draw the the 'makeVennDiagram()' plot or not. Default: FALSE
#'
#'
#' @export


getVennCounts <- function(peaks,
                          conds = names(peaks),
                          conds_order = conds,
                          stranded = FALSE,
                          plot = FALSE){

  # Load required packages
  (require(dplyr, quietly = T))
  (require(purrr, quietly = T))
  (require(plyranges, quietly = T))
  (require(magrittr, quietly = T))
  (require(ChIPpeakAnno, quietly = T))

  if(!is.list(peaks)){ stop("'peaks' must be a (named) list of dataframes with the columns 'seqnames', 'start' and 'end'.") }
  else if(is.null(conds)){ stop("'conds' must a not-NULL character vector with the conditions of the data frames in 'peaks.") }
  else if(length(peaks) != length(conds)){ stop("'peaks' and 'conds' must have the same length.") }

  len <- length(peaks)

  peaks <- peaks %>% purrr::set_names(nm = conds)

  if(!stranded){

    overlaps <- suppressMessages(peaks %>%
                                   purrr::map(~plyranges::as_granges(.x)) %>%
                                   makeVennDiagram(Peaks = ., NameOfPeaks = conds, plot = plot, ignore.strand = !stranded))

    overlaps <- overlaps$vennCounts

  } else {

    # Separate peaks by strand
    peaks_plus     <- peaks %>% purrr::map(~dplyr::filter(.x, strand == "+"))
    peaks_minus    <- peaks %>% purrr::map(~dplyr::filter(.x, strand == "-"))
    peaks_nostrand <- peaks %>% purrr::map(~dplyr::filter(.x, strand %in% c(".", "*"))) %>% purrr::map(~dplyr::mutate(.x, strand = "*"))

    # Overlap peaks by strand
    overlaps_plus     <- suppressMessages(peaks_plus %>% purrr::map(~plyranges::as_granges(.x)) %>%
                                        makeVennDiagram(Peaks = ., NameOfPeaks = conds, plot = plot))

    overlaps_minus    <- suppressMessages(peaks_minus %>% purrr::map(~plyranges::as_granges(.x)) %>%
                                        makeVennDiagram(Peaks = ., NameOfPeaks = conds, plot = plot))

    overlaps_nostrand <- suppressMessages(peaks_nostrand %>% purrr::map(~plyranges::as_granges(.x)) %>%
                                        makeVennDiagram(Peaks = ., NameOfPeaks = conds, plot = plot))

    overlaps_plus     <- overlaps_plus$vennCounts
    overlaps_minus    <- overlaps_minus$vennCounts
    overlaps_nostrand <- overlaps_nostrand$vennCounts

    overlaps <- overlaps_plus
    for (i in 1:nrow(overlaps)){

      overlaps[i,ncol(overlaps)] <- sum(overlaps_plus[i,ncol(overlaps)], overlaps_minus[i,ncol(overlaps)], overlaps_nostrand[i,ncol(overlaps)])

    }
  }

  # Create matrix of overlapping peaks
  matrix <- matrix(data = rep(0, len), ncol = len, byrow=T) %>% set_colnames(conds)
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
