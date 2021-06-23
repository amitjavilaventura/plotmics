# getVennCounts() -------------------------

#' @title getVennCounts
#' @author amitjavilaventura
#'
#' Function that calls 'ChIPpeakAnno::makeVennDiagram()' to obtain the venn counts.
#' It also does a list of peaks and with several variables corresponding to the presence of these peaks in each set.
#' It's a helper for other functions to draw Venn diagrams or UpSet plots for peak intersections.
#'
#' @seealso ChIPpeakAnno
#' @seealso plyranges
#'
#' @param peaks List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param conds Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param conds_order Character of the same length as 'conds'. Same values as in 'conds' but with the desired order of priority. Default: 'conds'.
#' @param plot Logical of lenght 1. Whether to draw the the 'makeVennDiagram()' plot or not. Default: FALSE
#'
#' @export


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

  peaks <- peaks %>% purrr::set_names(nm = conds)

  overlaps <- peaks %>%
    purrr::map(~plyranges::as_granges(.x)) %>%
    makeVennDiagram(Peaks = ., plot = plot) %>%
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
