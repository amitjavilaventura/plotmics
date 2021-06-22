# upsetPeaks() -------------------------

#' @title upsetPeaks
#' @author amitjavilaventura
#'
#' This function calls 'getVennCounts()' and draws and UpSet plot using the UpSetR package.
#'
#' @seealso plyranges
#' @seealso UpSetR
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
