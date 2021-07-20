# plotDendogram() -------------------------

#' @title plotDendogram
#' @author amitjavilaventura
#'
#' @description
#' Helper function that takes clustering data (hclust()) as input and
#' plots a dendogram using ggplot and ggdendro
#'
#' @seealso `ggpubr`
#' @seealso `ggdendro`
#'
#' @param hclust_data Output from `hclust()` function.
#' @param axis Characteer of length 1. Either "rows" or "cols". If rows, the dendogram is made to be placed in the left of the Y axis; if cols, the dendogram is made to be placed in the top of the X axis.
#'
#' @export

plotDendogram <- function(hclust_data, axis = "rows"){

  # Load required packages
  require(dplyr)
  require(ggdendro)
  require(ggplot2)
  require(ggpubr)

  dend_init    <- hclust_data %>% as.dendrogram()
  dend_data    <- dendro_data(dend_init)
  segment_data <- with( segment(dend_data), data.frame(x = y, y = x, xend = yend, yend = xend))
  pos_table    <- with( dend_data$labels, data.frame(y_center = x, label = as.character(label), height = 1))
  dend   <- ggplot(segment_data) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_y_discrete(breaks = pos_table$y_center,
                     labels = pos_table$label,
                     expand = c(0, 0.5)) +
    theme_pubr(border = F) +
    theme(text = element_blank(), line = element_blank())

  if(axis == "rows"){ dend <- dend + scale_x_reverse() }
  else if(axis == "cols"){ dend <- dend + coord_flip() }

  return(dend)

}
