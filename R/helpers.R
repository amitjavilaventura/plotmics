
# HELPER FUNCTIONS                              #
# ============================================= #


# theme_barDEGs()
# ============================================= #

#' @title theme_barDEGs
#' @author amitjavilaventura
#'
#' Helper for the 'barDEGs' function.
#'
#' @param xaxis Logical of length 1. Whether to draw the text of the X axis (number of DE) or not. Default: FALSE.
#' @param yaxis Logical of length 1. Whether to draw the text of the Y axis (contrast names) or not. Default: FALSE.

theme_barDEGs <- function(xaxis = F, yaxis = F){

  # Load required packages
  require(ggpubr)

  # Custom theme
  t <-
    theme_pubr(border = F, margin = T, legend = "bottom") +
    theme(axis.title = element_blank(),
          legend.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())

  if(!xaxis){ t <- t + theme(axis.text.x = element_blank()) }
  if(!yaxis){ t <- t + theme(axis.text.y = element_blank()) }

  t
}

# theme_chromReads()
# ============================================= #

#' @title theme_chromReads
#' @author amitjavilaventura
#'
#' Helper for the 'chromReads' function.
#'
#' @param xaxis Logical of length 1. Whether to draw the text of the X axis (number of DE) or not. Default: FALSE.
#' @param yaxis Logical of length 1. Whether to draw the text of the Y axis (contrast names) or not. Default: FALSE.

theme_chromReads <- function(main.size = 11, sub.size = 10, axis.size = 9){

  # Load required packages
  require(ggpubr)

  # Custom theme
  t <-
    theme_pubr(border = T, margin = T, legend = "none") +
    theme(plot.title = element_text(face = "bold", size = main.size, hjust = .5),
          plot.subtitle = element_text(face = "italic", size = sub.size, hjust = .5),
          axis.title = element_text(face = "bold", size = axis.size),
          legend.title = element_blank())

  t
}
