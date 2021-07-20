# DEcompare()
# ============================================= #

#' @title DEcompare
#' @author amitjavilaventura
#'
#' @description
#' Function that takes a list of 2 dataframes with expression data
#' and compares the log2FoldChange values of the two data sets.
#'
#' @param deg_list List of lenght 2. Must contain data frames with the columns Geneid, log2FoldChange, padj and DEG. Values in 'Geneid' must be the same in both data frames.+
#' @param threshold Numeric of lenght 1. Threshold for the log2FoldChange. It is used to draw horizontal and vertical lines. Default: 1.5.
#' @param genes Charachter of indefined lenght or NULL. Genes to be written in the plot. Default: NULL.
#' @param xlim Numeric of length 2. The limits of the X axis. Default: c(-10, 10)
#' @param xlim Numeric of length 2. The limits of the Y axis. Default: c(-10, 10)
#' @param xlab Character of lenght 1. Name of the X axis. Default: "Contrast1"
#' @param ylab Character of lenght 1. Name of the Y axis. Default: "Contrast2"
#' @param title Character of lenght 1. Title of the plot. Default: "Comparison of Log2FC"
#' @param subtitle Character of lenght 1. Subtitle of the plot. Default: paste(xlab, "vs", ylab)
#' @param colors Character of length 1 to 4. If too short, it will be recycled. Colors for each of the corners: top-left, bottom-left, top-right, bottom-right. Default: c("pink", "lightgreen", "cornflowerblue", "yellow")
#' @param alpha_corners Numeric of length 1 to 4. If too short, it will be recycled. Alpha (transparency) value for the vars. Default: 0.7
#'
#' @export

DEcompare <- function(deg_list, threshold = 1.5,
                      genes = NULL,
                      xlim = c(-10, 10), ylim = c(-10, 10),
                      xlab = "Contrast1", ylab = "Contrast2",
                      main = "Comparison of Log2FC",
                      subtitle = paste(xlab, "vs", ylab),
                      color_corners = c("pink", "lightgreen", "cornflowerblue", "yellow"),
                      alpha_corners = c(.7)){
  # Load packages
  require(dplyr)
  require(tibble)
  require(ggplot2)
  require(ggrepel)
  require(ggpubr)

  # Check that inputs are OK
  if(!is.list(deg_list) | length(deg_list) != 2){ stop("'deg_list' must be a list of 2 data frames with the columns 'Geneid', 'log2FoldChange', 'padj' and 'DEG'.") }
  else if(length(xlim) != 2 | length(ylim) != 2 ){ stop("Both 'xlim' and 'ylim' must be numeric vectors of lenght 2.") }


  data <- inner_join(x = deg_list[[1]] %>% select(Geneid, log2FoldChange, padj, DEG),
                     y = deg_list[[2]] %>% select(Geneid, log2FoldChange, padj, DEG),
                     by="Geneid")

  data <- data %>%
    mutate(color = if_else(log2FoldChange.x >= threshold & log2FoldChange.y >= threshold, "positivecorr",
                           if_else(log2FoldChange.x >= threshold & log2FoldChange.y <= -threshold, "negativecorr",
                                   if_else(log2FoldChange.x <= -threshold & log2FoldChange.y <= -threshold, "positivecorr",
                                           if_else(log2FoldChange.x <= -threshold & log2FoldChange.y >= threshold, "negativecorr", "ns")))))

  g <- ggplot(data, aes(log2FoldChange.x, log2FoldChange.y, color = color))

  if(!is.null(color_corners)){
    g <- g +
      annotate(geom = "rect",
               xmin = c(xlim[1], xlim[1], 0, 0),
               xmax = c(0, 0, xlim[2], xlim[2]),
               ymin = c(0,  ylim[1],  0,  ylim[1]),
               ymax = c(ylim[2], 0, ylim[2], 0),
               fill = color_corners, alpha = alpha_corners)
  }

  g <- g  +
    geom_point(alpha = 1) +

    coord_cartesian(xlim = xlim, ylim = ylim, expand = F) +
    geom_hline(yintercept = threshold, linetype = 2, color = "black") +
    geom_hline(yintercept = -threshold, linetype = 2, color = "black") +
    geom_vline(xintercept = threshold, linetype = 2, color = "black") +
    geom_vline(xintercept = -threshold, linetype = 2, color = "black") +
    scale_colour_manual(values=c("positivecorr" = "darkred", "negativecorr" = "darkgreen", "ns" = "gray50"),
                        drop = T) +

    ggtitle(main, subtitle) +
    xlab(xlab) + ylab(ylab) +

    theme_pubr(border = T, legend = "none") +
    theme(legend.title = element_blank(),
          plot.title = element_text(face="bold"),
          plot.subtitle = element_text(face="italic"),
          axis.title = element_text(face="bold"))

  if(!is.null(genes)){
    g <- g +
      geom_text_repel(data = data %>% filter(Geneid %in% genes),
                      mapping = aes(label = Geneid, x = log2FoldChange.x, y = log2FoldChange.y),
                      color = "black", size = 4)
  }


  # return
  return(g)
}

