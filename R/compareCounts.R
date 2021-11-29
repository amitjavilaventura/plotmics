

### compareCounts()
### ===============================================================================================

#' @title compareCounts()
#' @author amitjavilaventura
#'
#' @description
#' Draws a scatter plot of gene counts from two samples.
#'
#' @param df Dataframe with at least three columns: Geneid and one column for the counts of each sample.
#' @param cols Character of length 2. The names of the columns where the counts to be plotted are. Default: colnames(df)[2:3].
#' @param lims Numerical of length 2. Lower and upper limits of the axes. The same limits are used for both axes. Default: c(0,max(df[,cols[1]]*0.75)).
#' @param diag Logical of length 1. Whether to draw a dashed line in the diagonal or not. Default: TRUE.
#' @param writeN. Logical of length 1. Whether to write the total number of genes in the top left corner of the plot. Default: TRUE.
#' @param highlightGenes Character of undefined length or NULL. The names of the genes (values of Geneid column) to highlight. Default: NULL.
#' @param writeGeneLabels Logical of length 1. Whether to write the labels of the highlighted genes. Only if highlightGenes is not null. Default: TRUE.
#' @param geneLabelSize Numerical of length 1. Size of the gene label. Default: 3.
#' @param title Character of length 1. Title of the plot. Default: NULL.
#' @param subtitle Character of length 1. Subtitle of the plot. Default: NULL
#'
#'
#' @export

compareCounts <- function(df,
                          cols = colnames(df)[2:3],
                          lims = c(0,max(df[,cols[1]]*0.75)),
                          diag = TRUE,
                          title = NULL,
                          subtitle = NULL,
                          writeN = TRUE,
                          highlightGenes = NULL,
                          writeGeneLabels = TRUE,
                          geneLabelSize = 3){

  # Load required packages
  require(dplyr, quietly = T)
  require(magrittr, quietly = T)
  require(ggplot2, quietly = T)
  require(ggpubr, quietly = T)

  # Format dataframe
  d <- df %>%
    dplyr::select(c("Geneid", cols)) %>%
    dplyr::mutate(shape = "circle",
                  shape = ifelse(!!sym(cols[1]) > lims[2], "triangle", shape),
                  shape = ifelse(!!sym(cols[2]) > lims[2], "triangle", shape),
                  x = ifelse(!!sym(cols[1]) > lims[2], lims[2], !!sym(cols[1])),
                  y = ifelse(!!sym(cols[2]) > lims[2], lims[2], !!sym(cols[2])))

   # Initialize ggplot
   g <- ggplot(data = d, aes(x = x , y = y, color = Geneid))

   # Draw points and highlight them or not and draw labels if necessary
   if(!is.null(highlightGenes)){
     # Load gghighlight
     require(gghighlight)
     # Draw and highlight points
     g <- g + geom_point(aes(shape=shape), show.legend = F, alpha = .9) +
       gghighlight::gghighlight(Geneid %in% highlightGenes, label_key = Geneid, use_direct_label = writeGeneLabels,
                                label_params = list(size = geneLabelSize, box.padding = 0.15, label.padding = 0.15,
                                                    force = 1, force_pull = 0.5, fill = NA, label.size = 0))
   } else {
     g <- g + geom_point(aes(shape=shape), show.legend = F, alpha = .9, color = "black")
   }

   # Draw diagonal. Default TRUE
   if(diag){ g <- g + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black") }

   # Write total number of genes. Default: TRUE
   if(writeN){ g + annotate(geom = "text", label = paste("N=",nrow(d),sep=""), x=0, y=lims[2], hjust = 0, fontface = "italic") }

   # Format
   g <- g +
     # Set shape and alpha
     scale_shape_manual(values = c("circle", "triangle")) +
     scale_alpha_manual(values = c(.9, .7)) +
     # Set limits
     coord_cartesian(ylim = lims, xlim = lims, clip = "off") +
     # Draw title, labels...
     labs(x = cols[1], y = cols[2], title = title, subtitle = subtitle) +
     # Set theme: draw border and set panel grid major
     ggpubr::theme_pubr(legend = "none", border = F) +
     theme(panel.grid.major = element_line(colour = "gray70", linetype = 2, size = .4))



   return(g)

}
