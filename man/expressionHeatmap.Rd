% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expressionHeatmap.R
\name{expressionHeatmap}
\alias{expressionHeatmap}
\title{expressionHeatmap}
\usage{
expressionHeatmap(expr_df, genes, clust_rows = T, clust_cols = F, show_dend_rows = F, show_dend_cols = F, dist_method = "euclidean", hclust_method = "ward.D", write_label = T, label_size = 4, label_color = "black", label_digits = 2, hm_height = length(genes)*10, hm_width = (ncol(df)-1)*10, hm_colors = c("cornflowerblue", "white", "gold3"), legend_scale = NULL, legend_breaks_num = 5, legend_midpoint = 0, legend_height = hm_height, legend_title = NULL, dend_cols_prop = .1, dend_rows_prop = .2, title = "", subtitle = "", caption = NULL, xlab = "", ylab = NULL, axis_text_size = 10, x_axis_angle = 90)
}
\arguments{
\item{expr_df}{Dataframe with a 'Geneid' column and several columns with numerical expression data for different samples, such as TPMs or Log2FC.}

\item{genes}{Character. Names of the genes to be plotted. They must present in the column 'Geneid' of 'df'.}

\item{clust_rows}{Logical of length 1. Whether to cluster the rows (TRUE) or not (FALSE) using \code{hclust}. Default: T.}

\item{clust_cols}{Logical of length 1. Whether to cluster the columns (TRUE) or not (FALSE) using \code{hclust}. Default: F.}

\item{show_dend_rows}{Logical of length 1. Whether to draw the dendogram of the rows clustering (TRUE) or not (FALSE). It only works with 'clust_rows = T'. Default: F.}

\item{show_dend_cols}{Logical of length 1. Whether to draw the dendogram of the columns clustering (TRUE) or not (FALSE). It only works with 'clust_cols = T'. Default: F.}

\item{dist_method}{Character of length 1. Distance calculation. One of \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, "\verb{canberra"}, \code{"binary"} or \code{"minkowski"}. Default: "euclidean".}

\item{hclust_method}{Character of length 1. Method for hierarchical clustering. One of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"} (= UPGMA), \code{"mcquitty"} (= WPGMA), \code{"median"} (= WPGMC) or \code{"centroid"} (= UPGMC). Default: "ward.D"}

\item{write_label}{Logical of length 1. Whether to write the expression values in the heatmap (TRUE) or not (FALSE). Default: T.}

\item{label_size}{Numerical of length 1. Size of the expression values written in each cell of the heatmap. Default: 3.}

\item{label_color}{Character of length 1. Color of the expression values written in each cell of the heatmap. Default: "black".}

\item{label_digits}{Numerical of length 1. Number of digits the expression values are rounded to. Defalut: 2.}

\item{hm_height}{Numerical of length 1. Height of the heatmap in mm. Default: length(genes)*10.}

\item{hm_width}{Numerical of length 1. Width of the heatmap in mm. Default: (ncol(df)-1)*10.}

\item{hm_colors}{Character of length 3. Colors in the lower limit, midpoint (defined by 'legend_midpoint') and higher limit, respectively. Default: c("cornflowerblue", "white", "gold3").}

\item{na_color}{Charachter of length 1. Color for the NA values. Default: "gray".}

\item{border_color}{Character of length 1 or NA. Color of the border of each cell in the heatmap. Default: NA.}

\item{legend_scale}{Numerical of length 2 or NULL. If NULL, the color scale of the heatmap will take the minimum and the maximum values as limits. If numerical, the color scale will take the first as the lower limit and the second element as the higher limit. Default: NULL.}

\item{legend_breaks_num}{Numerical of length 1. Only 'legend_scale' is NULL. The number of breaks you want in the legend. Default: 5.}

\item{legend_breaks_by}{Numerical of length 1. Only 'legend_scale' is numerical (e.g. \code{c(-1,1)}). The distance between the breaks of the legend. Default: .5.}

\item{legend_midpoint}{Numerical of length 1. Only if 'legend_scale' is not NULL. Point where the central color of the legend will placed. Default: 0}

\item{legend_height}{Numerical of length 1. Height of the legend which, by default, is the height of the heatmap. Default: hm_height.}

\item{legend_title}{Character of length 1 or NULL. Title of the legend, placed in the right part of it. Default: NULL.}

\item{dend_cols_prop}{Numerical of length 1. Proportion (from 0 to 1) of the columns dendogram compared to the heatmap height. Default: 0.1.}

\item{dend_rows_prop}{Numerical of length 1. Proportion (from 0 to 1) of the rows dendogram compared to the heatmap width Default: 0.2.}

\item{title}{Character of length 1 or NULL. Title of the plot. Default: "".}

\item{subtitle}{Character of length 1 or NULL. Subtitle of the plot. Default: "".}

\item{caption}{Character of length 1 or NULL. Caption of the plot; placed at the bottom. Default: NULL.}

\item{xlab}{Character of length 1 or NULL. Title of the X axis. Default: "}

\item{ylab}{Character of length 1 or NULL. Title of the Y axis. If row dendogram is plotted, the Y axis is placed to the right. Default: NULL.}

\item{title_size}{Numerical of length 1. Size of the plot title. Default: 13.}

\item{subtitle_size}{Numerical of length 1. Size of the plot subtitle. Default: 12.}

\item{caption_size}{Numerical of length 1. Size of the plot caption. Default: 7.}

\item{axis_text_size}{Numerical of length 1. Size of the text in the axes. Default: 10.}

\item{axis_title_size}{Numerical of length 1. Size of the axes titles. Default: 11.}

\item{title_hjust}{Numerical of length 1. Horizontal justification for the title and the subtitle. Default: 0.}

\item{scale}{Character of length 1 or NULL. One of c("rows", "cols"). If "rows", it scales data by row; if "cols", it scales data by columns; if NULL (the default), it does not scale the data. Default: NULL.}
}
\description{
Function that takes a data frame with a Geneid column and several columns for expression data (e.g. TPMs, Log2FC...)
and draws a heatmap using ggplot.
}
\seealso{
\code{plotDendogram}

\code{stats::hclust}

\code{stats::dist}
}
\author{
amitjavilaventura
}
