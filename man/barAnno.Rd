% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barAnno.R
\name{barAnno}
\alias{barAnno}
\title{barAnno}
\usage{
barAnno(
  anno_list,
  anno_names = names(anno_list),
  names_order = unique(anno_names),
  protein = NULL,
  protein_order = unique(protein),
  main = NULL,
  subtitle = NULL,
  ylab = NULL,
  xlab = NULL,
  color_palette = "Set2",
  legend_position = "right",
  anno_num = 2,
  fill_position = T,
  xangle = 20,
  width = 0.6,
  counts_label = F,
  counts_angle = 90,
  counts_size = 3
)
}
\arguments{
\item{anno_list}{Named list of annotation objects that come from annotatePeak().}

\item{anno_names}{Charachter vector of the same length as 'anno_list'. Names that will be given to each of the objects in anno_list. Not that will be the names plotted in the bargraph}

\item{names_order}{Character vector with the same entries as 'names' with the order wanted to plot the data. Defaults to: unique(names).}

\item{protein}{NULL or haracter vector of the same length as 'anno_list' with the protein chipped in each dataframe of anno_list. If not NULL, these names will be passed through \code{facet_wrap()}. Default: NULL.}

\item{protein_order}{Character vector with the same entries as 'protein' with the order wanted to plot the data. Defaults to: unique(protein).}

\item{main}{Character of lenght 1. Title of the plot. Default: NULL.}

\item{subtitle}{Character of lenght 1. Subtitle of the plot. Default: NULL.}

\item{ylab}{Character of lenght 1. Title of the Y axis. Default: NULL.}

\item{xlab}{Character of lenght 1. Title of the X axis. Default: NULL}

\item{legend_position}{Character of lenght 1. Position of the legend. One of c("none", "bottom", "right", "left," "top"). Default: "right"}

\item{anno_num}{Numerical or character of length 1. Number of annotations to plot, either 2 (Promoter/Distal), 3 (Promoter/Gene body/Distal) or 'all' (the annotatePeak() default). Default: 2.}

\item{fill_position}{Logical. If TRUE (default), it the plotted bars will represent proportion of peaks in each feature. If FALSE, the bars will have the height of the total number of peaks with the correspondent feature color.}

\item{xangle}{Numerical of length 1. Angle of the text in the X axis. Default: 20.}

\item{width}{Numerical of length 1. Width of the bar in relative units. Default: 0.6}

\item{counts_label}{Logical of length 1. Whether to write the number and percentage of observations in each annotation. Default: FALSE.}

\item{counts_angle}{Numerical of length 1. Angle of the number and percentage labels. Default: 90}

\item{counts_size}{Numerical of length 1. Size of the number and percentage labels. Default: 3}

\item{palette}{Character of lenght 1. Color palette used to color the bars through the function \code{scale_fill_brewer()}. Default: "Set1".}
}
\description{
Function for ChIP-seq and ATAC-seq.
It must be used after the function annotatePeak() from the R package ChIPseeker. @seealso \code{\link{`annotatePeak`}}
It takes a list of annotation objects that come as output of annotatePeak() and changes the features to "Promoter", "Distal" and "Gene body" (or to "Promoter" and "Distal"). Finally it plots a bargraph with the distribution of all the proportions
As a ggplot2-based function, it allows to add more layers to format the plot.
}
\seealso{
\code{ChIPseeker::annotatePeak}
}
\author{
amitjavilaventura
}
