### chromReads()
### ===============================================================================================

#' @title chromReads
#' @author amitjavilaventura
#'
#' @usage chromReads(bamfile, chr.filt = c("Un", "random") main = NULL, main.size = 13, subtitle = NULL, sub.size = 11, xlab = "Mapped reads", ylab = "Chromosome", axis.size = 9, percent.label = T, percent.size = 3)
#' Function counts the reads mapped to each chromosome and plots a bar graph.
#' It uses the function idxstatsBam() from the Rsamtools R package. Run help("idxstatsBam") for more information.
#'
#' @param bamfile Character of lenght 1. Path and name of the BAM file whose reads are to be mapped.
#' @param chr.filt Character. Vector of undefined length with strings to filter chromosomes (i.e. "Un" would filter all chromosomes containing "Un" in their name). Default: c("Un", "random").
#' @param main Character of lenght 1. Title of the pie chart. Default: NULL.
#' @param main.size Numeric of length 1. Font size of the title. It works only if main is not NULL. Default: 13.
#' @param subtitle Character of lenght 1. Subtitle of the bar plot. It works only if main is not NULL. Default: NULL.
#' @param sub.size Numeric of length 1. Font size of the subtitle. It works only if main is not NULL. Default: 11.
#' @param xlab Character of lenght 1. Title of the X axis. Default: "Mapped reads".
#' @param ylab Character of lenght 1. Title of the Y axis. Default: "Chromosome".
#' @param axis.size Numeric of length 1. Size of the axis labels. Default: 9
#' @param percent.label Logical of lenght 1. If TRUE, it draws the percentage of reads in each chromosome.
#' @param percent.size Character of lenght 1. Size of the percentage drawn in the plot.
#'
#' @export

chromReads <- function(bamfile, chr.filt = c("Un", "random"),
                       main = NULL, main.size = 13, subtitle = NULL, sub.size = 11,
                       xlab = "Mapped reads", ylab = "Chromosome", axis.size = 9,
                       percent.label = T, percent.size = 3){

  # Load required packages
  require(Rsamtools)
  require(ggplot2)
  require(dplyr)
  require(stringr)

  # Calculate the number of reads mapping to each chromosome with Rsamtools::idxstatsBam()
  chromReads <- idxstatsBam(bamfile)

  # Remove names of strange chromosomes
  chromReads <- chromReads[grep("chr", chromReads$seqnames),]

  # Calculate the total number of mapped reads in "good" chromosomes
  totalReads <- sum(chromReads$mapped)

  # Calculate percentage of mapped reads in each chromosome against all mapped reads
  chromReads$percentage <- chromReads$mapped/totalReads*100

  # Filter strange chromosomees (i.e. "chrUn...")
  for(i in chr.filt){
    chromReads <- chromReads[(str_detect(chromReads$seqnames, pattern = i, negate = T)),]
  }


  # Draw a bar graph
  b <- ggplot(data = chromReads, mapping = aes(mapped, seqnames, fill = seqnames)) +
    geom_bar(stat = "identity", show.legend = legend, colour = "Gray15") +

    ylab(label = ylab) + xlab(label = xlab) +
    ggtitle(main, subtitle) +

    # General formatting
    theme_chromReads(main.size = main.size, sub.size = sub.size, axis.size = axis.size)

  # Annotate labels (percentages)
  if(percent.label == T){
    b <- b + geom_text(aes(label = round(percentage, 2), x = mapped+max(mapped)*0.05), size = percent.size, hjust = 0) +
      xlim(0, max(chromReads$mapped*1.20))

  }


  # Return bar graph
  return(b)
}

