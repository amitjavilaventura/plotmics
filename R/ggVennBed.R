
# ggVennBed() -------------------------

#' @title ggVennBed
#' @author amitjavilaventura
#'
#' @description
#' Does the intersection of two BED files (A and B) and
#' draws the intersection with the unique peaks of A and B and
#' the intersections of A with B and B with A.
#'
#' @seealso `bedtoolsr::bt.intersect`
#' @seealso `ggvenn::ggvenn`
#'
#' @param a Data frame or character of length 1. Path to a BED file or data frame with BED format. If data frame, it will be written to a temporary file.
#' @param b Data frame or character of length 1. Path to a BED file or data frame with BED format. If data frame, it will be written to a temporary file.
#' @param setnames Character of length 2. Names of the data sets of 'a' and 'b'. Default: c("A", "B").
#' @param stranded Logical of length 1. If TRUE, the intersections are done regarding the strand. Default = TRUE.
#' @param color Character of length 2. Color of the Venn circles and the labels indicating the number of regions. Default: c("blue", "gold3").
#' @param namesize Numeric of length 1. Size of the names of each set. Default: 7.
#' @param labsize Numeric of length 1. Size of the labels indicating the number of unique and overlapping regions.
#'
#' @export

ggVennBed <- function(a,
                      b,
                      setnames = c("A", "B"),
                      stranded = T,
                      color    = c("blue", "gold3"),
                      namesize = 7,
                      labsize  = 5){

  # Load required packages
  suppressWarnings(require(dplyr, quietly = T, warn.conflicts = F))
  suppressWarnings(require(bedtoolsr, quietly = T, warn.conflicts = F))
  suppressWarnings(require(magrittr, quietly = T, warn.conflicts = F))
  suppressWarnings(require(ggplot2, quietly = T, warn.conflicts = F))
  suppressWarnings(require(ggvenn, quietly = T, warn.conflicts = F))
  suppressWarnings(require(colorspace, quietly = T, warn.conflicts = F))

  # Check that inputs are OK
  tmpfile_a <- "a.tmp.bed"
  tmpfile_b <- "b.tmp.bed"
  if(!is.character(a)) {
    if(is.data.frame(a)) {  write.table(a, tmpfile_a, quote = F, col.names = F, row.names = F, sep = "\t") }
    else { stop("'a' must be a path to a BED file or a data frame") }
    a <- tmpfile_a
  }
  if(!is.character(b)) {
    if(is.data.frame(b)) { write.table(b, tmpfile_b, quote = F, col.names = F, row.names = F, sep = "\t") }
    else { stop("'b' must be a path to a BED file or a data frame") }
    b <- tmpfile_b
  }

  # Count total number of regions in A and B
  aa <- read.delim(a, header = F) %>% nrow()
  bb <- read.delim(b, header = F) %>% nrow()

  # Count the regions of A overlapping with B, and vicecersa
  ab <- bedtoolsr::bt.intersect(a = a, b = b, s = stranded, u = T) %>% nrow()
  ba <- bedtoolsr::bt.intersect(a = b, b = a, s = stranded, u = T) %>% nrow()

  # Count the unique peaks in A and B
  aa <- aa-ab
  bb <- bb-ba

  # Create a pseudo data for the backbone of the Venn diagram
  data <- list(rep("A",2), rep("B",2)) %>% magrittr::set_names(setnames)
  # Draw the Venn diagram
  venn <- ggvenn::ggvenn(data, show_percentage = F)
  # Remove labels from Venn diagram
  venn$layers[[4]] <- NULL

  # Annotate the unique peaks and overlaps in the Venn diagram
  venn <- venn +
    scale_x_discrete(expand = c(0,0.5)) +
    scale_y_discrete(expand = c(0,0.5)) +
    annotate("text", -.8, 0, label = as.character(aa), color = colorspace::darken(color[1], amount = 0.5), fontface = "bold", size = labsize) +
    annotate("text", .8, 0, label = as.character(bb),  color = colorspace::darken(color[2], amount = 0.5), fontface = "bold", size = labsize) +
    annotate("text", 0, 0.15, label = as.character(ab), color = colorspace::darken(color[1], amount = 0.5), fontface = "plain", size = labsize) +
    annotate("text", 0, -0.15, label = as.character(ba), color = colorspace::darken(color[2], amount = 0.5), fontface = "plain", size = labsize)

  # Remove temporary files
  if(file.exists(tmpfile_a)) { file.remove(tmpfile_a) }
  if(file.exists(tmpfile_b)) { file.remove(tmpfile_b) }

  # Return the Venn diagram
  return(venn)

}
