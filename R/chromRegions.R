### chromRegions()
### ===============================================================================================

#' @title chromRegions
#' @author amitjavilaventura
#'
#' @description
#' Function that draws the chromosomes and the position of the desired regions.
#'
#' @usage chromRegions(chrom_sizes, regions, regions2 = NULL, chr_exclude = c("Un", "Random", "JH", "GL", "\\."), chr_order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT"), strand_color = c("Darkred", "Blue4"), regions_names = c("Regions1", "Regions2"), regions_colors = c("Black", "Red"), title = NULL, subtitle = NULL, caption = NULL, xlab = "Position (bp)", ylab = "Chromosome", legend = "right", draw_points = TRUE, coord_flip = TRUE)
#'
#' @param chrom_sizes Character of lenght 1 or dataframe. Path to the chrom.sizes file or data frame with the names of the chromosomes and their sizes. The file must not contain column names.
#' @param regions Character of length 1 or dataframe. Path to BED/TSV file or data frame with chromosome names, start, end, region id, length and strand. Strand must be "+", "-" or ".".
#' @param regions2 Character of length 1, dataframe or NULL. Path to BED/TSV file or data frame with chromosome names, start, end, region id, length and strand. Strand must be "+", "-" or ".". Default: NULL.
#' @param chr_exclude Character. Regular expressions to match the chromosomes names to exclude . Default: c("Un", "Random", "JH", "GL", "\\.").
#' @param chr_order Character. Chromosome names written in the order to be plotted. The names must matchthe chromosome names in chrom_sizes and regions. Default: c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT")
#' @param strand_color Character of lenght 2 or NULL. Colors for each strand of "regions". If regions2 is not NULL or there is a "." in the regions strand, this will be converted NULL. Default: c("Darkred", "Blue4").
#' @param regions_names Character of length 2. Names of the input regions. If regions2 is NULL, only the first value will be used. Default: c("Regions1", "Regions2").
#' @param regions_colors Character of length 2. Colors of the input regions. Only if strand_color is NULL. If regions2 is NULL, only the first value will be used. Default: c("Black", "Red").
#' @param title Character of length 1 or NULL. Title of the plot. Default: NULL
#' @param subtitle Charachter of length 1 or NULL. Subtitle of the plot. Default: NULL.
#' @param caption Character of length 1, TRUE or NULL. Caption to be written in the bottom-right corner. If TRUE, it will write the number of regions in the input; if charachter, it will write the written caption. Default: NULL.
#' @param ylab Character of length 1. Title of the X axis. Default: "Position (bp)"
#' @param xlab Character of lenght 1. Title of the Y axis. Default: "Chromosome".
#' @param legend Character of length 1. Position of the legend, passed through ggpubr::theme_pubr(). One of c("top", "bottom", "left", "right", "none"). Default: "right".
#' @param draw_points Logical of lenght 1. If TRUE, draws points in each of the drawn regions. Default: TRUE
#' @param coord_flip Logical of lenght 1. If TRUE, change the position of the axes (Y axis is position and X is chromosome). Default: FALSE
#'
#' @export


chromRegions <- function(chrom_sizes,
                         regions,
                         regions2       = NULL,
                         chr_exclude    = c("Un", "Random", "JH", "GL", "\\."),
                         chr_order      = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT"),
                         strand_color   = c("Darkred", "Blue4"),
                         regions_names  = c("Regions1", "Regions2"),
                         regions_colors = c("Black", "Red"),
                         title          = NULL,
                         subtitle       = NULL,
                         caption        = NULL,
                         ylab           = "Position (bp)",
                         xlab           = "Chromosome",
                         legend         = "bottom",
                         draw_points    = TRUE,
                         coord_flip     = FALSE) {

  # Load packages -----
  require(dplyr)
  require(magrittr)
  require(ggplot2)
  require(ggpubr)

  # Check that inputs are OK -----
  ## Chrom sizes
  if(!class(chrom_sizes) %in% c("character", "data.frame", "tibble")) { stop("'chrom_sizes' must be a character vector with the path to a chrom.sizes file or a data.frame with the chromosomes and their sizes") }
  else if(class(chrom_sizes) == "character") { chromsizes <- read.delim(chrom_sizes, header = F, col.names = c("seqnames", "end")) }
  else if(class(chrom_sizes) %in% c("data.frame", "tibble")) { chromsizes <- chrom_sizes %>% magrittr::set_colnames(c("seqnames", "end")) }

  ## Regions
  if(!class(regions) %in% c("character", "data.frame", "tibble")) { stop("'regions' must be a character vector with the path to a TSV/BED file or a data.frame with the genomic coordinates, the id, the length and the strand for each region") }
  else if(class(regions) == "character") { regions <- read.delim(regions, header = F, col.names = c("seqnames", "start", "end", "id", "length", "strand"))}
  else if(class(regions) %in% c("data.frame", "tibble")) { regions <- regions %>% magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand")) }

  ## Seqnames in regions = seqnames in chrom sizes
  if(!all(regions$seqnames %in% chromsizes$seqnames)) { stop("The chromosome names in 'regions' must be equal to the chromosome names in 'chrom_sizes'") }

  ## No strandness
  if(any(regions$strand == ".")) { strand_color = NULL }


  ## Regions2 stuff
  if(!is.null(regions2)){
    if(!class(regions2) %in% c("character", "data.frame", "tibble")) { stop("'regions2' must be a character vector with the path to a TSV/BED file or a data.frame with the genomic coordinates, the id, the length and the strand for each region") }
    else if(class(regions2) == "character") { regions2 <- read.delim(regions2, header = F, col.names = c("seqnames", "start", "end", "id", "length", "strand"))}
    else if(class(regions2) %in% c("data.frame", "tibble")) { regions2 <- regions2 %>% magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand")) }

    ## Seqnames in regions = seqnames in chrom sizes
    if(!all(regions2$seqnames %in% chromsizes$seqnames)) { stop("The chromosome names in 'regions2' must be equal to the chromosome names in 'chrom_sizes'") }
  }

  ## Chromosomes to exclude
  if(!class(chr_exclude) == "character") { stop("'chr_exclude' must be a character vector with regular expressions that match the chromosomes to exclude") }

  ## Order of the chromosomes
  if(!class(chr_order) == "character") { stop("'chr_order' must be a character vector with the names of the chromosomes in the desired order") }

  # Filter chromsizes and regions ------
  chrom_filt <- chromsizes %>% dplyr::filter(!stringr::str_detect(seqnames, pattern = paste(chr_exclude, collapse = "|"))) %>% dplyr::mutate(seqnames = factor(seqnames, levels = chr_order))
  regions_filt <- regions %>% dplyr::filter(!stringr::str_detect(seqnames, pattern = paste(chr_exclude, collapse = "|"))) %>% dplyr::mutate(seqnames = factor(seqnames, levels = chr_order))

  # Make plot -----
  # Initialize chromosomes
  g <- ggplot() + geom_col(data = chrom_filt, mapping = aes(y = end, x = seqnames), color = "black", fill = NA, width = 0.7)

  # Draw regions
  ## If there are 2 regions to plot
  if(!is.null(regions2)){

    # Nullify strand color
    strand_color <- NULL

    # Filter regions2
    regions2_filt <- regions2 %>% dplyr::filter(!stringr::str_detect(seqnames, pattern = paste(chr_exclude, collapse = "|"))) %>% dplyr::mutate(seqnames = factor(seqnames, levels = chr_order))

    # Merge regions1 and regions2
    regions_merged <- regions_filt %>% dplyr::mutate(region = regions_names[1]) %>%
      bind_rows(regions2_filt %>% dplyr::mutate(region = regions_names[2]))

    # Draw regions 1 and 2
    g <- g + geom_tile(data = regions_merged, mapping = aes(x = seqnames, y = start + length/2, height = length, width = .7, fill = region, color = region))
    # Draw points if desired
    if(draw_points) { g <- g + geom_point(data = regions_merged, mapping = aes(x = seqnames, y = start + length/2, fill = region, color = region), size = 1) }
    # Change colors
    g <- g + scale_fill_discrete(type = regions_colors) + scale_color_discrete(type = regions_colors)

    # If there is only 1 set of regions and the strand_color is not null
  } else if(is.null(regions2) & !is.null(strand_color)){

    # Draw regions
    g <- g + geom_tile(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, height = length, width = .7, fill = strand, color = strand))
    # Draw points if desired
    if(draw_points) { g <- g + geom_point(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, fill = strand, color = strand), size = 1) }
    # Change colors
    g <- g + scale_fill_discrete(type = strand_color) + scale_color_discrete(type = strand_color)

    # If there is only 1 set of regions and the strand_color is null
  } else if(is.null(regions2) & is.null(strand_color)){

    # Draw regions
    g <- g +  geom_tile(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, height = length, width = .7), fill = regions_colors[1], color = regions_colors[1])
    # Draw points if desired
    if(draw_points) { g <- g + geom_point(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2), fill = regions_colors[1], color = regions_colors[1], size = 1) }

  }

  # Write axis labels
  g <- g + xlab(xlab) + ylab(ylab)

  # Write title and subtitle
  if(!is.null(title)){  g <- g + ggtitle(label = title, subtitle = subtitle) }

  # Write caption
  if(!is.null(caption)){

    # If caption is TRUE, write number of regions (only regions1 or regions 1 and 2)
    if(is.logical(caption)) {
      if(caption) {
        if(is.null(regions2)) { g <- g + labs(caption = paste("Number of regions:", nrow(regions_filt))) }
        else { g <- g + labs(caption = paste(regions_names[1],": ", nrow(regions_filt), "; ", regions_names[2], ": ", nrow(regions2_filt), sep = "")) }
      }
    }

    # If caption is character, write the caption.
    if(is.character(caption)) { g <- g + labs(caption = caption) }

  }

  # Customize plot
  g <- g +
    ## Remove space between axis and plot
    scale_y_continuous(expand = c(0,0)) +
    ## Theme pubr to remove border and choose the legend position
    theme_pubr(border = F, margin = T, legend = legend) +
    # Further costumization
    theme(plot.title = element_text(hjust = .5, face = "bold"), # title
          plot.subtitle = element_text(hjust = .5, face = "italic"), # subtitle
          axis.title = element_text(hjust = .5, face = "bold"), # axis titles
          legend.title = element_blank(), # remove legend title
          legend.key.size = unit(4, "mm"), #
          axis.ticks = element_blank(), axis.line = element_blank())

  # Flip coordinates
  if(coord_flip) { g <- g + coord_flip() }

  # Return plot
  return(g)
}
