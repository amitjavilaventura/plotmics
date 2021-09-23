### chromRegions()
### ===============================================================================================

#' @title chromRegions
#' @author amitjavilaventura
#'
#' @description
#' Function that draws the chromosomes and the position of the desired regions.
#'
#' @usage chromRegions(chrom_sizes, regions, chr_exclude = c("Un", "Random", "JH", "GL", "\\."), chr_order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT"), col_by_strand = FALSE, regions_names = names(regions_list), regions_colors = sample(x = colors(), size = length(regions_list), replace = F), title = NULL, subtitle = NULL, caption = NULL, xlab = "Position (bp)", ylab = "Chromosome", legend = "right", draw_points = TRUE, coord_flip = TRUE)
#'
#' @param chrom_sizes Character of lenght 1 or dataframe. Path to the chrom.sizes file or data frame with the names of the chromosomes and their sizes. The file must not contain column names.
#' @param regions_list (Named) List of character vectors or dataframes. List with paths to BED/TSV files or dataframes with chromosome names, start, end, region id, length and strand. Strand must be "+", "-" or ".".
#' @param chr_exclude Character. Regular expressions to match the chromosomes names to exclude . Default: c("Un", "Random", "JH", "GL", "\\.").
#' @param chr_order Character. Chromosome names written in the order to be plotted. The names must matchthe chromosome names in chrom_sizes and regions. Default: c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT")
#' @param col_by_strand Logical of length 1. If TRUE, it colors the drawn regions by their strandness. If strand fieald contains any "." (i.e. ChiP-seq data), it is internally turned to FALSE. Default = FALSE.
#' @param regions_names Character of length equal to regions_list. Names of the input regions. Default: names(regions_list).
#' @param regions_colors Character of length 2 or of length equal to regions_list. Colors of the input regions. If col_by_strand is FALSE, the length must be equal to the input regions, if col_by_strand is TRUE, the length must be 2. Default: sample(x = colors(), size = length(regions_list), replace = F).
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
                         regions_list,
                         regions_names   = names(regions_list),
                         colors          = c("Black", "Red", "Blue"),
                         col_by_strand   = FALSE,
                         chr_exclude     = c("Un", "Random", "JH", "GL", "\\."),
                         chr_order       = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT"),
                         title           = NULL,
                         subtitle        = NULL,
                         caption         = NULL,
                         ylab            = "Position (bp)",
                         xlab            = "Chromosome",
                         legend          = "bottom",
                         draw_points     = TRUE,
                         coord_flip      = FALSE) {

  # Load packages -----
  require(dplyr)
  require(magrittr)
  library(purrr)
  require(ggplot2)
  require(ggpubr)

  # Check that inputs are OK and read data -----
  ## Chrom sizes
  if(!class(chrom_sizes) %in% c("character", "data.frame", "tibble")) { stop("'chrom_sizes' must be a character vector with the path to a chrom.sizes file or a data.frame with the chromosomes and their sizes") }
  else if(class(chrom_sizes) == "character") { chromsizes <- read.delim(chrom_sizes, header = F, col.names = c("seqnames", "end")) }
  else if(class(chrom_sizes) %in% c("data.frame", "tibble")) { chromsizes <- chrom_sizes %>% magrittr::set_colnames(c("seqnames", "end")) }

  ## Length of regions names/colors are equal to length of regions list
  if(length(regions_list) != length(regions_names)){ stop("'regions_names' must be a character vector with the names of the sets in 'regions_list' and their length must be the same.")}
  #else if(length(regions_list) != length(regions_names)){ stop("'regions_colors' must be a character vector with the colors of the sets in 'regions_list' and their length must be the same.")}

  ## Regions
  if(!class(regions_list) == "list") { stop("'regions_list' must be a list with paths to regions files or data frames with the positions of the regions") }
  if(all(regions_list %>% purrr::map(~class(.x)) %>% unlist()=="data.frame")) {

    # Format regions
    regions <- regions_list %>%
      purrr::set_names(regions_names) %>%
      purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
      dplyr::bind_rows() %>%
      magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand", "region"))

  } else if(all(regions_list %>% purrr::map(~class(.x)) %>% unlist()=="character")) {

    # Read and format regions
    regions <- purrr::map(regions_list, read.delim, header = F) %>%
      purrr::set_names(regions_names) %>%
      purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
      dplyr::bind_rows() %>%
      magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand", "region"))

  } else if(any(regions_list %>% purrr::map(~class(.x)) %>% unlist()) %in% c("data.frame", "character")) { stop("'regions_list' must be a list with paths to regions files or data frames with the positions of the regions") }

  ## Seqnames in regions = seqnames in chrom sizes
  if(!all(regions$seqnames %in% chromsizes$seqnames)) { stop("The chromosome names in 'regions' must be equal to the chromosome names in 'chrom_sizes'") }

  ## No strandness --> turn color by strand to false
  if(any(regions$strand == ".")) {
    warning("Not all regions have a strand accepted value (i.e. '+' or '-'), converting 'col_by_strand' to FALSE.")
    col_by_strand = F
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
  if(!col_by_strand) {
    # Draw lines
    g <- g +  geom_tile(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, height = length, width = .7, fill = region, color = region))
    # Draw points if desired
    if(draw_points) { g <- g + geom_point(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, fill = region, color = region), size = 1) }
    # Change colors
    g <- g + scale_fill_discrete(type = colors) + scale_color_discrete(type = colors)
  } else {
    # Draw lines
    g <- g +  geom_tile(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, height = length, width = .7, fill = strand, color = strand))
    # Draw points if desired
    if(draw_points) { g <- g + geom_point(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, fill = strand, color = strand), size = 1) }
    # Change colors
    g <- g + scale_fill_discrete(type = colors) + scale_color_discrete(type = colors)
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
        regions_num <- regions_filt %>% dplyr::count(region) %>% dplyr::mutate(Num = paste(region, ": ", n, sep = "")) %>% dplyr::pull(Num)
        g <- g + labs(caption = paste(regions_num, collapse = "; "))
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
