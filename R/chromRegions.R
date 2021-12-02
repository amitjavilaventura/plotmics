### chromRegions()
### ===============================================================================================

#' @title chromRegions
#' @author amitjavilaventura
#'
#' @description
#' Function that draws the chromosomes and the position of the desired regions.
#'
#' @param chrom_sizes Character of lenght 1 or dataframe. Path to the chrom.sizes file or data frame with the names of the chromosomes and their sizes. The file must not contain column names.
#' @param regions_list (Named) List of character vectors or dataframes. List with paths to BED/TSV files or dataframes with chromosome names, start, end, region id, length and strand. Strand must be "+", "-" or ".".
#' @param chr_exclude Character. Regular expressions to match the chromosomes names to exclude . Default: c("Un", "Random", "JH", "GL", "\\.").
#' @param chr_order Character. Chromosome names written in the order to be plotted. The names must matchthe chromosome names in chrom_sizes and regions. Default: c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT")
#' @param color_by Character of length 1. One of 'region', 'strand' or 'extra'. If 'region', colours by the regions in 'regions_list'. If 'strand', colours by strand. If 'extra' colours by the grouping provided field in 'extra_info' (if provided).
#' @param regions_names Character of length equal to regions_list. Names of the input regions. Default: names(regions_list).
#' @param regions_order Character of length equal to regions_list. Names of the input regions in the desired order. Default: names(regions_list).
#' @param colors Character of length 2 or of length equal to regions_list. Colors of the input regions. If col_by_strand is FALSE, the length must be equal to the input regions, if col_by_strand is TRUE, the length must be 2. Default: sample(x = colors(), size = length(regions_list), replace = F).
#' @param title Character of length 1 or NULL. Title of the plot. Default: NULL
#' @param subtitle Charachter of length 1 or NULL. Subtitle of the plot. Default: NULL.
#' @param caption Character of length 1, TRUE or NULL. Caption to be written in the bottom-right corner. If TRUE, it will write the number of regions in the input; if charachter, it will write the written caption. Default: NULL.
#' @param ylab Character of length 1. Title of the X axis. Default: ""
#' @param xlab Character of lenght 1. Title of the Y axis. Default: "Chromosome".
#' @param legend Character of length 1. Position of the legend, passed through ggpubr::theme_pubr(). One of c("top", "bottom", "left", "right", "none"). Default: "right".
#' @param draw_points Logical of lenght 1. If TRUE, draws points in each of the drawn regions. Default: TRUE
#' @param coord_flip Logical of lenght 1. If TRUE, change the position of the axes (Y axis is position and X is chromosome). Default: FALSE
#' @param extra_info NULL or (named) list of length equal to 'regions_list', and with the same order. List of path to files or list of data.frames with two columns containing the id column of the regions contained in each set of 'regions_list' and a column with extra (discrete) information. Default: NULL.
#' @param cyto_bands NULL, data.frame or character. Path to the file (character) or imported file (dataframe) of the positions of the cytogenetic bands. Default: NULL
#' @param size_y_text NULL or numeric of length 1. If NULL, do not draw the text on the Y axis (position in bp). If numeric, the size of the text in the Y axis. Default: NULL
#'
#' @export


chromRegions <- function(chrom_sizes,
                         regions_list,
                         regions_names   = names(regions_list),
                         regions_order   = names(regions_list),
                         colors          = c("Black", "Red", "Blue"),
                         color_by        = "region",
                         chr_exclude     = c("Un", "Random", "JH", "GL", "\\."),
                         chr_order       = c(seq(from = 1, to = 23, by = 1), "X", "Y", "M", "MT", paste("chr", c(seq(from = 1, to = 23, by = 1), "X", "Y", "M", "MT"), sep = "")),
                         title           = NULL,
                         subtitle        = NULL,
                         caption         = NULL,
                         ylab            = "",
                         xlab            = "Chromosome",
                         legend          = "bottom",
                         draw_points     = TRUE,
                         coord_flip      = FALSE,
                         extra_info      = NULL,
                         cyto_bands      = NULL,
                         y_text_size     = NULL) {

  # Load packages -----
  require(dplyr)
  require(magrittr)
  library(purrr)
  require(ggplot2)
  require(ggpubr)
  require(scales)
  require(ggnewscale)

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

  ## Check color_by
  if(tolower(color_by) %in% c("regions", "region")){ color_by <- "region" }
  else if(tolower(color_by) %in% c("strands", "strand")){ color_by <- "strand" }
  else if(!(tolower(color_by) %in% c("strands", "strand","regions", "region")) & !is.null(extra_info)){ color_by <- "extra" }
  else { stop("'color_by' must be one of 'region' or 'strand' and, if an 'extra_info' list is provided, it can be 'extra'.") }

  ## No strandness --> turn color by strand to false
  if(any(regions$strand == ".")) {
    if(color_by == "strand"){
      warning("Not all regions have an accepted strand value (i.e. '+' or '-'), converting 'color_by' to 'region'.")
      color_by <- "region"
    }
  }

  ## Chromosomes to exclude
  if(!class(chr_exclude) == "character") { stop("'chr_exclude' must be a character vector with regular expressions that match the chromosomes to exclude") }

  ## Order of the chromosomes
  if(!class(chr_order) == "character") { stop("'chr_order' must be a character vector with the names of the chromosomes in the desired order") }

  ## Extra information
  if(!is.null(extra_info)){
    if(class(extra_info) != "list"){
      stop("'extra_info' must be a list with paths to files or list of dataframes with the id of the regions/peaks in each element of 'region_list' and a column for extra information")
    }
    if(all(extra_info %>% purrr::map(~class(.x)) %>% unlist()=="data.frame")) {

      # Format regions
      extrainfo <- extra_info %>%
        purrr::set_names(sets_names) %>%
        purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
        dplyr::bind_rows() %>%
        magrittr::set_colnames(c("id", "extra", "region")) %>%
        dplyr::mutate(region = factor(region, levels = sets_names))

    } else if(all(extra_info %>% purrr::map(~class(.x)) %>% unlist()=="character")) {

      # Read and format regions
      extrainfo <- purrr::map(extrainfo, read.delim, header = F) %>%
        purrr::set_names(sets_names) %>%
        purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
        dplyr::bind_rows() %>%
        magrittr::set_colnames(c("id", "extra", "region")) %>%
        dplyr::mutate(region = factor(region, levels = regions_names))

    } else if(any(extrainfo %>% purrr::map(~class(.x)) %>% unlist()) %in% c("data.frame", "character")) {
      stop("'extrainfo' must be a list with paths to regions files or data frames with the id and extra information to color the regions")
    }
  }


  # Filter chromsizes and regions ------
  chrom_filt <- chromsizes %>% dplyr::filter(!stringr::str_detect(seqnames, pattern = paste(chr_exclude, collapse = "|"))) %>%
    dplyr::mutate(seqnames = factor(seqnames, levels = chr_order))
  regions_filt <- regions %>% dplyr::filter(!stringr::str_detect(seqnames, pattern = paste(chr_exclude, collapse = "|"))) %>%
    dplyr::mutate(seqnames = factor(seqnames, levels = chr_order),
                  region = factor(region, levels = regions_order))

  # Join extra info to regions
  if(!is.null(extra_info)){
    regions_filt <- dplyr::left_join(regions_filt, extrainfo, by = c("id", "region")) %>%
      dplyr::mutate(extra = factor(extra))
  }

  # Make plot -----
  # Initialize chromosomes with or without cytobands
  if(!is.null(cyto_bands)){
    if(is.character(cyto_bands)) { cyto_bands <- read.delim(cyto_bands, header = F) }
    else if(is.data.frame(cyto_bands)) { cyto_bands <- cyto_bands }

    cytobands <- cyto_bands %>%
      magrittr::set_colnames(c("seqnames","start","end","name","gieStain")) %>%
      dplyr::filter(stringr::str_detect(string = seqnames, pattern = paste(chr_exclude, collapse = "|"),negate = T)) %>%
      dplyr::group_by(seqnames) %>%
      dplyr::mutate(size = ifelse(end == max(end), sum(end-start), NA)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(seqnames = factor(seqnames, chr_order))

    g <- ggplot() +
      geom_col(data = cytobands, mapping = aes(y = end-start, x = seqnames, group = name, fill = gieStain), width = 0.65, show.legend = F) +
      scale_fill_manual(values = c(NA, "Gray90", "Gray80", "Gray70", "Gray60"), na.value = NA) +
      ggnewscale::new_scale_fill() +
      geom_col(data = cytobands, mapping = aes(y = size, x = seqnames), color = "black", fill = NA, width = 0.7)

  } else {
    g <- ggplot() + geom_col(data = chrom_filt, mapping = aes(y = end, x = seqnames), color = "black", fill = NA, width = 0.7)
  }

  # Draw regions
  # Draw lines
  g <- g +  geom_tile(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, height = length, width = .7, fill = .data[[color_by]], color = .data[[color_by]]))
  # Draw points if desired
  if(draw_points) { g <- g + geom_point(data = regions_filt, mapping = aes(x = seqnames, y = start + length/2, fill = .data[[color_by]], color = .data[[color_by]]), size = 1) }
  # Change colors
  g <- g + scale_fill_discrete(type = colors) + scale_color_discrete(type = colors)

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
    ## Turn scientific notation of
    scale_y_continuous(expand = c(0,0), labels = function(x) format(x, scientific = FALSE)) +
    ## Theme pubr to remove border and choose the legend position
    theme_pubr(border = F, margin = T, legend = legend) +
    # Further costumization
    theme(plot.title = element_text(hjust = .5, face = "bold"), # title
          plot.subtitle = element_text(hjust = .5, face = "italic"), # subtitle
          axis.title = element_text(hjust = .5, face = "bold"), # axis titles
          legend.title = element_blank(), # remove legend title
          legend.key.size = unit(4, "mm"), #
          axis.ticks = element_blank(), axis.line = element_blank())

  if(is.null(y_text_size)){
    g <- g + theme(axis.text.y = element_blank())
  } else { g <- g + theme(axis.text.y = element_text(size = y_text_size)) }


  # Flip coordinates
  if(coord_flip) { g <- g + coord_flip() }

  # Return plot
  return(g)
}

