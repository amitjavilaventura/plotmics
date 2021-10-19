### circleRegions()
### ===============================================================================================

#' @title circleRegions
#' @author amitjavilaventura
#'
#' @description
#' Function that draws the chromosomes and the position of the desired regions in a circular plot
#' It draws each regions set provided as input in separate circles, so you can input chromosome sizes for different assemblies.
#' It also allows to connect "paired" regions across each provided regions set and circle.
#' It allows to color regions by "region", "strand" or "extra", which is provided through extra_info.
#'
#' @param chromsizes_sets (Named) List of character vectors or dataframes. List with paths to chrom.sizes file or list of data frames with the names of the chromosomes and their sizes. The file must not contain column names. Length and order must be the same as 'regions_sets'.
#' @param regions_sets (Named) List of character vectors or dataframes. List with paths to BED/TSV files or dataframes with chromosome names, start, end, region id, length and strand. Strand must be "+", "-" or ".".  Length and order must be the same as 'chromsizes_sets'.
#' @param chr_exclude Character. Regular expressions to match the chromosomes names to exclude . Default: c("Un", "Random", "JH", "GL", "\\.").
#' @param chr_order Character. Chromosome names written in the order to be plotted. The names must matchthe chromosome names in chrom_sizes and regions. Default: c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","MT")
#' @param chr_label Character of length 1. Color of the chromosome labels. If all the elements in chromsizes_sets are equal, only prints the labels once. Default: "Black"
#' @param chr_line Logical of length 1. If TRUE, it draws a dashed line that separates a the chromosomes. This line is done only for the first chromsizes set.
#' @param sets_names Character of length equal to 'chromsizes_sets' and 'regions_sets'. Names of the input regions in the correct order. They are used to name the regions_sets, the chromsizes_sets and the extra_info lists. Default: names(regions_sets).
#' @param colors Character. Colors of the input regions. Depends on 'color_by': if it is 'region', the length must be equal or greater to the lenght of 'regions_sets'; if it is 'strand' the lenght must be 2 or more; if it is 'extra', the length must be equal or greater to the groups provided in extra_info.
#' @param color_by Character of length 1. One of 'region', 'strand' or 'extra'. If 'region', colors by the regions in 'regions_sets'. If 'strand', colors by strand. If 'extra' colors by the grouping provided in 'extra_info' (if provided).
#' @param title Character of length 1 or NULL. Title of the plot. Default: NULL
#' @param subtitle Charachter of length 1 or NULL. Subtitle of the plot. Default: NULL.
#' @param caption Character of length 1, TRUE or NULL. Caption to be written in the bottom-right corner. If TRUE, it will write the number of regions in the input; if charachter, it will write the written caption. Default: NULL.
#' @param legend Character of length 1. Position of the legend, passed through ggpubr::theme_pubr(). One of c("top", "bottom", "left", "right", "none"). Default: "bottom".
#' @param draw_points Logical of lenght 1. If TRUE, draws points in each of the drawn regions. Default: TRUE.
#' @param paired Logical of lenght 1. If TRUE, draws a line between the paired regions in each set/circle. Default: FALSE.
#' @param paried_color Character of length 1. Color of the lines that connect the paired regions.
#' @param extra_info NULL or (named) list of length equal to 'chromsizes_sets' and 'regions_sets', and with the same order. List of path to files or list of data.frames with two columns containing the id column of the regions contained in each set of 'regions_sets' and a column with extra (discrete) information. Default: NULL.
#'
#' @export

circleRegions <- function(chromsizes_sets,
                          regions_sets,
                          chr_exclude     = c("Un", "Random", "JH", "GL", "\\."),
                          chr_order       = c(seq(from = 1, to = 23, by = 1), "X", "Y", "M", "MT", paste("chr", c(seq(from = 1, to = 23, by = 1), "X", "Y", "M", "MT"), sep = "")),
                          chr_label       = "Black",
                          chr_line        = FALSE,
                          sets_names      = names(regions_sets),
                          colors          = c("Black", "Red", "Blue"),
                          color_by        = "region", # one of "region", "extra", "strand"
                          title           = NULL,
                          subtitle        = NULL,
                          caption         = NULL,
                          legend          = "bottom",
                          draw_points     = TRUE,
                          paired          = FALSE,
                          paired_color    = "blue",
                          extra_info      = NULL
){


  # Load packages -----
  require(dplyr)
  require(magrittr)
  library(purrr)
  require(ggplot2)
  require(ggpubr)

  # Check that inputs are OK and read data --------------------------------------------------------

  ## Chromosomes to exclude
  if(!class(chr_exclude) == "character") { stop("'chr_exclude' must be a character vector with regular expressions that match the chromosomes to exclude") }

  ## Order of the chromosomes
  if(!class(chr_order) == "character") { stop("'chr_order' must be a character vector with the names of the chromosomes in the desired order") }


  ## Colors to factor
  colors <- factor(colors, levels = colors)

  ## Length of regions and chromsizes sets must be equal
  if(length(regions_sets) != length(chromsizes_sets)){
    stop( "'regions_sets' and 'chromsizes_sets' must be lists with the same length.
          They should have the desired regions to plot and the chromosome size information of the corresponding assembly, respectively." )
  }

  ## Length of regions names/colors are equal to length of regions list
  if(length(regions_sets) != length(sets_names)){
    stop("'sets_names' must be a character vector with the names of the sets in 'chromsizes_sets' and 'regions_sets' and their length must be the same.")
  }
  if(length(chromsizes_sets) != length(sets_names)){
    stop("'sets_names' must be a character vector with the names of the sets in 'chromsizes_sets' and 'regions_sets' and their length must be the same.")
  }

  ## Chrom sizes
  if(!class(chromsizes_sets) == "list") { stop("'chromsizes_sets' must be a list with paths to chrom.sizes files or list of data frames with the names and size of the chromosomes") }
  if(all(chromsizes_sets %>% purrr::map(~class(.x)) %>% unlist()=="data.frame")) {

    # Format regions
    chromsizes <- chromsizes_sets %>%
      purrr::set_names(c(sets_names)) %>%
      purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
      dplyr::bind_rows() %>%
      magrittr::set_colnames(c("seqnames", "size", "region")) %>%
      dplyr::mutate(region = factor(region, levels = sets_names))

  } else if(all(chromsizes_sets %>% purrr::map(~class(.x)) %>% unlist()=="character")) {

    # Read and format regions
    chromsizes <- purrr::map(chromsizes_sets, read.delim, header = F) %>%
      purrr::set_names(c(sets_names)) %>%
      purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
      dplyr::bind_rows() %>%
      magrittr::set_colnames(c("seqnames", "size", "region")) %>%
      dplyr::mutate(region = factor(region, levels = sets_names))

  } else if(any(chromsizes_sets %>% purrr::map(~class(.x)) %>% unlist()) %in% c("data.frame", "character")) {
    stop("'regions_sets' must be a list with paths to regions files or data frames with the positions of the regions")
  }

  ## Regions
  if(!class(regions_sets) == "list") { stop("'regions_sets' must be a list with paths to regions files or list of data frames with the positions (BED6-like format) of the regions") }
  if(all(regions_sets %>% purrr::map(~class(.x)) %>% unlist()=="data.frame")) {

    # Format regions
    regions <- regions_sets %>%
      purrr::set_names(c(sets_names)) %>%
      purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
      dplyr::bind_rows() %>%
      magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand", "region")) %>%
      dplyr::mutate(region = factor(region, levels = sets_names))

  } else if(all(regions_sets %>% purrr::map(~class(.x)) %>% unlist()=="character")) {

    # Read and format regions
    regions <- purrr::map(regions_sets, read.delim, header = F) %>%
      purrr::set_names(c(sets_names)) %>%
      purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
      dplyr::bind_rows() %>%
      magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand", "region")) %>%
      dplyr::mutate(region = factor(region, levels = sets_names))

  } else if(any(regions_sets %>% purrr::map(~class(.x)) %>% unlist()) %in% c("data.frame", "character")) {
    stop("'regions_sets' must be a list with paths to regions files or data frames with the positions of the regions")
  }

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

  ## Extra information
  if(!is.null(extra_info)){
    if(class(extra_info) != "list"){
      stop("'extra_info' must be a list with paths to files or list of dataframes with the id of the regions/peaks in each element of 'region_list' and a column for extra information")
    }
    if(all(extra_info %>% purrr::map(~class(.x)) %>% unlist()=="data.frame")) {

      # Format regions
      extrainfo <- extra_info %>%
        purrr::set_names(c(sets_names)) %>%
        purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
        dplyr::bind_rows() %>%
        magrittr::set_colnames(c("id", "extra", "region")) %>%
        dplyr::mutate(region = factor(region, levels = sets_names))

    } else if(all(extra_info %>% purrr::map(~class(.x)) %>% unlist()=="character")) {

      # Read and format regions
      extrainfo <- purrr::map(extrainfo, read.delim, header = F) %>%
        purrr::set_names(c(sets_names)) %>%
        purrr::imap(~dplyr::mutate(.x, region = .y)) %>%
        dplyr::bind_rows() %>%
        magrittr::set_colnames(c("id", "extra", "region")) %>%
        dplyr::mutate(region = factor(region, levels = sets_names))

    } else if(any(extrainfo %>% purrr::map(~class(.x)) %>% unlist()) %in% c("data.frame", "character")) {
      stop("'extrainfo' must be a list with paths to regions files or data frames with the id and extra information to color the regions")
    }
  }



  # Filter and format chromsizes and regions ------------------------------------------------------
  chrom_filt <- chromsizes %>%
    # Filter chromosomes to exclude undesired ones.
    dplyr::filter(!stringr::str_detect(seqnames, pattern = paste(chr_exclude, collapse = "|"))) %>%
    # Convert seqnames to factor with desired order of chromosomes
    dplyr::mutate(seqnames = factor(seqnames, levels = chr_order)) %>%
    # Arrange by chromsosomes
    dplyr::arrange(seqnames) %>%
    # Calculate sizes
    dplyr::group_by(region) %>%
    dplyr::mutate(cum_size = cumsum(as.numeric(size)),
                  prev_size = lag(as.numeric(size)),
                  prev_size = ifelse(is.na(prev_size), 0, prev_size),
                  total_size = sum(as.numeric(size))) %>%
    dplyr::ungroup()

  regions_filt <- regions %>%
    # Filter regions to exclude undesired chromosomes.
    dplyr::filter(!stringr::str_detect(seqnames, pattern = paste(chr_exclude, collapse = "|"))) %>%
    # Convert seqnames and regions to factor with the desired order.
    dplyr::mutate(seqnames = factor(seqnames, levels = chr_order),
                  region   = factor(region, levels = sets_names))

  # Join extra info to regions
  if(!is.null(extra_info)){
    regions_filt <- dplyr::left_join(regions_filt, extrainfo, by = c("id", "region")) %>%
      dplyr::mutate(extra = factor(extra))
  }

  ## Seqnames in regions = seqnames in chrom sizes
  if(!all(regions_filt$seqnames %in% chrom_filt$seqnames)) { stop("The chromosome names in 'regions' must be equal to the chromosome names in 'chrom_sizes' (after filtering chromosomes)") }


  # Left join of chrom_filt and regions_filt -----
  regions_chrom <- dplyr::left_join(chrom_filt, regions_filt, by = c("seqnames", "region")) %>%
    # Group by chromosome
    # remove all the redundant "size" values
    # remove all the redundant "prev_size" values
    # ungroup
    # convert sequnames to factor and arrange by seqnames
    dplyr::group_by(seqnames, region) %>%
    dplyr::mutate(row = row_number(),
                  size2 = max(size),
                  size = if_else(row == 1, size, NULL),
                  prev_size2 = max(prev_size),
                  prev_size = if_else(row == 1, prev_size, NULL),
                  num_regions = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = group_indices(., .dots = "region")) %>%
    dplyr::mutate(seqnames = factor(seqnames, levels = chr_order)) %>%
    dplyr::arrange(seqnames) %>%
    dplyr::mutate(region = factor(region, levels = sets_names))

  regions_chrom2 <- regions_chrom %>% dplyr::mutate(start = start/total_size,
                                                    end = end/total_size,
                                                    length = length/total_size,
                                                    size = size/total_size,
                                                    size2 = size2/total_size,
                                                    cum_size = cum_size/total_size,
                                                    prev_size = prev_size/total_size,
                                                    prev_size2 = prev_size2/total_size,
                                                    total_size = total_size/total_size)

  if(length(chromsizes_sets)>1){
    if(identical(chromsizes_sets[[1]], chromsizes_sets[[2]])){
      regions_chrom_filt <- regions_chrom2 %>%
        dplyr::mutate(seqnames = if_else(region == sets_names[1], seqnames, NULL))
    } else{
      regions_chrom_filt <- regions_chrom2
    }
  } else {
    regions_chrom_filt <- regions_chrom2
  }

  # DRAW PLOT -------------------------------------------------------------------------------------
  # Initialize plot
  g <- ggplot(data = regions_chrom2) +
    # Draw chromosomes
    geom_col(aes(x = 11-group*2, y = size, group = seqnames), color = "black", fill = NA, width = .5, na.rm = T, position = "stack")

  # Draw line to separate chromosomes
  if(chr_line){
    g <- g + geom_hline(data = regions_chrom2 %>% dplyr::filter(region == sets_names[1]), aes(yintercept = total_size-cum_size+size2), linetype = 2, size = 0.1, na.rm = T)
  }
  # Draw segments and points for each region
  g <- g + geom_segment(data = regions_chrom2, aes(x = 11-group*2-.25, xend = 11-group*2+.25, y = total_size-cum_size+size2-start,
                                                  yend = total_size-cum_size+size2-end, color = .data[[color_by]]), na.rm = T)
  # Draw points for each region
  if(draw_points){ g <- g + geom_point(data = regions_chrom2, aes(x = 11-group*2, y = total_size-cum_size+size2-start-length/2, color = .data[[color_by]], group = id), size = .4, na.rm = T) }
  else{  g <- g + geom_point(data = regions_chrom2, aes(x = 11-group*2, y = total_size-cum_size+size2-start-length/2, group = id), size = .4, na.rm = T, color = NA) }

  # Draw paired lines
  if(paired){
    g <- g + geom_line(data = regions_chrom2, aes(x = 11-group*2, y = total_size-cum_size+size2-start-length/2, group = id), size = .2, color = paired_color, alpha = .7)
  }

  # Write chromosome labels
  g <- g  + geom_text(data = regions_chrom_filt, aes(x = 11-group*2+1, y = size, label = seqnames, group = seqnames),
                      position = position_stack(vjust = .5), show.legend = F,  na.rm = T, size = 2.5, hjust = .5, color = chr_label, alpha = 0.8)

  # Change colors and remove NAs from legend
  g <- g + scale_color_manual(values = colors, na.translate = F)

  # Make plot circular
  g <- g + coord_polar(theta = "y", direction = -1, start = 0) + xlim(0,10)


  # Add labels and titles
  if(!is.null(title)){ g <- g + ggtitle(label = title, subtitle = subtitle) }
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
  # Customize plot theme
  g <- g +
    theme_pubr(legend = legend) + # Add theme pubr from ggpubr with options for legend
    theme(axis.text = element_blank(), # Remove axis text
          axis.ticks = element_blank(), # Remove axis ticks
          axis.title = element_blank(), # Remove axis titles
          axis.line = element_blank(), # Remove axis lines
          legend.title = element_blank(), # Remove legend title
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5, face = "italic"))

  # Return plot
  return(g)
}
