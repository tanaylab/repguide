#' Imports USCS repeatmasker file
#'
#'
#'
#' @export
importTEs <- function(filepath = 'path to UCSC RMSK file')
{
  if (!is.null(filepath))
  {
    # Import
    ucsc_rmsk <-
      fread(filepath) %>% 
      as_tibble %>%
      dplyr::rename(chrom = genoName, start = genoStart, end = genoEnd, bin = '#bin') %>%
      arrange(chrom, start, end)
    
    # Adjust 0/1-based coordinates
    ucsc_rmsk_adj <-
      ucsc_rmsk %>%
      mutate(start = start + 1)
      
    colnames(ucsc_rmsk_adj) <- tolower(colnames(ucsc_rmsk_adj))
    
    # Make GRanges
    ucsc_rmsk_gr <- GenomicRanges::makeGRangesFromDataFrame(ucsc_rmsk_adj, keep.extra.columns = TRUE)
    
    # Filter overlapping coords 
    ucsc_rmsk_filt_gr <- ucsc_rmsk_gr[!GenomicRanges::countOverlaps(ucsc_rmsk_gr) > 1, ]
    
    # Add unique TE ID
    ucsc_rmsk_filt_gr$te_id <- 0:(length(ucsc_rmsk_filt_gr) - 1)
  }
  
  return(ucsc_rmsk_filt_gr)
}

importKmers <- function(filepath = 'path to bowtie mapped kmers')
{
  # Import
  kmers <- 
    fread(filepath) %>% 
    as_tibble %>%
    dplyr::rename(kmer_id = V1, strand = V2, chrom = V3, start = V4, seq = V5, n_mappings = V6, mismatches = V7) %>%
    mutate(start = start + 1, # adjust for difference between 0-based bowtie and 1-based GRanges
           end = start + nchar(seq) - 1,
           n_mappings = n_mappings + 1,
           n_mismatches = stringr::str_count(mismatches, stringr::fixed(':'))) %>%
    replace_na(list(n_mismatches = 0))
  
  if (nrow(kmers) > 5e7 & parallel::detectCores() > 20)
  {
    print ('Parallel GRanges conversion')
    kmers_gr <- makeGRangesFromDataFramePar(kmers, keep.extra.columns = TRUE)  
  } else {
    kmers_gr <- makeGRangesFromDataFrame(kmers, keep.extra.columns = TRUE)
  }
 
  return(kmers_gr)
}