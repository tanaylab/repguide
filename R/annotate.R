.annoGuides <- function(guideSet, blacklist_penalty, consensus_range = NULL)
{
  guide_length <- guideSet@guide_length
  genomic_bins <- binGenome(guideSet@genome, bin_width = 500)
  blacklisted <- guideSet@blacklist
  whitelisted <- guideSet@whitelist
  kmers <- guideSet@kmers 
  anno_te <- guideSet@tes
  targets <- guideSet@targets
  cis <- guideSet@cis
  
  if (!is.null(guideSet@kmers$te_id)) # unique in case of repeated calling of annoGuides (kmers_exp may duplicate rows)
  {
    print('Overwriting existing annotation')
    mcols(kmers) <- mcols(kmers)[, c('kmer_id', 'seq', 'n_mappings', 'mismatches', 'n_mismatches', 'guide_seq')] 
  } else { print('Annotating guides') }
  
  # remove duplicated binding sites
  kmers_unique <- unique(kmers)
  mcols(kmers_unique) <- NULL
  kmers_unique$queryHits <- 1:length(kmers_unique)
  
  # Add genomic bin ID
  kmers_unique$genomic_bin <- findOverlaps(kmers_unique, genomic_bins, minoverlap = floor(guide_length/2), select = 'first')
  
  # Add black/whitelisted
  kmers_unique$blacklisted <- 1:length(kmers_unique) %in% findOverlaps(kmers_unique, blacklisted, ignore.strand = TRUE)@from
  kmers_unique$whitelisted <- 1:length(kmers_unique) %in% findOverlaps(kmers_unique, whitelisted, ignore.strand = TRUE)@from
  
  # Add distance to nearest cis
  kmers_unique$cis_dist <- NA
  hits <- distanceToNearest(kmers_unique, cis, ignore.strand = TRUE)
  kmers_unique[hits@from]$cis_dist <- mcols(hits)$distance
  
  # Add TE id, repname, and on_target
  hits_te         <- as_tibble(findOverlaps(kmers_unique, anno_te, ignore.strand = TRUE))
  hits_te$repname <- anno_te[hits_te$subjectHits]$repname
  hits_te$te_id   <- anno_te[hits_te$subjectHits]$te_id
  hits_te         <- select(hits_te, -subjectHits)
  kmers_unique    <- full_join(as_tibble(kmers_unique), hits_te, by = 'queryHits') %>% 
                     select(-queryHits) %>%
                     mutate(unique_id = ifelse(is.na(te_id), paste0('g', genomic_bin), paste0('t', te_id)),
                            on_target = as.numeric(te_id %in% guideSet@targets$te_id),
                            on_target = ifelse(on_target <= 0, -1, on_target))
  
  # Add anno to full kmer data.frame (could be improved data.table)
  kmers <- full_join(as_tibble(kmers), 
                      kmers_unique, by = c('seqnames', 'start', 'end', 'strand', 'width'))
  
  # Add GC conent
  kmer_gc <- as.numeric(letterFrequency(DNAStringSet(kmers$guide_seq), "GC"))
  kmers$gc <- round(kmer_gc / guide_length, 2) 
  
  # Export results
  guideSet@kmers <- makeGRangesFromDataFrame(kmers, keep.extra.columns = TRUE) # Could be improved
  
  # Add position on consensus for kmers
  if (length(guideSet@alignments) != 0)
  {
    guideSet <- .annoConPos(guideSet)  
  } else {
    guideSet@kmers$con_pos <- NA
  }
  
  # Set on-target 0 if mapping is outside of targeted consensus range
  if (!is.null(consensus_range))
  {
    if (sum(colnames(consensus_range) %in% c('repname', 'start', 'end')) != 3) { stop ('Consensus range requires repname, start, end columns') }
    
    guideSet@kmers <- 
      guideSet@kmers %>%
      as_tibble %>%
      right_join(consensus_range, 
                 ., 
                 by = 'repname', suffix = c('_con', '')) %>%
      mutate(outside = con_pos < start_con | con_pos > end_con,
             on_target = ifelse(!is.na(outside) & outside, 0, on_target)) %>%
      select(-outside, -start_con, -end_con) %>%
      makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  }

  # Compute guide scores
  guideSet <- .compGuideScores(guideSet, blacklist_penalty)
  
  return(guideSet)
}

.annoConPos <- function(guideSet) # Strandedness!!! To DO
{
  print ('Adding position on consensus')
  n_cores <- guideSet@.n_cores
  families <- names(guideSet@consensus)
  
  foo <- function(df)
  {
    res <- df %>% mutate(con_pos = start(pairwiseAlignment(guide_seq, consensus, type = 'global-local')@subject))
    return(res)
  }
  
  cons <- guideSet@consensus
  kmers_full <- guideSet@kmers %>% as_tibble
  
  cons_df <- 
    cons %>% 
    as.data.frame %>% 
    rownames_to_column %>% 
    as_tibble %>% 
    dplyr::rename(repname = rowname, consensus = x)
    
  kmers_slim <- 
    kmers_full %>% 
      filter(repname %in% families) %>% 
      select(kmer_id, strand, guide_seq, repname) %>% 
      distinct
    
  kmers_cons_df <- 
    left_join(kmers_slim, cons_df, by = 'repname') %>% 
    mutate(chunk = ntile(kmer_id, n_cores))
  
  kmers_anno <- plyr::ddply(kmers_cons_df, .variable = 'chunk', .fun = foo, .parallel = TRUE) %>% 
    as_tibble %>%
    select(kmer_id, strand, repname, con_pos) %>%
    right_join(., kmers_full, by = c('kmer_id', 'strand', 'repname'))
    
  guideSet@kmers$con_pos <- kmers_anno$con_pos
  
  return(guideSet)
}
