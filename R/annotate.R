annoGuides <- function(guideSet)
{
  guide_length <- guideSet@guide_length
  genomic_bins <- binGenome(guideSet@genome, bin_width = 250)
  kmers <- guideSet@kmers 
  anno_te <- guideSet@tes
  targets <- guideSet@targets
  cis <- guideSet@cis
  PAM <- guideSet@PAM
  
  if (!is.null(guideSet@kmers$te_id)) # unique in case of repeated calling of annoGuides (kmers_exp may duplicate rows)
  {
    print('Overwriting existing annotation')
    kmers <-
      kmers %>% 
      as_tibble %>%
      select(-cis_dist, -gc, -Sbind, -Scis, -Soff, -Son, -te_id, -unique_id, -genomic_bin) %>%
      distinct %>%
      makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  } else { print('Annotating guides') }
  
  
  # Add genomic bin ID
  kmers$genomic_bin <- GenomicRanges::findOverlaps(kmers, genomic_bins, minoverlap = floor(guide_length/2), select = 'first')
    
  kmers$te_id <- as.numeric(NA)
  kmers$repname <- as.character(NA)
  kmers$on_target <- as.character(NA)
	
  # Find Overlaps between kmers and te coords
  hits_te <- GenomicRanges::findOverlaps(kmers, anno_te) %>% # Sometimes more than 1 TE hit per binding site gs@
    as_tibble %>%
    right_join(., tibble(queryHits = 1:length(kmers)))
	kmers_exp <- kmers[hits_te$queryHits]
  hits_te_exp <- hits_te %>% mutate(queryHits = 1:nrow(hits_te))
  # Add universal TE ID, family name, and the #bases into the te locus (required to map onto MSA)
  
  kmer_indeces <- hits_te_exp %>% filter(!is.na(subjectHits)) %>% pull(queryHits)
  te_indeces <- na.omit(hits_te_exp$subjectHits)
  
  kmers_exp[kmer_indeces]$te_id <- anno_te[te_indeces]$te_id
  kmers_exp[kmer_indeces]$repname <- anno_te[te_indeces]$repname
  kmers_exp$unique_id <- ifelse(is.na(kmers_exp$te_id), kmers_exp$genomic_bin, kmers_exp$te_id)
  
  kmers_exp$on_target <- as.numeric(1:length(kmers_exp) %in% GenomicRanges::findOverlaps(kmers_exp, guideSet@targets)@from & kmers_exp$repname %in% guideSet@families)
  kmers_exp$on_target[kmers_exp$on_target <= 0] <- -1
  
  # Add distance to nearest cis
  kmers_exp$cis_dist <- NA
  hits <- distanceToNearest(kmers_exp, cis)
  kmers_exp$cis_dist[hits@from] <- mcols(hits)$distance
  
  # Add GC conent
  kmer_gc <- as.numeric(letterFrequency(DNAStringSet(substring(kmers_exp$guide_seq, 1, guide_length)), "GC"))
  kmers_exp$gc <- round(kmer_gc / guide_length, 2) 
  
  guideSet@kmers <- kmers_exp
  
  # Add position on consensus for kmers
  if (length(guideSet@alignments) != 0)
  {
    guideSet <- .annoConPos(guideSet)  
  } else {
    guideSet@kmers$con_pos <- NA
  }

  # Compute guide scores
  guideSet <- .compGuideScores(guideSet)
  
  return(guideSet)
}

.annoConPos <- function(guideSet) # Strandedness!!! To DO
{
  print ('Adding position on consensus')
  n_cores <- guideSet@.n_cores
  families <- names(guideSet@consensus)
  
  foo <- function(df)
  {
    res <- df %>% mutate(con_pos = start(pairwiseAlignment(seq, consensus, type = 'global-local')@subject))
    return(res)
  }
  
  cons <- guideSet@consensus
  kmers_full <- guideSet@kmers %>% as_tibble
  
  cons_df <- 
    cons %>% 
    as.data.frame %>% 
    rownames_to_column %>% 
    as_tibble %>% 
    rename(repname = rowname, consensus = x)
    
  kmers_slim <- 
    kmers_full %>% 
      filter(repname %in% families) %>% 
      select(kmer_id, strand, seq, repname) %>% 
      distinct
    
  kmers_cons_df <- 
    left_join(kmers_slim, cons_df) %>% 
    mutate(chunk = ntile(kmer_id, n_cores))
  
  kmers_anno <- plyr::ddply(kmers_cons_df, .variable = 'chunk', .fun = foo, .parallel = TRUE) %>% 
    as_tibble %>%
    select(kmer_id, strand, repname, con_pos) %>%
    right_join(., kmers_full)
    
  guideSet@kmers$con_pos <- kmers_anno$con_pos
  
  return(guideSet)
}
