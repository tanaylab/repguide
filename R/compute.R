.compGuideScores <- function(guideSet) # Sbind can be interprated as the likelyhood of binding to a site (0 - 1)
{
  print('Computing guide scores')
  kmers <- as_tibble(guideSet@kmers) %>% select(-matches('Sbind'))
  PAM <- guideSet@PAM
  guide_length <- guideSet@guide_length + nchar(PAM)
  max_mismatches <- max(kmers$n_mismatches)
  
  ###################
  # Calculate Sbind #
  ###################
  if (max_mismatches > 0)
  {
    mismatch_universe <- .compMismatchUniverse(guide_length, max_mismatches) # guidelength must be 1 based
    mismatch_universe_scored <- .compMismatchScore(mismatch_universe, guide_length) # guidelength must be 1 based

    # get kmer mismatches
    mismatches <- stringr::str_extract_all(kmers$mismatches, '[0-9]+', simplify = TRUE) # could be parallelized!, quite fast already
    colnames(mismatches) <- c('first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eights')[1:ncol(mismatches)]
    
    # add mismatches to kmer data.frame
    kmers <- bind_cols(kmers, 
                      (as_tibble(mismatches) %>% mutate_all(as.numeric)))
    
    # add mismatch score (Sbind), could be improved with data.table
    kmers <- 
      left_join(kmers, 
                mismatch_universe_scored) %>%
      mutate(Sbind = ifelse(n_mismatches == 0, 1, Sbind))
  } else {
    kmers$Sbind <- 1
  }  
  
  ###################
  # Calculate Scis #
  ################### 

  # Compute cis score (Boltzmann sigmoid)
  .scoreCis <- function(x, slope = 5, midpoint = 2000)
  {
    s <- 1 + ((0.01 - 1) / (1 + exp((midpoint - x)/(slope * 100)))) # Boltzmann
    #s <- pmin(100, s)
    return(s)  
  } 
  kmers$Scis = .scoreCis(kmers$cis_dist)
  kmers <- kmers %>% mutate(Scis = ifelse(is.na(Scis), 0.01, Scis))
  #######################################
  
  # Compute Soff
  kmers$Soff <- kmers$Sbind * kmers$Scis
  kmers <- kmers %>% mutate(Soff = ifelse(on_target >= 0, 0, Soff))
  #kmers$Soff <- round(kmers$Soff, 2)
  
  # Compute Son
  kmers$Son <- kmers$Sbind
  kmers <- kmers %>% mutate(Son = ifelse(on_target <= 0, 0, Son))
  #kmers$Son <- round(kmers$Son, 2)
  
  # Add results to guideSet
  guideSet@kmers$Sbind <- kmers$Sbind
  guideSet@kmers$Scis  <- kmers$Scis
  guideSet@kmers$Soff  <- kmers$Soff
  guideSet@kmers$Son   <- kmers$Son

  return(guideSet)
}

.compMismatchUniverse <- function(guide_length, max_mismatches)
{
  column_labels <- c('first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eights')
  
  liste <- list()
  for (n_mism in 1:max_mismatches)
  {
    liste[[n_mism]] <-
      combn(0:(guide_length -1), n_mism) %>% 
      t %>%
      as.data.frame
  }
  
  mismatch_universe <- as_tibble(do.call(bind_rows, liste))
  colnames(mismatch_universe) <- column_labels[1:max_mismatches]
  return(mismatch_universe)
}

.compMismatchScore <- function(mismatch_universe, guide_length)
{
  guide_length <- guide_length - 1
  
  # Formula adapted from Breaking-Cas: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987939/
  .scoring <- function(weight, d, m, guide_length)
  {
    s <- weight * (1 / (((guide_length - d) / guide_length) * 1 + 1)) * (1 / (m * m))
    return(s)
  }
  
  # clean data.frame in case of repeated usage
  mismatch_universe <- mismatch_universe %>% select(-which(colnames(mismatch_universe) %in% c('m', 'max', 'min', 'd', 'weight', 'Sbind')))
  
  # total mismatches
  mismatch_universe$m <- rowSums(!is.na(mismatch_universe %>% select(-which(colnames(mismatch_universe) %in% c('m')))))
  
  # approximate of average pairwise distance between mismatches
  mismatch_universe <- 
    mismatch_universe %>% 
    mutate(max = apply(mismatch_universe[, !colnames(mismatch_universe) %in% 'm'], 1, max, na.rm = TRUE),
           min = apply(mismatch_universe[, !colnames(mismatch_universe) %in% 'm'], 1, min, na.rm = TRUE),
           d = ifelse(m > 1, (max - min) / (m - 1), 0))               
  
  # custom weighting function based on position in guide
  mismatch_universe <-
    mismatch_universe %>%
    mutate(weight = apply(mismatch_universe[, !colnames(mismatch_universe) %in% c('m', 'max', 'min', 'd', 'weight')], 1, function(x)
                    {
                      x <- 1 - x / guide_length
                      weight <- prod(x, na.rm = TRUE)
                    })
          )
  
  mismatch_universe <- 
    mismatch_universe %>% 
    mutate(Sbind = .scoring(weight, d, m, guide_length)) %>%
    select(-m, -max, -min, -d, -weight)

  return(mismatch_universe)
}


.compCombinations <- function(guideSet,
                              n_guides_max = 5,
                              method = 'cluster')
{
  if (method == 'cluster')
  {
    if (is.null(guideSet@kmers$best)) { stop ('No best Kmer selection found') }
    
    kmers <- as_tibble(guideSet@kmers)
    kmers_best_per_clust <- kmers %>% filter(best) %>% pull(kmer_id) %>% unique
    
    # Calc all combinations for up to max guides
    combinatorics_df <- 
      lapply(
        sapply(1:n_guides_max, function(x) 
        { 
          a <- combn(kmers_best_per_clust, x); b <- as_tibble(a); return(b) 
        }), 
      function(y) 
      { 
        y %>% gather('combi_id', 'kmer_id') %>% mutate(n_guides = nrow(y)) 
      }) %>% 
      do.call(rbind, .) %>%
      mutate(combi_id = paste0(n_guides, combi_id)) %>%
      data.table::as.data.table(key = 'kmer_id')
    # ===========================================
      
    kmers_slim <- # max Son/Soff per kmer and unique_id
      kmers %>%
      filter(kmer_id %in% kmers_best_per_clust) %>%
      select(kmer_id, unique_id, Son, Soff) %>% 
      distinct %>% 
      data.table::as.data.table(key = 'kmer_id')
      
    # Calculate on/off targets per combo
    combination_hits <- kmers_slim[combinatorics_df, on = 'kmer_id', allow.cartesian = TRUE, nomatch = 0]
    combination_stats <- combination_hits[, .(Son = max(Son), 
                                              Soff = max(Soff)),
                                              by = c('combi_id', 'unique_id')
                                         ][, .(Son_tot = sum(Son), 
                                               Soff_tot = sum(Soff),
                                               on_tot = length(unique_id[Son > 0]),
                                               off_tot = length(unique_id[Soff > 0])),
                                               by = 'combi_id'] %>%
      as_tibble %>%
      left_join(combinatorics_df) %>%
      arrange(n_guides) %>%
      nest(kmer_id, .key = 'kmer_id') %>%
      mutate(enr = ifelse(Soff_tot == 0, (Son_tot) / (Soff_tot + 0.01), Son_tot / Soff_tot)) # adding a pseudo-count to avoid Inf)
    # ===================================
    guideSet@combinations <- combination_stats
  } 
   
  if (method == 'rank')
  {
    a <- kmers
    kmers_best <- list()
    for (i in 1:n_guides)
    {
     best_kmer <-
      a %>%
      count(kmer_id, on_target) %>%
      spread(on_target, n, fill = 0, sep = '_') %>%
      filter(on_target_FALSE <= max_off_target) %>%
      top_n(1, on_target_TRUE) %>%
      top_n(-1, on_target_FALSE) %>%
      sample_n(1) %>%
      pull(kmer_id)
      
    covered_loci <- a %>% filter(kmer_id == best_kmer & on_target) %>% pull(te_id)
    
    a <-
      a %>%
      mutate(on_target = ifelse(te_id %in% covered_loci, FALSE, on_target))
    
      kmers_best[[i]] <- best_kmer 
    print(i)
    }
    
    kmers_final <- unlist(kmers_best)
  }
  return(guideSet)
}

                         


compStats <- function(kmers)
{
  stats <-
    kmers %>% 
      filter(selected) %>% 
      summarise(total_hits = sum(on_target), 
                unique_hits = sum(on_target[!duplicated(te_id)]), 
                total_off = sum(!on_target), 
                total_blacklisted = sum(blacklisted))
                
  print(stats)
}

# Check multi processor performance
# Think of providing substituiton matrix 
.compMSA <- function(seqs,
                     max_gap_freq = 0.8,
                     iterations,
                     refinements,
                     kmer_length = 7,
                     n_clust = 10,
                     clust_perc = 1,
                     seed = NULL) # ultrafast but rough MSA
{
  if (class(seqs) == 'DNAStringSet') 
  {
    if (is.null(names(seqs))) 
    {
      names(seqs) <- 1:length(seqs)
    }
  }
  
  if (class(seqs) == 'GRanges')
  {
    seqs_gr <- seqs
    seqs <- seqs_gr$seq
    names(seqs) <- seqs_gr$te_id
  }
  
  if (length(seqs) > 10000) { seqs <- sample(seqs, 10000) }
  
  seqs <- DNAStringSet(seqs)
  seqs <- seqs[width(seqs) < 9e4]
  
  seqs_bin <- ape::as.DNAbin(seqs)
  clusts <- 
    tibble(te_id = names(seqs), 
                 clust = cutree(fastcluster::hclust(tgstat::tgs_dist(kmer::kcount(seqs_bin, k = kmer_length)), 'ward.D2'), n_clust))
  clust_sizes <- count(clusts, clust)
    
  # sample ids per clust proportional to clust size
  ids_sel <- 
    clusts %>% 
    nest(te_id) %>% 
    mutate(n = ceiling(clust_sizes$n / 100 * clust_perc),
           samp = purrr::map2(data, n, sample_n)) %>% 
    select(clust, samp) %>% 
    unnest %>% 
    pull(te_id)

  message(paste0('Aligning ', length(ids_sel), ' sequences'))
  
  alignment <- DECIPHER::AlignSeqs(seqs[ids_sel], 
                                   #processors = n_cores, 
                                   iterations = iterations, 
                                   refinements = refinements, 
                                   useStructures = FALSE,
                                   verbose = FALSE)
  alignment_wo_gaps <- .rmGaps(alignment, max_gap_freq = max_gap_freq)
  
  return(alignment_wo_gaps)
}





