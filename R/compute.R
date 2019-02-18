#' Cluster guideRNAs
#'
#' Clusters guides based on their target binding profile
#'
#' @param guideSet guideSet containing guide mappings
#' @param min_Son Numeric from 0 through 1. Only considers genomic target binding sites above \code{min_Son} score.
#' @param n_clust Integer from 1 to 20. Number of clusters to group guides into. Passed to \code{cutree()} function.
#' @return guideSet object with clustered guides.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13A')
#' gs <- addGuides(gs, guide_length = 16)
#' gs <- compClusts(gs, min_Son = 0.25, n_clust = 10)
#' gs <- plotGuides(gs)
#' }
#' @seealso [addGuides()], and [plotGuides()]
#' @export
clustGuides <- function(guideSet, 
                        min_Son = 0,
                        n_clust = 15)
{
  if(n_clust > 20) { stop('Maximal 20 clusters currently supported') }
  message('Clustering kmers')
  set.seed(guideSet@.seed)
  kmers <- as_tibble(guideSet@kmers) %>% select(-matches('kmer_clust|te_clust'))
  kmers_filt <- 
    kmers %>%
    filter(Son > min_Son & on_target == 1) %>%
    filter(valid)
    
  if (nrow(kmers_filt) == 0) { stop ('No valid guides found, try relaxing selection parameters of addGuides function') }
  if (length(unique(kmers_filt$kmer_id)) < n_clust) { stop ('Less valid guides than number of clusters') }
  
  mat_full <- 
    kmers_filt %>%
    select(kmer_id, te_id, Son) %>%
    .tidyToSparse()
    
  #mat_full = log2(mat_full+1)
   
  # mat_slim <- kmers %>%
    # filter(on_target >= 0) %>%
    # mutate(on_target = on_target * Sbind) %>%
    # select(kmer_id, te_id, on_target) %>%
    # tidyToSparse()  

  print(paste0('Clustering ', nrow(mat_full), ' kmers into ', n_clust, ' groups'))
  kmer_cors <- as.matrix(qlcMatrix::cosSparse(t(mat_full)))
  #kmer_cors <- tgs_cor(as.matrix(t(mat_full)), spearman = TRUE)
  kmer_clusts <- tibble(kmer_id = as.numeric(rownames(mat_full)),
                        kmer_clust = as.numeric(cutree(fastcluster::hclust(tgstat::tgs_dist(kmer_cors), 'ward.D2'), n_clust)))
                          
  
  if (ncol(mat_full) > 2e4) 
  {
    message ('Downsampling the matrix')
    #vars <- matrixStats::colVars(mat_full)
    vars <- apply(mat_full, 2, var)
    mat_full <- mat_full[, tail(order(vars), 2e4)]
  } else {
    mat_full <- mat_full
  }
  
  print(paste0('Clustering ', ncol(mat_full), ' loci into ', n_clust, ' groups'))
  loci_cors <- as.matrix(qlcMatrix::cosSparse(mat_full))
  loci_clusts <- tibble(te_id = as.numeric(colnames(mat_full)),
                        te_clust = as.numeric(cutree(fastcluster::hclust(tgstat::tgs_dist(loci_cors), 'ward.D2'), n_clust)))   

  kmers <- left_join(kmers, kmer_clusts, by = 'kmer_id') %>% left_join(., loci_clusts, by = 'te_id')                        
  guideSet@kmers$kmer_clust <- kmers$kmer_clust
  guideSet@kmers$te_clust <- kmers$te_clust
  
  return(guideSet)
}

.compGuideScores <- function(guideSet,
                             blacklist_penalty = 10) # Sbind can be interprated as the likelyhood of binding to a site (0 - 1)
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
      suppressMessages(left_join(kmers, 
                mismatch_universe_scored)) %>%
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
  kmers <- 
    kmers %>% 
    mutate(Soff = ifelse(on_target >= 0, 0, Soff),
                         Soff = ifelse(blacklisted, Soff * blacklist_penalty, Soff),
                         Soff = ifelse(is.na(Soff), 0, Soff), # if Soff(0) x blacklist_penalty(Inf) = NaN
                         Soff = ifelse(whitelisted, 0, Soff))
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
      left_join(., combinatorics_df, by = 'combi_id') %>%
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

# Check multi processor performance
# Think of providing substituiton matrix 
.compMSA <- function(seqs,
                     max_gap_freq = 0.8,
                     iterations,
                     refinements,
                     kmer_length = 7,
                     n_clust = 10,
                     clust_perc = 1,
                     seed = 19) # ultrafast but rough MSA
{
  set.seed(seed)
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

.compGreedy <- function(guideSet,
                        iterations = 10)
{
  set.seed(guideSet@.seed)
  kmers <- as_tibble(guideSet@kmers) %>% filter(valid) %>% mutate(kmer_id = as.character(kmer_id))
  combinations <- guideSet@combinations %>% unnest
  max_n_guides <- max(combinations$n_guides)
  
  if (max_n_guides > 1)
  {
    message('Running greedy optimization')
    # Create the score matrix
    mat <- 
      kmers %>%
      mutate(score = ifelse(on_target == 1, Son, Soff)) %>%
      select(kmer_id, unique_id, score) %>%
      .tidyToSparse() %>%
      as.matrix
    on_indeces <- which(colnames(mat) %in% unique(kmers$unique_id[kmers$on_target == 1]))
    off_indeces <- which(colnames(mat) %in% unique(kmers$unique_id[kmers$on_target == -1]))
    mat[mat > 1] <- 1 # ceil indv loci score at 1 
    
    report <- foreach::foreach (nguides = 2:max_n_guides, .combine = rbind) %dopar%
    {
      # Get current best combination and calc stats
      kmers_best <- combinations %>% filter(n_guides == nguides & best) %>% pull(kmer_id) %>% as.character()
      #kmers_best <- sample(rownames(mat), nguides)
      score_all_best <- matrixStats::colMaxs(mat[kmers_best, ])
      score_on_best  <- sum(score_all_best[on_indeces])
      score_off_best <- sum(score_all_best[off_indeces])
      
      # Run greedy search
      df <- tibble(iterations = 1:iterations,
                   Son_tot = score_on_best,
                   Soff_tot = score_off_best,
                   n_guides = nguides,
                   kmer_id = list(kmers_best))
           
      for (i in 1:iterations)
      {
        #print(i)
        
        #print(structure(c(nguides, i, score_on_best, score_off_best), 
        #                names = c('N_guides', 'Iteration', 'Son_tot', 'Soff_tot')))
       
        # Throw one kmer randomly
        kmers_subs <- sample(kmers_best, length(kmers_best) -1)
        
        # Score kmer subset
        kmers_subs_score <- if (length(kmers_subs) == 1) { mat[kmers_subs,] } else { matrixStats::colMaxs(mat[kmers_subs, ]) }
       
        # Score delta against all kmers
        score_delta <- t(t(mat) - kmers_subs_score)
        score_delta[score_delta < 0] <- 0
        
        # Add best other kmer
        kmers_new <- c(kmers_subs, names(which.max(rowSums(score_delta))))
        
        # Score new subset
        score_all_new <- matrixStats::colMaxs(mat[kmers_new, ])
        score_on_new  <- sum(score_all_new[on_indeces])
        score_off_new <- sum(score_all_new[off_indeces])
            
        # Update kmers if optimized
        if (score_on_new > score_on_best)
        {
          kmers_best <- kmers_new
          score_on_best <- score_on_new
          score_off_best <- score_off_new
        }
        if (score_on_new == score_on_best & score_off_new < score_off_best)
        {
          kmers_best <- kmers_new
          score_off_best <- score_off_new   
        }
        
        # Update results df
        df[i, 'Son_tot']  <- score_on_best
        df[i, 'Soff_tot'] <- score_off_best
        df[i, 'kmer_id'][[1]][[1]]  <- kmers_best
      }
      
      return(df)
    }
      
    #report %>% ggplot(aes(iterations, Son_tot)) + geom_point() + facet_wrap(~n_guides, scales = 'free')
    
    # format report to match combinations for binding
    greedy_res <-
      report %>% 
      unnest %>% 
      group_by(iterations, n_guides) %>% 
        mutate(on_tot = sum(colSums(mat[kmer_id, on_indeces] != 0)!=0),
               off_tot = sum(colSums(mat[kmer_id, off_indeces] != 0)!=0)) %>%
      ungroup %>%
      mutate(kmer_id = as.double(kmer_id),
             enr = ifelse(Soff_tot == 0, (Son_tot) / (Soff_tot + 0.01), Son_tot / Soff_tot),
             best = iterations == max(iterations),
             combi_id = paste0(n_guides, 'V', iterations, '_greedy'))

    # Update combinations
    combinations <- 
      combinations %>% 
      filter(!best | n_guides == 1) %>%
      bind_rows(., greedy_res) %>%
      nest(kmer_id, .key = 'kmer_id')
  }
  guideSet@combinations <- combinations
  return(guideSet)
}
  


