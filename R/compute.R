.compGuideScores <- function(guideSet) # Sbind can be interprated as the likelyhood of binding to a site (0 - 100)
{
  print('Computing guide scores')
  kmers <- guideSet@kmers
  kmers$Sbind <- 1
  PAM <- guideSet@PAM
  guide_length <- guideSet@guide_length + nchar(PAM)
    
  mismatches <- stringr::str_extract_all(kmers$mismatches, '[0-9]+', simplify = TRUE)
  rownames(mismatches) <- 1:nrow(mismatches)

  mismatches_long <- 
    mismatches %>% 
    as.data.frame %>% 
    rownames_to_column %>% 
    gather(V1, V2, -rowname) %>% 
    as_tibble %>% 
    select(-V1) %>% 
    filter(V2 != '') %>%
    mutate(V2 = as.numeric(V2)) %>%
    data.table::setDT(key = 'rowname')
    
  # Calculate Score
  # Formula adapted from Breaking-Cas: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987939/
  .scoring <- function(weight, d, m, guide_length)
  {
    s <- weight * (1 / (((guide_length - d) / guide_length) * 1 + 1)) * (1 / (m * m))
    return(s)
  }
  
  mismatches_stats <- 
    mismatches_long[ , ({m = .N;                                # total mismatches
                         d = (max(V2) - min(V2)) / (m - 1);     # approximate of average pairwise distance between mismatches
                         weight = prod((1 - V2 / guide_length));            # custom weighting function based on position in guide
                         list(m = m, d = d, weight = weight, pos = V2)}),
                         by = rowname]
  mismatches_stats[is.na(d), d := 0]
  mismatches_stats <- subset(unique(mismatches_stats, by = c('rowname')), select = -pos)
  mismatches_stats[, Sbind := .scoring(weight, d, m, guide_length)]
  ######################################################################################

  kmers[as.numeric(mismatches_stats$rowname), ]$Sbind <- mismatches_stats$Sbind                      

  # Compute cis score (Boltzmann sigmoid)
  .scoreCis <- function(x, slope = 5, midpoint = 2000)
  {
    s <- 1 + ((0.01 - 1) / (1 + exp((midpoint - x)/(slope * 100)))) # Boltzmann
    #s <- pmin(100, s)
    return(s)  
  } 
  kmers$Scis = .scoreCis(kmers$cis_dist)
  kmers[is.na(kmers$Scis)]$Scis <- 1
  #######################################
  
  # Compute Soff
  kmers$Soff <- kmers$Sbind * kmers$Scis
  kmers[kmers$on_target >= 0]$Soff <- 0
  #kmers$Soff <- round(kmers$Soff, 2)
  
  # Compute Son
  kmers$Son <- kmers$Sbind
  kmers[kmers$on_target <= 0]$Son <- 0
  #kmers$Son <- round(kmers$Son, 2)
  
  guideSet@kmers <- kmers
  return(guideSet)
}

# Calculates given Kmer input the coverage along given MSA and outputs its coverage in long format
.compCov <- function(msa,
                     kmers,
                     GUIDE_LENGTH = 19)
{
  if(class(msa)[1] == 'DNAStringSet')
  {
    print('Converting MSA to long format')
    msa <- msaToLong(msa)
    print('done')
  }
  
  n_loci = unique(length(msa$te_id))
  
  loci_kmers_cov <- 
    kmersToCov(kmers, guide_length = GUIDE_LENGTH) %>%
    dplyr::rename(pos_wo_gaps = pos)

  msa_kmers_cov <-
    left_join(msa, loci_kmers_cov, by = c('te_id', 'pos_wo_gaps')) %>%
    replace_na(list(n = 0)) # add 0 coverage to non-covered pos
    
  return(msa_kmers_cov)
}

.compCombinations <- function(guideSet,
                              n_guides_max = 10,
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
                                              Soff = max(Soff),
                                              on_tot = sum(Son > 0),
                                              off_tot = sum(Soff > 0)),
                                              by = c('combi_id', 'unique_id')
                                         ][, .(Son_tot = sum(Son), 
                                               Soff_tot = sum(Soff)), 
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

.compMSA <- function(seqs,
                     max_gap_freq = 0.8,
                     kmer_length = 7,
                     n_clust = 10,
                     n_cores = 1,
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
                 clust = cutree(fastcluster::hclust(tgs_dist(kmer::kcount(seqs_bin, k = kmer_length)), 'ward.D2'), n_clust))
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

  print(paste0('Aligning ', length(ids_sel), ' sequences'))
  
  alignment <- DECIPHER::AlignSeqs(seqs[ids_sel], 
                                   processors = n_cores, 
                                   iterations = 1, 
                                   refinements = 0, 
                                   useStructures = FALSE,
                                   verbose = FALSE)
  alignment_wo_gaps <- .rmGaps(alignment, max_gap_freq = max_gap_freq)
  
  return(alignment_wo_gaps)
}

# .compConsensus <- function(msa)
# {
  # if (class(msa) == 'DNAStringSet')
  # {
    # msa <- msaToLong(msa)
  # }
  
  # con_seq <-
    # msa %>% 
    # filter(base != '-') %>% 
    # count(pos, base) %>% 
    # group_by(pos) %>% 
      # top_n(1) %>% 
      # sample_n(1) %>% 
    # ungroup %>% 
    # pull(base) %>% 
    # glue::glue_collapse(.) %>% 
    # DNAStringSet
    
  # return(con_seq)
# }






