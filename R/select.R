#' @export
selGuides <- function(guideSet,
                      min_Son = 10,
                      max_Soff = Inf,
                      consensus_range = NULL,
                      gc_content = c(0.4, 0.8)
                      )
{
  print('Blacklisting kmers')
  kmers <- as_tibble(guideSet@kmers)
    
  if (!is.null(consensus_range))
  {
    if (sum(colnames(consensus_range) %in% c('repname', 'start', 'end')) != 3) { stop ('Consensus range requires repname, start, end columns') }
    
    kmers <- 
      right_join(consensus_range, 
                 kmers, 
                 by = 'repname', suffix = c('_con', '')) %>%
      mutate(outside = con_pos < start_con | con_pos > end_con,
             on_target = ifelse(!is.na(outside) & outside, 0, on_target)) %>%
      select(-start_con, -end_con, -outside)
  }
   
  kmers_gc_valid <-
    kmers %>%
    filter(between(gc, gc_content[1], gc_content[2])) %>%
    pull(kmer_id) %>%
    unique

  kmers_score_valid <-
    kmers %>%
    .kmerStats %>%
    filter(Son_tot >= min_Son & Soff_tot <= max_Soff) %>%
    pull(kmer_id) %>%
    unique
    
  guideSet@kmers$on_target <- kmers$on_target
  guideSet@kmers$valid_gc <- guideSet@kmers$kmer_id %in% kmers_gc_valid
  guideSet@kmers$valid_score <- guideSet@kmers$kmer_id %in% kmers_score_valid
  guideSet@kmers$valid <- guideSet@kmers$valid_gc + guideSet@kmers$valid_score == 2

  return(guideSet)  
}

.selBestKmers <- function(guideSet,
                          method = 'regression',
                          coeff = 1)
{
  kmers <- as_tibble(guideSet@kmers) %>% filter(valid)

  n_iter <- c('Kmer' = !is.null(guideSet@kmers$kmer_clust), 'Combi' = length(guideSet@combinations) > 0)
  
  if (sum(n_iter) == 0) { stop ('No Kmer clustering or combinations found') }
  if (n_iter['Kmer'])
  {
    # Calc Kmer stats
    kmer_stats <- 
      .kmerStats(kmers) %>%
      left_join(., kmers %>% select(kmer_id, kmer_clust) %>% distinct) %>%
      rename(grouping_var = kmer_clust, id = kmer_id) %>%
      mutate(enr = ifelse(Soff_tot == 0, (Son_tot) / (Soff_tot + 0.01), Son_tot / Soff_tot)) # adding a pseudo-count to avoid Inf
    # ================
  }
  
  for (iter in names(n_iter)[n_iter]) 
  {
    if (iter == 'Kmer') { stats <- kmer_stats }
    if (iter == 'Combi') { stats <- guideSet@combinations %>% rename(grouping_var = n_guides, id = combi_id) }
    # Mark best kmer per clust based on method
    if (method == 'max')
    {
      ids_best <- 
        stats %>%
        group_by(grouping_var) %>%
        top_n(1, Son_tot) %>%
        top_n(-1, Soff_tot) %>%
        sample_n(1) %>%
        ungroup %>%
        pull(id)
    }
    
    if (method == 'min')
    {
      ids_best <- 
        stats %>%
        group_by(grouping_var) %>%
        top_n(-1, Soff_tot) %>%
        top_n(1, Son_tot) %>%
        sample_n(1) %>%
        ungroup %>%
        pull(id)
    }
    
    if (method == 'ratio')
    {
      ids_best <- 
        stats %>%
        group_by(grouping_var) %>%
        top_n(1, enr) %>%
        top_n(1, Son_tot) %>%
        sample_n(1) %>%
        ungroup %>%
        pull(id)  
    }
    
    if (method == 'regression')
    {
      # Calc linear regression per kmer clust
      clust_lm <- 
        stats %>% 
        group_by(grouping_var) %>% 
        do(clust_fit = lm(Son_tot ~ Soff_tot, data = .)) %>%
        broom::tidy(., clust_fit) %>%
        ungroup %>%
        select(grouping_var, term, estimate) %>%
        spread(term, estimate) %>%
        rename(intercept = `(Intercept)`, slope = Soff_tot) %>%
        mutate(slope = slope * coeff)
        
      ids_best <-
        left_join(stats, clust_lm) %>%
        mutate(delta = Son_tot - (intercept + Soff_tot * slope)) %>%
        group_by(grouping_var) %>%
        top_n(1, delta) %>%
        top_n(1, Son_tot) %>%
        sample_n(1) %>%
        ungroup %>%
        pull(id)
    }
    if (iter == 'Kmer') { guideSet@kmers$best <- guideSet@kmers$kmer_id %in% ids_best }
    if (iter == 'Combi') { guideSet@combinations$best <- guideSet@combinations$combi_id %in% ids_best }
  }
  # =====================================
  
  return(guideSet)  
}  
  

                         