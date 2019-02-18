#' Blacklist guideRNAs
#'
#' Adds a blacklist flag to guides not fulfilling defined criteria
#'
#' @param min_Son Numeric. Guides below \code{min_Son} are blacklisted. If \code{NULL} (the default), \code{min_Son} is calculated based on the target score distribution (95th percentile).
#' @param max_Soff Numeric. Guides exceeding \code{min_Son} are blacklisted. If \code{NULL} (the default), \code{max_Soff} is calculated based on the off-target score distribution (99th percentile).
#' @param consensus_range Data.frame with repname, start, and end columns. Scores guideRNA target binding sites outside of provided consensus range neutrally.
#' @param gc_content Numeric vector of length two. Elements must be from 0 through 1. For example, c(0.4, 0.8) blacklists guides with GC content not between 40 and 80 percent
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13')
#' gs <- addGuides(gs, guide_length = 16, n_mismatches = 0, gc_content = c(0.25, 0.9), n_clust = 12)
#' gs <- selGuides(gs, min_Son = 100, gc_content = c(0.1, 0.3))
#' gs <- plotGuides(gs)
#' }
#' @export
selGuides <- function(guideSet,
                      min_Son = NULL,
                      max_Soff = NULL,
                      consensus_range = NULL,
                      gc_content = c(0.4, 0.8)
                      )
{
  print('Blacklisting kmers')
  kmers <- as_tibble(guideSet@kmers)
  
  # Compute Son Soff scores if not given  
  if (sum(is.null(min_Son), is.null(max_Soff)) > 0)
  {
    message ('Computing score blacklisting threshold')
    kmer_stats <- guideSet@kmers %>% as_tibble %>% .kmerStats
  
    if (is.null(min_Son))
    {
      min_Son <- round(quantile(kmer_stats$Son_tot, 0.95), 2)
      message ('Set min_Son to ', min_Son)
    }
    if (is.null(max_Soff)) 
    {
      max_Soff <- round(quantile(kmer_stats$Soff_tot, 0.99), 2)
      message ('Set max_Soff to ', max_Soff)
    }
  }

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
      left_join(., kmers %>% select(kmer_id, kmer_clust) %>% distinct, by = 'kmer_id') %>%
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
  

                         