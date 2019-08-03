#' Blacklist guideRNAs
#'
#' Adds a blacklist flag to guides not fulfilling defined criteria
#'
#' @param min_Son Numeric. Guides below \code{min_Son} are blacklisted. If \code{NULL} (the default), \code{min_Son} is calculated based on the target score distribution (75th percentile).
#' @param max_Soff Numeric. Guides exceeding \code{min_Son} are blacklisted. If \code{NULL} (the default), \code{max_Soff} is calculated based on the off-target score distribution (99.9th percentile).
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
    kmer_stats <- kmers %>% as_tibble %>% .kmerStats
  
    if (is.null(min_Son))
    {
      min_Son <- round(quantile(kmer_stats$Son_tot, 0.75), 2)
      message ('Set min_Son to ', min_Son)
    }
    if (is.null(max_Soff)) 
    {
      max_Soff <- round(quantile(kmer_stats$Soff_tot, 0.999), 2)
      message ('Set max_Soff to ', max_Soff)
    }
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
  guideSet@kmers$kmer_clust <- NA
  guideSet@kmers$te_clust <- NA
  guideSet@kmers$best <- FALSE
  
  return(guideSet)  
}

.selBestKmers <- function(guideSet,
                          alpha = 1,
                          type = 'kmer')
                          #method = 'deprecated')
                          #coeff = 1)
{
  if (alpha == Inf) { stop ('Alpha must be a finite number') }
  if (alpha <0) { stop ('Alpha coefficient must be >= 0') }
  
  if (type == 'kmer')
  {
    kmers <- as_tibble(guideSet@kmers) %>% filter(valid)
    
    # Calc Kmer stats
    stats <- 
      .kmerStats(kmers) %>%
      left_join(., kmers %>% select(kmer_id, kmer_clust) %>% distinct, by = 'kmer_id') %>%
      dplyr::rename(grouping_var = kmer_clust, id = kmer_id) %>%
      mutate(enr = ifelse(Soff_tot == 0, (Son_tot) / (Soff_tot + 0.01), Son_tot / Soff_tot)) # adding a pseudo-count to avoid Inf
    # ================
    
    ids_best <-
        stats %>%
        group_by(grouping_var) %>%
        mutate(alpha_score = Son_tot - Soff_tot * alpha,
               best = alpha_score == max(alpha_score)) %>%
        ungroup %>%
        filter(best) %>%
        pull(id)
        
    guideSet@kmers$best <- guideSet@kmers$kmer_id %in% ids_best
  }
  
  if (type == 'combination')
  {
    stats <- guideSet@combinations %>% dplyr::rename(grouping_var = n_guides, id = combi_id)
    ids_best <-
        stats %>%
        group_by(grouping_var) %>%
        mutate(alpha_score = Son_tot - Soff_tot * alpha,
               best = id == id[which.max(alpha_score)]) %>%
        ungroup %>%
        filter(best) %>%
        pull(id)
    
    guideSet@combinations$best <- guideSet@combinations$combi_id %in% ids_best
  }
  
  return(guideSet)  
}  

# Shrinks target regions overlapping blacklisted domains
# Currently doesn't deal with blacklisted regions in center of repeat
.selBlacklistRegions <- function(targets, blacklist)
{
    orig_size <- width(targets)
    targets_bp_res <- unlist(tile(targets, width = 1))
    targets_bp_res$te_id <- rep(targets$te_id, orig_size)
    
    targets_bp_res_filt <- subsetByOverlaps(targets_bp_res, blacklist, ignore.strand = TRUE, invert = TRUE)
    
    targets_resized  <- unlist(GenomicRanges::reduce(split(targets_bp_res_filt, ~te_id)))
    targets_resized$te_id <- as.numeric(names(targets_resized))
    
    targets_resized_anno <- left_join(as_tibble(targets_resized), 
                                      as_tibble(targets) %>% select(-seqnames, - start, -end, -width, -strand), 
                                      by = 'te_id') %>%
                            makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    return(targets_resized_anno)
}
  

                         