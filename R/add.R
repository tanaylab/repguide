#' Add multiple sequence alignment to guideSet object
#'
#' Calculates rough multiple sequence alignment and consensus sequence or imports from external file
#'
#' @param guideSet guideSet containing targets
#' @param files data.frame. Requires columns 'repname' and 'path' storing family identifiers and paths to the multiple sequence alignment, respectively.
#' @param max_gap_freq Numeric. Removes positions on alignment with higher gap (-) frequency than \code{max_gap_freq}.
#' @param iterations Numeric. Passed to AlignSeqs from DECIPHER package.
#' @param refinements Numeric. Passed to AlignSeqs from DECIPHER package.
#' @param force Logical. If TRUE, overwrites existing alignments.
#' @return Returns a guideSet object with alignments.
#' @export
addAlignments <- function(guideSet,
                          files = NULL, # Either single character vector or df
                          max_gap_freq = 0.8,
                          iterations = 2,
                          refinements = 1,
                          force = FALSE
                          )
{
  if(length(guideSet@alignments) != 0 & !force) { stop('guideSet already contains alignments. Use force = TRUE to overwrite (will remove all downstream results)')}
  if(length(guideSet@targets == 0)) { stop('Add targets to guideSet using addTargets function before calling addAlignments') }
    
  # Remove downstream results
  slot(guideSet, name = 'kmers') <- GRanges()
  slot(guideSet, name = 'combinations') <- tibble()
  slot(guideSet, name = 'plots') <- list('targets' = list(), 'guides' = list(), 'combinations' = list())
  slot(guideSet, name = 'guide_length') <- numeric()
  slot(guideSet, name = 'PAM') <- character()
  slot(guideSet, name = 'alignments') <- DNAStringSetList()
  slot(guideSet, name = 'consensus') <- DNAStringSet() 
  
  n_cores <- guideSet@.n_cores
  seed <- guideSet@.seed
  repnames_all <- guideSet@families
  repnames_comp <- sort(repnames_all[!repnames_all %in% files$repname])
  repnames_import <- sort(repnames_all[repnames_all %in% files$repname])
  targets <- guideSet@targets
  targets_comp <- targets[targets$repname %in% repnames_comp]
  
  # Compute alignments
  message(paste0('Computing alignments for ', paste(repnames_comp, collapse = ', ')))
  alignments_comp <- 
    DNAStringSetList(unlist(tapply(1:length(targets_comp), targets_comp$repname, function(x)
    {
      .compMSA(targets_comp[x], max_gap_freq, iterations, refinements, seed = seed)  
    })))
  
  # Import alignments
  if(!is.null(files)) 
  {
    message(paste0('Importing alignments for ', paste(repnames_import, collapse = ', ')))
    files$path <- as.character(files$path)
    alignments_import <- 
      DNAStringSetList(sapply(files$path, function(x)
      {
        msa <- Biostrings::readDNAStringSet(x)
        msa_wo_gaps <- .rmGaps(msa, max_gap_freq)
        return(msa_wo_gaps)
      }))
    names(alignments_import) <- files$repname
    guideSet@alignments <- c(alignments_comp, alignments_import)  
  } else { 
    guideSet@alignments <- alignments_comp
  }
   
  # Compute consensus seqs
  guideSet@consensus <-   
    unlist(DNAStringSetList(sapply(guideSet@alignments, 
                                   DECIPHER::ConsensusSequence, 
                                      noConsensusChar = 'N', 
                                      threshold = 0.99, 
                                      ignoreNonBases = TRUE, 
                                      minInformation = 0.1
                                  )
                            )
          )

  guideSet@calls$alignments <- match.call()
  return(guideSet)
}

#' Add guideRNA combinations to guideSet object.
#'
#' Selects set of non-redundant guideRNAs and computes all potential combinations against targets of interest.
#' 
#' @param guideSet guideSet object containing guides.
#' @param max_guides Numeric. Maximum number of distinct guides to consider when calculating combinations.
#' @param greedy Logical. If \code{TRUE} (the default), performs an additional greedy optimization for selected guide combinations.
#' @param iterations Integer > 0. Number of greedy search iterations. Usually 10 iterations are sufficient for convergence.
#' @param method String. Method of how to pick the best guideRNA per cluster. 
#' @param coeff Integer. If \code{method} is \code{regression}, \code{coeff} modifies the slope of the linear regression.
#' @param force Logical. If \code{TRUE}, overwrites existing results
#' @return Returns a guideSet object with combinations. 
#' @export 
addCombinations <- function(guideSet, 
                            max_guides = 5,
                            greedy = TRUE,
                            iterations = 10,
                            method = c('max', 'min', 'ratio', 'regression'),
                            coeff = 1,
                            force = FALSE)
{
  if (length(guideSet@kmers) == 0) { stop('Add guides to guideSet using addGuides function before calling addCombinations') }
  if (length(guideSet@combinations) != 0 & !force) { stop('guideSet already contains combinations. Use force = TRUE to overwrite (will remove QC plots)') }
  if (max_guides > length(unique(gs@kmers$kmer_id[gs@kmers$best]))) { stop ('max_guides cannot exceed number of selected guides') }
  
  # Remove downstream results
  slot(guideSet, name = 'plots') <- list('targets' = guideSet@plots$targets, 'guides' = guideSet@plots$guides, 'combinations' = list())
  
  guideSet <- .compCombinations(guideSet, n_guides_max = max_guides)
  guideSet <- .selBestKmers(guideSet, method = method, coeff = coeff)
  if (greedy) { guideSet <- .compGreedy(guideSet, iterations) }
  
  guideSet@calls$combinations <- match.call()
  return(guideSet)
}

#' Add guides to a guideSet object.
#'
#' Adds all potential guideRNAs against the targets of interest. 
#' 
#' @param guideSet guideSet object with targets. 
#' @param guides Character. Optional vector of pre-computed guideRNA sequences to map and annotate. Will restrict downstream analysis to this set of guides.
#' @param n_mismatches Single integer of either 0, 1, 2, or 3. Maximal number of tolerated mismatches when assessing guideRNA binding targets. Defaults to 0.
#' @param guide_length Single integer between 15 and 25. Basepair size of the guideRNAs. Defaults to 19. 
#' @param gc_content Numeric vector. Allowed GC content range of guides, e.g. c(0.4, 0.8) blacklists guides with GC content lower and higher than 40% and 80%, respectively 
#' @param min_Son Numeric. Minimal on target score of guides. Defaults to 10.
#' @param min_Soff Numeric. Maximal off target score of guides. Defaults to 50. 
#' @param n_clust Single positive integer <= 20. Number of groups to cluster guideRNAs into. 
#'   Higher \code{n_clust} usually gives better results but comes with a speed penalty. 
#' @param consensus_range DataFrame with repname, start, and end columns. Discards guideRNAs targeting parts outside of \code{consensus_range} on the consensus.
#' @param PAM Character. Currently only 'NGG' PAM is supported.
#' @param lower_count Numeric. Passed to jellyfish kmer counting. Only kmers occuring at least \code{lower_count} times are considered for mapping.
#' @param force Logical. If \code{TRUE}, overwrite existing guides.
#' @return Returns a guideSet object containing guides.
#' @examples
#' @export
addGuides <- function(guideSet, 
                      guides = NULL,
                      n_mismatches = 0, 
                      guide_length = 19, 
                      gc_content = c(0.4, 0.8),
                      min_Son = 10,
                      max_Soff = 50,
                      consensus_range = NULL,
                      n_clust = 11,
                      method = 'max',
                      coeff = 1,
                      PAM = 'NGG',
                      lower_count = 5,
                      force = FALSE)
{
  if (!is.null(guides) & class(guides) != 'character') { stop ('Provided guides must be a character vector') }
  if (length(unique(nchar(guides))) > 1) { stop ('guides must be of same length') }
  if (!n_mismatches %in% c(0, 1, 2, 3)) { stop('Mismatches must be 0, 1, 2, or 3') }
  if (guide_length < 15 | guide_length > 25 ) { stop('Guide length must be between 15 and 25') }
  if (length(guideSet@targets) == 0) { stop('Add targets to guideSet using addTargets function before calling addGuides') }
  if (!is.null(consensus_range) & length(guideSet@alignments) == 0) { stop ('No consensus model found. Call addAlignments on guideSet or omit consensus range') }
  
  # Check if slot already exists
  if (length(guideSet@kmers) != 0 & !force)
  {
    stop('guideSet already contains guides. Use force = TRUE to overwrite (will remove all downstream results)')
  }
  
  # Remove downstream results
  slot(guideSet, name = 'combinations') <- tibble()
  slot(guideSet, name = 'plots') <- list('targets' = guideSet@plots$targets, 'guides' = list(), 'combinations' = list())
 
  guideSet@guide_length <- ifelse (is.null(guides), guide_length, unique(nchar(guides)))
  guideSet@PAM <- PAM
  
  # Add results to guideSet
  if (is.null(guides)) 
  { 
    guideSet <- .jellyfish(guideSet, lower_count) 
  } else { 
    guideSet@kmers <- GRanges(seqnames = 1:length(guides), ranges = 1:length(guides), seq_guide = guides)
  }
  guideSet <- .bowtie(guideSet, n_mismatches)
  guideSet <- annoGuides(guideSet)
  guideSet <- selGuides(guideSet, min_Son = min_Son, max_Soff = max_Soff, consensus_range = consensus_range, gc_content = gc_content)
  if (is.null(guides) | length(guides) > 20) 
  { 
    guideSet <- clustGuides(guideSet, n_clust = n_clust) 
    guideSet <- .selBestKmers(guideSet, method = method, coeff = coeff)
  } else { 
    guideSet@kmers$kmer_clust <- NA
    guideSet@kmers$te_clust <- NA
    guideSet@kmers$best = ifelse(guideSet@kmers$valid, TRUE, FALSE)
  }
    
  guideSet@calls$guides <- match.call()
  return(guideSet)
}   

#' Add targets to a guideSet object.
#'
#' @param guideSet guideSet object containing genome annotation. 
#' @param targets Character vector or GRanges object. Either identifier of repeat families or GRanges object with target coordinates.
#' @param blacklist Boolean. Should targets overlapping cis regulatory regions be blacklisted. False by default.
#' @param blacklist_dist Integer >= 0. Minimum basepair distance of target loci to cis regulatory regions.
#' @param force. If \code{TRUE}, overwrites existing results.
#' @return Returns a guideSet object containing targets
#' @examples
#' guideSet <- addTargets(guideSet, targets = c('LTR12C', 'LTR12E'))
#' @export
addTargets <- function(
                       guideSet,
                       targets = NULL, # either repnames, GRanges with coords, or seqs
                       force = FALSE
                       )
                       
{
  # Check if slot already exists
  if(length(guideSet@targets) != 0 & !force)
  {
    stop('guideSet already contains targets. Use force = TRUE to overwrite (will remove all downstream results)')
  }
   
  # Remove downstream results
  slot(guideSet, name = 'kmers') <- GRanges()
  slot(guideSet, name = 'combinations') <- tibble()
  slot(guideSet, name = 'plots') <- list('targets' = list(), 'guides' = list(), 'combinations' = list())
  slot(guideSet, name = 'guide_length') <- numeric()
  slot(guideSet, name = 'PAM') <- character()
  slot(guideSet, name = 'alignments') <- DNAStringSetList()
  slot(guideSet, name = 'consensus') <- DNAStringSet()  
  
  # Add targets based on family names
  if (class(targets) != 'GRanges')
  {
    tes <- guideSet@tes
    if (length(setdiff(targets, tes$repname)) > 0) { warning('Target(s) ', setdiff(targets, tes$repname), ' not found in guideSet') }
    targets <- tes[tes$repname %in% targets]  
  }
  
  if (class(targets)!= 'GRanges') { stop ('Targets must be family name(s) or a GRanges object') }
  
  # Remove loci overlapping blacklisted regions from targets
  if (length(guideSet@blacklist) > 0)
  {
    blacklist <- guideSet@blacklist
    start(blacklist) <- start(blacklist)
    end(blacklist) <- end(blacklist)
    targets <- targets[-findOverlaps(targets, blacklist, ignore.strand = TRUE)@from]
  }
  
  # Add locus seq column
  targets$seq <- getSeq(guideSet@genome, targets)
  
  guideSet@targets <- targets
  guideSet@families <- unique(guideSet@targets$repname)
  
  guideSet@calls$targets <- match.call()
  return(guideSet)
}	
