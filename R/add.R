#' Add multiple sequence alignment to guideSet object
#'
#' Either calculates multiple sequence alignment and consensus sequence or imports from external file
#'
#' @param guideSet guideSet with targets
#' @param file String
#'
#'
#'
#' @export
addAlignments <- function(guideSet,
                          files = NULL, # Either single character vector or df
                          max_gap_freq = 0.8,
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
      .compMSA(targets_comp[x], n_cores = n_cores, max_gap_freq)  
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
#' @param guideSet guideSet object with kmers.
#' @param n_clust Single positive integer <= 20. Number of groups to cluster guideRNAs into. 
#'   Higher \code{n_clust} usually gives better results but comes with a speed penalty.
#' @param method String. Method of how to pick the best guideRNA per cluster. 
#' @param coeff Integer. If \code{method} is \code{regression}, \code{coeff} modifies the slope of the linear regression.
#'
#' @return Returns a guideSet object 
#' @examples
#' gs <- createGuideSet(Hsapiens)
#' gs <- addTargets(targets = c('LTR12C', 'LTR12E'), alignment = TRUE)
#' @export 
addCombinations <- function(guideSet, 
                            n_guides = 5,
                            method = c('max', 'min', 'ratio', 'regression'),
                            coeff = 1,
                            force = FALSE)
{
  if(length(guideSet@kmers) == 0) { stop('Add guides to guideSet using addGuides function before calling addCombinations') }
  if(length(guideSet@combinations) != 0 & !force) { stop('guideSet already contains combinations. Use force = TRUE to overwrite (will remove QC plots)') }
    
  # Remove downstream results
  slot(guideSet, name = 'plots') <- list('targets' = guideSet@plots$targets, 'guides' = guideSet@plots$guides, 'combinations' = list())
  
  guideSet <- .compCombinations(guideSet, n_guides_max = n_guides)
  guideSet <- .selBestKmers(guideSet, method = method, coeff = coeff)
  
  guideSet@calls$combinations <- match.call()
  return(guideSet)
}

#' Add guides to a guideSet object.
#'
#' Adds all potential guideRNAs against the targets of interest. 
#' 
#' @param guideSet guideSet object with targets. 
#' @param n_mismatches Single integer of either 0, 1, 2, or 3. Maximal number of tolerated mismatches when assessing guideRNA binding targets. Defaults to 0.
#' @param guide_length Single integer between 10 and 30. Basepair size of the guideRNAs. 
#'   A spacer length of 19 (the default) has been optimal in 'Sequence determinants of improved CRISPR sgRNA design', 
#'     but shorter guides may result in fewer actual off-targets 'Improving CRISPR-Cas nuclease specificity using truncated guide RNAs'
#' @param min_targets Single integer >= 0. Discard guideRNAs with less than \code{min_targets} hits on target.
#' @param max_off_targets Single integer >= 0. Discard guideRNAs with more than \code{max_off_targets} off target hits.
#' @param consensus_range DataFrame with repname, start, and end columns. Discard guideRNAs targeting parts outside of \code{consensus_range} on the consensus.
#' @param PAM Character. Currently only 'NGG' PAMs are supported.
#' 
#' @return Returns a guideSet object 
#' @examples
#' gs <- createGuideSet(Hsapiens)
#' gs <- addTargets(targets = c('LTR12C', 'LTR12E'), alignment = TRUE)
#' @export
addGuides <- function(guideSet, 
                      n_mismatches = 0, 
                      guide_length = 19, 
                      gc_content = c(0.4, 0.8),
                      min_Son = 10,
                      max_Soff = Inf,
                      consensus_range = NULL,
                      n_clust = 15,
                      method = 'max',
                      coeff = 1,
                      PAM = 'NGG',
                      force = FALSE)
{
  if(!n_mismatches %in% c(0, 1, 2, 3)) { stop('Mismatches must be 0, 1, 2, or 3') }
  if(guide_length < 15 | guide_length > 25 ) { stop('Guide length must be between 15 and 25') }
  if(length(guideSet@targets) == 0) { stop('Add targets to guideSet using addTargets function before calling addGuides') }
  
  # Check if slot already exists
  if(length(guideSet@kmers) != 0 & !force)
  {
    stop('guideSet already contains guides. Use force = TRUE to overwrite (will remove all downstream results)')
  }
  
  # Remove downstream results
  slot(guideSet, name = 'combinations') <- tibble()
  slot(guideSet, name = 'plots') <- list('targets' = guideSet@plots$targets, 'guides' = list(), 'combinations' = list())
 
  guideSet@guide_length <- guide_length
  guideSet@PAM <- PAM
  
  guideSet <- .bowtie(guideSet, n_mismatches = n_mismatches)
  guideSet <- annoGuides(guideSet)
  guideSet <- selGuides(guideSet, min_Son = min_Son, max_Soff = max_Soff, consensus_range = consensus_range, gc_content = gc_content)
  guideSet <- clustGuides(guideSet, n_clust = n_clust)
  guideSet <- .selBestKmers(guideSet, method = method, coeff = coeff)
  
  guideSet@calls$guides <- match.call()
  return(guideSet)
}   

#' Add targets to a guideSet object.
#'
#' Adds target loci you want to design guideRNAs against.
#' 
#' @param guideSet guideSet object to add targets to. 
#' @param targets Character or GRanges. Either names of repeat families or GRanges object with target coordinates.
#' @param alignment Boolean. Should a rough multiple sequence alignment and consensus sequence be added? 
#'   Useful for inspecting or targeting specific parts of a family.
#' @param blacklist Boolean. Should targets overlapping cis regulatory regions be blacklisted. False by default.
#' @param blacklist_dist Integer >= 0. Minimum basepair distance of target loci to cis regulatory regions.
#' 
#' @return Returns a guideSet object 
#' @examples
#' gs <- createGuideSet(Hsapiens)
#' gs <- addTargets(targets = c('LTR12C', 'LTR12E'), alignment = TRUE)
#' @export
addTargets <- function(
                       guideSet,
                       targets = NULL, # either repnames, GRanges with coords, or seqs
                       max_gap_freq = 0.8,
                       blacklist = FALSE, 
                       min_dist = 0,
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
  
  if (class(targets) == 'GRanges')
  {
    guideSet@targets <- targets
  } else {
    tes <- guideSet@tes
    guideSet@targets <- tes[tes$repname %in% targets]  
  }
  
  if (blacklist)
  {
    blacklist <- guideSet@cis
    start(blacklist) <- start(blacklist) - min_dist
    end(blacklist) <- end(blacklist) + min_dist    

    #guideSet@blacklisted <- blacklist
    guideSet@targets$blacklisted <- 1:length(guideSet@targets) %in% findOverlaps(guideSet@targets, blacklist)@from
  }
  
  guideSet@targets$seq <- getSeq(guideSet@genome, guideSet@targets)
  guideSet@families <- unique(guideSet@targets$repname)
  
  guideSet@calls$targets <- match.call()
  return(guideSet)
}	
