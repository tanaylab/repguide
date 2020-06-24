#' Add multiple sequence alignment to guideSet object
#'
#' Calculates rough multiple sequence alignment and consensus sequence or imports from external fasta file
#'
#' @param guideSet. guideSet containing targets.
#' @param files data.frame. Imports multiple sequence alignments from file for provided families. For other families, the alignment is computed. Requires columns 'repname' and 'path' storing family identifiers and paths to the multiple sequence alignment, respectively.
#' @param max_gap_freq Numeric between 0 and 1. Removes positions on alignment with higher gap (-) frequency than \code{max_gap_freq}.
#' @param iterations Integer. Passed to AlignSeqs from DECIPHER package.
#' @param refinements Integer. Passed to AlignSeqs from DECIPHER package.
#' @param force Logical. If \code{TRUE}, overwrites existing alignments.
#' @return Returns a guideSet object with multiple sequence alignment and consensus model.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(genome = Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = c('THE1B', 'THE1C'))
#' # import multiple sequence alignment for THE1C and compute for THE1B
#' gs <- addAlignments(gs, files = data.frame('repname' = 'THE1C', path = 'path to multiple sequence alignment fasta file'))
#' gs@alignments
#' }
#' @seealso [DECIPHER::AlignSeqs()]
#' @export
addAlignments <- function(guideSet,
                          files = NULL, # Either single character vector or df
                          max_gap_freq = 0.8,
                          iterations = 2,
                          refinements = 1,
                          force = FALSE
                          )
{
  if (max_gap_freq < 0 | max_gap_freq > 1) { stop ('Max gap frequency must be between 0 and 1') }
  if(length(guideSet@alignments) != 0 & !force) { stop('guideSet already contains alignments. Use force = TRUE to overwrite (will remove all downstream results)')}
  if(length(guideSet@targets) == 0) { stop('Add targets to guideSet using addTargets function before calling addAlignments') }
    
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
  guideSet@alignments <- guideSet@alignments[order(names(guideSet@alignments))]
     
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
#' Finds the optimal combination of non-redundant guideRNAs.
#' 
#' @param guideSet. guideSet object containing guides.
#' @param max_guides Numeric. Maximum number of distinct guides to consider when calculating combinations. Do not use a higher number than experimentally feasible.
#' @param greedy Logical. If \code{TRUE} (the default), performs an additional greedy search optimization for selected guide combinations.
#' @param iterations Integer > 0. Number of greedy search iterations (10 by default).
#' @param alpha Numeric. Off-target score coefficient. Large \code{alpha} penalizes combinations with high off-target score while \code{alpha = 0} ignores off-targets and picks combinations with highest on-target binding. Defaults to 10.
#' @param method String. Method of how to pick the best guideRNA per cluster. 
#' @param coeff Integer. If \code{method} is \code{regression}, \code{coeff} modifies the slope of the linear regression.
#' @param force Logical. If \code{TRUE}, overwrites existing combinations.
#' @return Returns a guideSet object contaning combinations. 
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13A')
#' gs <- addGuides(gs, guide_length = 16)
#' gs <- addCombinations(gs, max_guides = 8)
#' gs <- plotCombinations(gs)
#' }
#' @seealso [plotCombinations()], [createGuideSet()], [addTargets()], and [addGuides()]
#' @export 
addCombinations <- function(guideSet, 
                            max_guides = 5,
                            greedy = TRUE,
                            iterations = 10,
                            alpha = 100,
                            force = FALSE)
{
  if (length(guideSet@kmers) == 0) { stop('Add guides to guideSet using addGuides function before calling addCombinations') }
  if (length(guideSet@combinations) != 0 & !force) { stop('guideSet already contains combinations. Use force = TRUE to overwrite (will remove QC plots)') }
  if (max_guides > length(unique(guideSet@kmers$kmer_id[guideSet@kmers$best]))) { stop ('max_guides cannot exceed number of selected guides') }
  
  # Remove downstream results
  slot(guideSet, name = 'combinations') <- tibble()
  slot(guideSet, name = 'plots') <- list('targets' = guideSet@plots$targets, 'guides' = guideSet@plots$guides, 'combinations' = list())
  
  guideSet <- .compCombinations(guideSet, n_guides_max = max_guides)
  guideSet <- .selBestKmers(guideSet, alpha = alpha, type = 'combination')
  if (greedy) { guideSet <- .compGreedy(guideSet, alpha, iterations) }
  
  guideSet@calls$combinations <- match.call()
  return(guideSet)
}

#' Add guideRNAs to guideSet object
#'
#' Computes target guideRNA universe and maps, annotates, and scores their genomic targets.
#' 
#' @param guideSet. guideSet object containing targets. 
#' @param guides Character. Optional vector of pre-computed guideRNA sequences to map and annotate. Will restrict downstream analysis to \code{guides} rather than computing all possible guides de novo.
#' @param n_mismatches Integer from 0 through 3. Maximal number of tolerated mismatches when assessing guideRNA binding targets. Defaults to 0.
#' @param blacklist_penalty. Numeric. Off-target score multiplier (10 by default) for blacklisted genomic regions. Ignored if \code{blacklisted(guideSet)} is empty.
#' @param guide_length Single integer between 12 and 26. Basepair size of the guideRNAs. Defaults to 16. 
#' @param gc_content Numeric vector of length two. Elements must be from 0 through 1. For example, c(0.4, 0.8) blacklists guides with GC content not between 40 and 80 percent. Passed to \code{selGuides()}.  
#' @param min_Son Numeric. Minimal on target score of guides. If \code{NULL} (the default), will estimate a value based on the on score distribution. Passed to \code{selGuides()}.
#' @param min_Soff Numeric. Maximal off target score of guides. If \code{NULL} (the default), will estimate a value based on the off score distribution. Passed to \code{selGuides()}.
#' @param n_clust Single positive integer <= 20. Number of groups to cluster guideRNAs into. Higher n_clust usually gives better results but comes with a speed penalty when computing guide combinations. Passed to \code{clustGuides()}.
#' @param consensus_range Data.frame with repname, start, and end columns. Scores guideRNA target binding sites outside of provided consensus range neutrally. Passed to \code{selGuides()} 
#' @param alpha Numeric. Off-target score coefficient. Large \code{alpha} penalizes guides with high off-target score while \code{alpha = 0} ignores off-targets and picks guides with highest on-target binding.
#' @param five_prime_seq Character. Sequence requirement for 5' start of guideRNAs, e.g. G nucleotide for transcription from U6 promoter.
#' @param PAM Character. Currently only 'NGG' PAM is supported.
#' @param lower_count Numeric. Passed to jellyfish kmer counting. Only kmers occuring at least \code{lower_count} times are considered for bowtie mapping.
#' @param force Logical. If \code{TRUE}, overwrite existing guides.
#' @return Returns a guideSet object containing guides.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13')
#' gs <- addGuides(gs, guide_length = 16, n_mismatches = 0, gc_content = c(0.25, 0.9), n_clust = 12)
#' gs <- plotGuides(gs)
#' }
#' @seealso [plotGuides()], [createGuideSet()], [addTargets()], [addCombinations()], [clustGuides()], and [selGuides]
#' @export
addGuides <- function(guideSet, 
                      guides = NULL,
                      n_mismatches = 0, 
                      blacklist_penalty = 10,
                      guide_length = 19, 
                      gc_content = c(0.4, 0.8),
                      min_Son = NULL,
                      max_Soff = NULL,
                      consensus_range = NULL,
                      alpha = 100,
                      n_clust = 11,
                      five_prime_seq = NULL,
                      PAM = 'NGG',
                      lower_count = 5,
                      force = FALSE)
{
  if (PAM != 'NGG') { stop ('Currently only NGG PAM is supported') }
  if (!is.null(guides) & class(guides) != 'character') { stop ('Provided guides must be a character vector') }
  if (length(unique(nchar(guides))) > 1) { stop ('guides must be of same length') }
  if (!n_mismatches %in% c(0, 1, 2, 3)) { stop('Mismatches must be 0, 1, 2, or 3') }
  if (guide_length < 12 | guide_length > 26 ) { stop('Guide length must be between 15 and 25') }
  if (length(guideSet@targets) == 0) { stop('Add targets to guideSet using addTargets function before calling addGuides') }
  if (!is.null(consensus_range) & length(guideSet@alignments) == 0) { stop ('No consensus model found. Call addAlignments on guideSet or omit consensus range') }
  
  # Check if slot already exists
  if (length(guideSet@kmers) != 0 & !force)
  {
    stop('guideSet already contains guides. Use force = TRUE to overwrite (will remove all downstream results)')
  }
  
  # Remove downstream results
  slot(guideSet, name = 'combinations') <- tibble()
  slot(guideSet, name = 'plots')        <- list('targets' = guideSet@plots$targets, 'guides' = list(), 'combinations' = list())
 
  guideSet@guide_length <- ifelse (is.null(guides), guide_length, unique(nchar(guides)))
  guideSet@PAM          <- PAM
  
  # Add results to guideSet
  if (is.null(guides)) 
  { 
    guideSet <- .jellyfish(guideSet, lower_count, five_prime_seq) 
  } else { 
    if (PAM == 'NGG')
    {
      guides <- unlist(lapply(guides, function(x) paste0(x, c('AGG'))))
    } else {
      warnings ('No valid PAM provided, (off)targets may be incorrect')
    }
    
    guideSet@kmers <- 
      GRanges(seqnames = 1:length(guides), 
              ranges = 1:length(guides), 
              guide_seq = guides)
  }
  guideSet <- .bowtie(guideSet, n_mismatches)
  guideSet <- .annoGuides(guideSet, blacklist_penalty, consensus_range = consensus_range)
  guideSet <- selGuides(guideSet, min_Son = min_Son, max_Soff = max_Soff, gc_content = gc_content)
  if (is.null(guides) | length(guides) > 20) 
  { 
    guideSet <- clustGuides(guideSet, n_clust = n_clust, alpha = alpha) 
  } else { 
    guideSet@kmers$kmer_clust <- NA
    guideSet@kmers$te_clust <- NA
    guideSet@kmers$best = ifelse(guideSet@kmers$valid, TRUE, FALSE)
  }
    
  guideSet@calls$guides <- match.call()
  return(guideSet)
}   

#' Add targets to a guideSet object
#'
#' @param guideSet guideSet object containing genome annotation. 
#' @param targets Character vector of family identifiers or GRanges object with coordinates and 'repname' metacolumn.
#' @param force Logical. If \code{TRUE}, overwrites existing results.
#' @return Returns a guideSet object containing targets.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_anno)
#' 
#' # Use family identifiers as targets
#' gs <- addTargets(guideSet, targets = c('LTR12C', 'LTR12E'))
#'
#' # Use custom regions as targets
#' target_coords <- GenomicRanges::makeGRangesFromDataFrame(data.frame('chrom'   = c('chr1', 'chr4'), 
#'                                                                     'start'   = c(1000, 10020),
#'                                                                     'end'     = c(1122, 10090),
#'                                                                     'strand'  = c('+', '-'),
#'                                                                     'repname' = c('ImaginativeTE', 'ImaginativeTE'),
#'                                                                     'te_id'  = c(1:2)),
#'                                                          keep.extra.columns = TRUE)
#' gs <- addTargets(guideSet, targets = target_coords)
#' }
#' @seealso [plotTargets()], and [GenomicRanges::makeGRangesFromDataFrame()] 
#' @export
addTargets <- function(
                       guideSet,
                       targets = NULL, # either repnames, GRanges with coords, or seqs
                       force = FALSE
                       )
                       
{
  # Check if slot already exists
  if (length(guideSet@targets) != 0 & !force) { stop('guideSet already contains targets. Use force = TRUE to overwrite (will remove all downstream results)') }
  if (is.null(targets)) { stop ('Please provide targets either as character vector of repname identifiers or GRanges object with coordinates and repname id') }
  if (class(targets) == 'character' & sum(!targets %in% guideSet@tes$repname) > 0) 
  {
    families <- unique(guideSet@tes$repname)
    targets_not_found <- setdiff(targets, families)
    if (length(targets_not_found) > 1)
    {
      suggestions <- families[apply(stringdist::stringdistmatrix(targets_not_found, families), 1, which.min)]
      suggestions <- paste(suggestions, collapse = ', ')
      targets_not_found <- paste(targets_not_found, collapse = ', ')
    } else {
      fam_dists <- structure(stringdist::stringdistmatrix(targets_not_found, families), names = families)
      suggestions <- paste(names(head(sort(fam_dists), 5)), collapse = ', ')
    }
    stop (glue::glue("{targets_not_found} not found in guideSet annotation. Did you mean any of {suggestions}?"))
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
    targets <- tes[tes$repname %in% targets]  
  }
  
  if (class(targets)!= 'GRanges') { stop ('Targets must be family name(s) or a GRanges object') }
  if (is.null(targets$repname) | is.null(targets$te_id)) { stop ('GRanges must contain repname and te_id metadata columns') }
  
  # Remove target domains overlapping blacklisted regions
  if (length(guideSet@blacklist) > 0)
  {
    targets <- .selBlacklistRegions(targets, guideSet@blacklist)
    if (length(targets) == 0) { stop ('No targets passed blacklisting, think of relaxing blacklisted regions') }
  }
  
  # Add locus seq column
  targets$seq <- getSeq(guideSet@genome, targets)
  
  guideSet@targets <- targets
  guideSet@families <- sort(unique(guideSet@targets$repname))
  
  if (length(guideSet@families) > 8) { warning ('Large number of families provided as targets. Think of reducing number!') }
  guideSet@calls$targets <- match.call()
  return(guideSet)
}	
