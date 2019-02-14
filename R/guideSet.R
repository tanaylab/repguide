.checkValidity <- function(object)
{
  if (length(subsetByOverlaps(object@whitelist, object@blacklist)) > 0) { stop ('Blacklisted and whitelisted intervals overlap') }
   
  return(TRUE)
} 

#' An S4 class to store all results
#'
#'
#'
#' @exportClass guideSet
guideSet <- 
  setClass('guideSet', 
           slots = c(
                                 genome = 'BSgenome',
                                 tes = 'GRanges',
                                 cis = 'GRanges',
                                 refdir = 'character',
                                 targets = 'GRanges',
                                 blacklist = 'GRanges',
                                 whitelist = 'GRanges',
                                 plots = 'list', 
                                 guide_length = 'numeric', 
                                 families = 'character',
                                 PAM = 'character',
                                 kmers = 'GRanges',
                                 alignments = 'DNAStringSetList',
                                 consensus = 'DNAStringSet',
                                 combinations = 'tbl_df',
                                 calls = 'list',
                                 .n_cores = 'numeric',
                                 .seed = 'numeric'
                                ),
           validity = .checkValidity)

setMethod(
  "initialize",
  signature = "guideSet",
  function(
           .Object, 
           genome,
           alt_chromosomes,
           tes,
           cis,
           refdir,
           blacklist,
           whitelist,
           n_cores,
           seed
           ) 
    {        
      # Filter alternative chromosomes from BSgenome
      if (!alt_chromosomes)
      {
        genome <- .keepBSgenomeSequences(genome, grep('_', seqnames(genome), value = TRUE, invert = TRUE))
      }
      
      # Check Bowtie reference directory
      if (fs::is_dir(refdir))
      {
        message ('Searching for bowtie indeces in ', 'refdir')
      } else {
        refdir <- paste0(system.file(package = 'Repguide'), '/bowtie_indeces')
        dir.create(refdir, showWarnings = FALSE)
        warning ('No bowtie index directory provided. Will search for indeces in ', refdir, ' and create them if not found')
      }
         
      # Detect number of cores
      if (is.null(n_cores)) { n_cores <- parallel::detectCores() - 5 }
      
      # Import 
      blacklist <- suppressWarnings(.importOrNot(blacklist, genome))
      whitelist <- suppressWarnings(.importOrNot(whitelist, genome))
      cis       <- suppressWarnings(.importOrNot(cis, genome))
      tes       <- suppressWarnings(.importOrNot(tes, genome))
      
      .Object@genome <- genome
      .Object@tes <- tes
      .Object@cis <- cis
      .Object@blacklist <- blacklist
      .Object@whitelist <- whitelist
      .Object@plots <- list('targets' = list(), 'guides' = list(), 'combinations' = list())
      .Object@refdir <- refdir
      .Object@.n_cores <- n_cores  
      .Object@.seed <- seed
      
      if (is.element('doMC', installed.packages())[1]) 
      { 
        doMC::registerDoMC(n_cores) 
      } else {
        warning ('Package doMC not found. Install package to support multicore usage')
        n_cores <- 1
      }
      
      timestamp(prefix = 'Created new guideSet at ', suffix = '', quiet = FALSE)
      return(.Object)
    }
)

setMethod("show", "guideSet", function(object) 
{ 
  genome_id <- object@genome@pkgname
  n_targets <- length(object@targets)
  n_guides <- length(unique(object@kmers$kmer_id))
  n_combinations <- ifelse(length(object@combinations) != 0, length(unique(object@combinations$combi_id)), 0)
  n_plots <- sum(sapply(object@plots, length))
  n_cores <- object@.n_cores
  found_index <- sum(grepl(genome_id, dir(object@refdir, pattern = '.ebwt'))) > 0
  
  message(glue::glue('guideSet object of {genome_id}'))
  message(glue::glue('with {n_targets} targets'))
  message(glue::glue('with {n_guides} guides'))
  message(glue::glue('with {n_combinations} combinations'))
  message(glue::glue('with {n_plots} QC plots'))
  message(glue::glue('registered {n_cores} cores'))
  message(glue::glue('matching bowtie index found: {found_index}'))
})

#' @export
setGeneric("mappings", function(object) standardGeneric("mappings"))
setMethod("mappings", signature("guideSet"), function(object) {
  out <- object@kmers
  return(out)
})

#' Exports results of guideSet
#' 
#' @param outdir String. Path to output directory
#'
#'
#' @export
setGeneric('export', function(guideSet, ...) standardGeneric('export'), signature = 'guideSet') 

#' Creates a new guideSet object
#'
#'
#'
#'
#' @export
createGuideSet <- function(genome, 
                           alt_chromosomes = FALSE,
                           tes = NULL, 
                           cis = NULL, 
                           blacklist = NULL,
                           whitelist = NULL,
                           n_cores = NULL, 
                           refdir = '', 
                           seed = 19)
{
  new('guideSet', genome, alt_chromosomes, tes, cis, refdir, blacklist, whitelist, n_cores, seed)
}                         