.checkValidity <- function(object)
{
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
                                 blacklisted = 'GRanges',
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
  function(.Object, 
           alt_chromosomes,
           genome, 
           tes,
           cis,
           refdir,
           n_cores,
           seed
           ) 
    {        
      if (fs::is_dir(refdir))
      {
        message ('Searching for bowtie indeces in ', 'refdir')
      } else {
        refdir <- paste0(system.file(package = 'Repguide'), '/bowtie_indeces')
        dir.create(index_dir, showWarnings = FALSE)
        warning ('No bowtie index directory provided. Will search for indeces in ', index_dir, ' and create them if not found')
      }
         
      if (is.null(n_cores)) { n_cores <- parallel::detectCores() - 5 }
      if (class(tes) != 'GRanges') { if (fs::is_file(tes)) { if (grepl('bed.gz$', tes) | grepl('.bed$', tes)) 
      { 
        tes <- rtracklayer::import.bed(tes) 
      } else {
        tes <- importTEs(tes) 
      }}}

      if (class(cis) != 'GRanges') { if (fs::is_file(cis)) { cis <- rtracklayer::import.bed(cis) }}
      
      if (!alt_chromosomes)
      {
        genome <- .keepBSgenomeSequences(genome, grep('_', seqnames(genome), value = TRUE, invert = TRUE))
      }
      
      .Object@genome <- genome
      .Object@tes <- tes
      .Object@cis <- cis
      .Object@plots <- list('targets' = list(), 'guides' = list(), 'combinations' = list())
      .Object@refdir <- refdir
      .Object@.n_cores <- n_cores  
      .Object@.seed <- seed
      
      if (is.element('doMC', installed.packages())[1]) 
      { 
        doMC::registerDoMC(n_cores) 
      } else {
        warning ('Package doMC not found. Install package to support multicore usage')
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
  message(glue::glue('with {n_cores} registered cores'))
  message(glue::glue('matching bowtie index found: {found_index}'))
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
                           tes = GRanges(), 
                           cis = GRanges(), 
                           n_cores = NULL, 
                           refdir = '', 
                           seed = 19)
{
  new('guideSet', genome, tes, cis, refdir, n_cores, seed)
}                         