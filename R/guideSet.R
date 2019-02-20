#' An S4 class to store genomic annotations and results of Repguide functions.
#'
#' @slot genome BSgenome object.
#' @slot tes GRanges object
#' @slot cis GRanges object
#' @slot refdir Character
#' @slot tempdir Character
#' @slot targets GRanges object
#' @slot blacklist GRanges object
#' @slot whitelist GRanges object
#' @slot plots List
#' @slot guide_length Integer
#' @slot families Character
#' @slot PAM Character
#' @slot kmers GRanges
#' @slot alignments DNAStringSetList
#' @slot consensus DNAStringSet
#' @slot combinations Data.frame
#' @slot calls List
#' @slot .n_cores Integer
#' @slot .seed Integer
#' @exportClass guideSet
guideSet <- 
  setClass('guideSet', 
           slots = c(
                                 genome = 'BSgenome',
                                 tes = 'GRanges',
                                 cis = 'GRanges',
                                 refdir = 'character',
                                 temp = 'character',
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
                                )
          )

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
           temp,
           blacklist,
           whitelist,
           n_cores,
           seed
           ) 
    {        
      if (class(genome) != 'BSgenome') { stop ('Provide a BSgenome object as genome') }
      if (is.null(tes)) { stop ('Provide transposable element annotation of genome') }
      # Filter alternative chromosomes from BSgenome
      if (!alt_chromosomes)
      {
        genome <- .keepBSgenomeSequences(genome, grep('_', seqnames(genome), value = TRUE, invert = TRUE))
      }
      if (is.null(temp)) { temp <- tempdir() }
      
      # Check Bowtie reference directory
      if (fs::is_dir(refdir))
      {
        message ('Searching for bowtie indeces in ', refdir)
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
      
      if (is.null(tes$repname)) { stop ('No repname column found in TE annotation') }
      
      .Object@genome <- genome
      .Object@tes <- tes
      .Object@cis <- cis
      .Object@blacklist <- blacklist
      .Object@whitelist <- whitelist
      .Object@plots <- list('targets' = list(), 'guides' = list(), 'combinations' = list())
      .Object@refdir <- refdir
      .Object@temp <- temp
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
  n_tes <- length(object@tes)
  n_fams <- length(unique(object@tes$repname))
  n_targets <- length(object@targets)
  n_blacklisted <- length(object@blacklist)
  n_guides <- length(unique(object@kmers$kmer_id[object@kmers$valid]))
  n_combinations <- ifelse(length(object@combinations) != 0, length(unique(object@combinations$combi_id)), 0)
  n_plots <- sum(sapply(object@plots, length))
  n_cores <- object@.n_cores
  found_index <- sum(grepl(genome_id, dir(object@refdir, pattern = '.ebwt'))) > 0
  
  message(glue::glue('guideSet object of {genome_id}'))
  message(glue::glue('with {n_blacklisted} blacklisted regions'))
  message(glue::glue('with {n_tes} loci from {n_fams} families'))
  message(glue::glue('with {n_targets} targets'))
  message(glue::glue('with {n_guides} valid guides'))
  message(glue::glue('with {n_combinations} combinations'))
  message(glue::glue('with {n_plots} QC plots'))
  message(glue::glue('registered {n_cores} cores'))
  message(glue::glue('matching bowtie index found: {found_index}'))
})

#' Exports results from guideSet
#' 
#' @param guideSet guideSet containing the results
#' @param outdir String. Creates a new directory with timestamp prefix in \code{outdir} and exports results. If \code{NULL} and \code{force = TRUE}, the folder is created in the current working directory.
#' @param force Logical. If \code{TRUE} and \code{outdir = NULL}, writes output to new folder in current working directory.
#' @param workspace Logical. If \code{FALSE} (the default), suppresses additional export of \code{guideSet} as .RData file.
#' @param dpi Integer. Resolution of exported images. 
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13')
#' gs <- addGuides(gs, guide_length = 16, n_mismatches = 0, gc_content = c(0.25, 0.9), n_clust = 12)
#' gs <- plotGuides(gs)
#' export(gs, outdir = NULL, force = TRUE) # Creates new folder in current working directory and exports results  
#' }
#' @export
setGeneric('export', function(guideSet, ...) standardGeneric('export'), signature = 'guideSet') 

#' Create new guideSet object
#'
#' @param genome BSgenome object (required). Target genome assembly stored as BSgenome object.
#' @param alt_chromosomes Logical. If \code{FALSE} (the default), restricts genome annotation to the main chromosome assembly.
#' @param tes Path to repeatmasker output file or GRanges object with 'repname' metacolumn (required).
#' @param cis Path to bed file with cis regulatory feature coordinates or GRanges object (optional).
#' @param blacklist Path to bed file with blacklisted regions or GRanges object (optional). Guides binding to \code{blacklist} regions are blacklisted.
#' @param whitelist Path to bed file with whitelisted regions or GRanges object (optional). Guide off-target binding to \code{whitelist} regions are scored neutrally.
#' @param temp Path to directory where temporary files are stored. Needs to be read- and writeable with sufficient storage space for large file sizes. If \code{NULL} (the default), is set to the return value of \code{tempdir()}.
#' @param n_cores Integer. Number of cores to use for downstream functions. If \code{NULL} (the default), detects the number of cores automatically. The [doMC](https://cran.r-project.org/web/packages/doMC/index.html) packge must be installed to register the cores.
#' @param refdir Path to search for bowtie index files. Will create new indeces in \code{refdir} if no corresponding files are found (i.e. do not match BSgenome prefix). If empty (the default), searches in the bowtie_indeces directory of the Repguide installation path.
#' @param seed Integer. Seed for the random number generator. 19 by default.
#' @return guideSet object.
#' @examples
#' \dontrun{
#' # Path to directory containing BSgenome.Hsapiens.UCSC.hg38 bowtie indeces (e.g. BSgenome.Hsapiens.UCSC.hg38.1.ebwt, ...)
#' indexdir <-  system.file(package = 'Repguide', 'bowtie_indeces') 
#' # Path to TE annotation file
#' te_anno <- system.file(package = 'Repguide', 'extdata', 'hg38_ucsc_rmsk_ltr.txt.gz')
#' 
#' gs <- createGuideSet(genome = BSgenome.Hsapiens.UCSC.hg38, tes = te_anno, refdir = indexdir)
#' gs  
#' }
#' @seealso [BSgenome::available.genomes()], [bowtie manual](http://bowtie-bio.sourceforge.net/manual.shtml), and [repeatmasker](http://www.repeatmasker.org/)
#' @export
createGuideSet <- function(genome, 
                           alt_chromosomes = FALSE,
                           tes = NULL, 
                           cis = NULL, 
                           blacklist = NULL,
                           whitelist = NULL,
                           temp = NULL,
                           n_cores = NULL, 
                           refdir = '', 
                           seed = 19)
{
  new('guideSet', genome, alt_chromosomes, tes, cis, refdir, temp, blacklist, whitelist, n_cores, seed)
}                         