#' @include guideSet.R
NULL

setMethod("export", signature('guideSet'), function(guideSet, 
                                                    outdir = NULL, 
                                                    full = FALSE,
                                                    force = FALSE, 
                                                    workspace = FALSE, 
                                                    dpi = 300) 
{ 
  if(is.null(outdir) & !force) 
  { 
    stop('No outdir provided. Either specify the outdir or use force = TRUE to save reports to the current working directory') 
  }
  oldwd <- getwd()
  
  # create output directory
  datum <- gsub('-', '', Sys.Date())
  zeit <- gsub(' ', '', gsub(':', '', format(Sys.time(), "%X")))

  suffix <- paste0('Repguide_report', '_', datum, '_', zeit, '/')
  outdir <- ifelse(is.null(outdir), suffix, paste0(outdir, '/', suffix))
  dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  dir.create(paste0(outdir, '/full_stats/'))
  setwd(outdir)
  message(paste0('Exporting guideSet to ', outdir))
  
  ##############
  # full report
  ##############
  if (full)
  {
    # export targets
    if(length(guideSet@targets) > 0)
    {
      data.table::fwrite(as_tibble(guideSet@targets), file = paste0('full_stats/', 'targets.txt'), sep = '\t')
    }

    # export kmers
    kmers <- as_tibble(guideSet@kmers)
    if(length(kmers) > 0)
    {
      data.table::fwrite(kmers, file = paste0('full_stats/', 'kmers.txt'), sep = '\t')
    }
     
    # export combinations
    if(length(guideSet@combinations) > 0)
    {
      data.table::fwrite(tidyr::unnest(guideSet@combinations), file = paste0('full_stats/', 'combinations.txt'), sep = '\t')
    }
  }
  
  ##############
  # slim report
  ##############
  
  # export guide seqs
  .exportGuides(guideSet)
  
  # export plots
  .exportPlots(guideSet, dpi = dpi)
  
  # export workspace
  if(workspace) { save(guideSet, file = 'guideSet') }
  
  # export history
  savehistory(file = 'command_history.txt')
  
  setwd(oldwd) 
})     

.exportPlots <- function(guideSet, dpi = 300)
{
  # n_fams <- length(guideSet@families)
  # n_guides <- max(guideSet@combinations$n_guides)
 
 # export target plots
 i <- 1
 for (p in guideSet@plots$targets)
 {
    fn <- paste0('targets', '_', i, '.png')
    ggsave(fn, p, device = 'png', dpi = dpi)
    i <- i + 1
 }
 
 # export guide plots
 i <- 1
 for (p in guideSet@plots$guides)
 {
    fn <- paste0('guides', '_', i, '.png')
    ggsave(fn, p, device = 'png', dpi = dpi)
    i <- i + 1
 }
 
 # export guide plots
 i <- 1
 for (p in guideSet@plots$combinations)
 {
    fn <- paste0('combinations', '_', i, '.png')
    ggsave(fn, p, device = 'png', dpi = dpi)
    i <- i + 1
 }
}

.exportGuides <- function(guideSet)
{
  combinations_subs <- guideSet@combinations %>% filter(best) %>% unnest
  kmers <- as_tibble(guideSet@kmers)
  kmers_subs <- inner_join(kmers, combinations_subs, by = 'kmer_id') 
    
  kmers_by_nguides <- 
    kmers_subs %>%
    group_by(n_guides) %>% 
    do(data.frame = as_tibble(.))
    
  kmer_seqs_by_nguides <-
    kmers_subs %>%
    group_by(n_guides) %>%
    do(guide_seq = unique(DNAStringSet(structure(.$guide_seq, names = as.character(.$kmer_id)), use.names = TRUE)))
  
  for (i in kmers_by_nguides$n_guides)
  {
    data.table::fwrite(kmers_by_nguides$data.frame[[i]], 
                       file = paste0(i, '_guides_binding.txt'), sep ='\t')
    Biostrings::writeXStringSet(kmer_seqs_by_nguides$guide_seq[[i]], 
                                filepath = paste0(i, '_guides_sequence.fasta'), format = 'fasta')  
  }
}