#' @include guideR.R
NULL

setMethod("export", signature('guideSet'), function(guideSet, outdir = NULL, force = FALSE, workspace = FALSE, dpi = 300) 
{ 
  if(is.null(outdir) & !force) 
  { 
    stop('No outdir provided. Either specify the outdir or use force = TRUE to save reports to the current working directory') 
  }
  oldwd <- getwd()
  
  # create output directory
  datum <- gsub('-', '', Sys.Date())
  zeit <- gsub(' ', '', gsub(':', '', format(Sys.time(), "%X")))

  suffix <- paste0('guideR_report', '_', datum, '_', zeit, '/')
  outdir <- ifelse(is.null(outdir), suffix, paste0(outdir, '/', suffix))
  dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  setwd(outdir)
  message(paste0('Exporting guideSet to ', outdir))
  
  # export targets
  if(length(guideSet@targets) > 0)
  {
    data.table::fwrite(as_tibble(guideSet@targets), file = 'targets.txt', sep = '\t')
  }

  # export kmers
  if(length(guideSet@kmers) > 0)
  {
    data.table::fwrite(as_tibble(guideSet@kmers), file = 'kmers.txt', sep = '\t')
  }
 
  # export combinations
  if(length(guideSet@combinations) > 0)
  {
    data.table::fwrite(tidyr::unnest(guideSet@combinations), file = 'combinations.txt', sep = '\t')
  }
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