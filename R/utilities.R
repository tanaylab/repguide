# Consider using 'I4!=======44444+++++++' weighted phred from https://github.com/mhorlbeck/CRISPRiaDesign/blob/master/Library_design_walkthrough.md
# add guide mother seq performance could be improved using data.table
.bowtie <- function(guideSet,
                    n_mismatches = 0,
                    lower_count = 2)
{
  genome <- guideSet@genome
  index_dir <- guideSet@refdir
  n_cores <- guideSet@.n_cores
  kmers_file <- tempfile()
  align_file <- tempfile()
  
  kmers <- .jellyfish(guideSet, lower_count)
  fwrite(list(kmers), kmers_file)
  
  index_id <- genome@pkgname
  index_path <- paste0(index_dir, '/', index_id)
  
  if (!file.exists(paste0(index_path, '.1.ebwt'))) { .bowtieIndex(guideSet) }
  
  print(paste0('Aligning kmers against: ', index_id))
  Rbowtie::bowtie(sequences = kmers_file,
                  index = index_path, 
                  type = 'single',
                  r = TRUE, tryhard = TRUE, a = TRUE, suppress = c(6), threads = n_cores, v = n_mismatches, # -a write all alignments, k max number valid alignments, -v determines n mismatches
                  outfile = align_file, 
                  sam = FALSE,
                  force = TRUE)
                  
  kmers_mapped <- importKmers(align_file)
  
  kmers_mapped$guide_seq <- # add guide mother seq
    left_join(kmers_mapped %>% as_tibble, 
              tibble(guide_seq = kmers, kmer_id = 0:(length(kmers) -1)), 
              by = 'kmer_id') %>%
    pull(guide_seq)
    
  guideSet@kmers <- kmers_mapped
              
  return(guideSet)
}

.bowtieIndex <- function(guideSet)
{
  genome <- guideSet@genome
  index_id <- genome@pkgname
  index_dir <- guideSet@refdir
    
  print ('Creating bowtie index, this may take a while')
  genome_fasta_fn <- tempfile()
    
  rtracklayer::export(genome, genome_fasta_fn, format = 'fasta', compress = FALSE)
  #system(glue('chmod 777 {genome_fasta_fn}'))
  
  # --noref don't buildt ebwt 3/4 - paired mapping | lower offrate increases performance for multi mappers
  tmp <- Rbowtie::bowtie_build(references = genome_fasta_fn, outdir = index_dir, prefix = index_id, force = TRUE, noref = TRUE, offrate = 3) 
}

.keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  return(genome)
}

.kmerStats <- function(kmers,
                       full = FALSE)
{
  kmers_dt <- as.data.table(kmers)
  
  score_stats <- kmers_dt[, .(Son = max(Son),
                              Soff = max(Soff)), by = c('kmer_id', 'unique_id')
                        ][, .(Son_tot = sum(Son),
                              Soff_tot = sum(Soff)), by = 'kmer_id']
  if (full)
  {
   score_stats <- score_stats[kmers_dt, on = 'kmer_id']
  }
  
  return(score_stats)
}

msaToLong <- function(msa)
{
  n_loci <- length(msa)
  
  if (is.null(names(msa)) | sum(duplicated(names(msa))) > 0)
  {
    row_labels <- as.character(1:n_loci)
    names(msa) <- row_labels
  } else {
    row_labels <- names(msa)
  }
  
  
  m <- as.matrix(msa)
  colnames(m) <- 1:ncol(m)
  
  #sm <- Matrix::Matrix(m, sparse = TRUE)
  msa_long <- 
    m %>% 
    as.data.frame %>% 
    tibble::rownames_to_column() %>% 
    gather(pos, base, -rowname) %>% 
    rename(te_id = rowname) %>% 
    as_tibble %>% 
    #filter(base != '-') %>%
    mutate(pos = as.numeric(pos),
           te_id = factor(te_id, levels = row_labels)) %>%
    arrange(te_id, pos)
  
  # msa_long <- 
    # summary(sm) %>% 
    # as_tibble %>%
    # dplyr::rename(te_id = i, pos = j, gap = x) %>%
    # mutate(te_id = row_labels[te_id], gap = FALSE) %>%
    # arrange(te_id, pos)
	
  # Add column with ungapped pos integer
  msa_long <-
    msa_long %>%
    group_by(te_id) %>%
    mutate(pos_wo_gaps = 1:n()) %>%
    ungroup

  # Add gap perc
  msa_long <-
    msa_long %>% 
    group_by(pos) %>% 
    mutate(gap_perc = 1 - sum(base != '-') / n_loci) %>% 
    ungroup
    
  # Add CpG track  
  msa_long <-
    msa_long %>% 
    filter(base != '-') %>%
    mutate(CpG = lag(base == 'C', n = 1) & base == 'G', 
           CpG = ifelse(lead(base == 'G', n = 1) & base == 'C', TRUE, CpG)) %>%
    select(te_id, pos, CpG) %>%
    left_join(msa_long, .) %>%
    replace_na(list(CpG = FALSE)) %>%
    mutate(CpG = ifelse(base == '-', NA, CpG))
  
  return(msa_long)
}

binGenome <- function(genome, bin_width = 250)
{
  genomic_bins <- tile(GRanges(seqinfo(genome)), width = bin_width) %>% unlist
  names(genomic_bins) <- 1:length(genomic_bins)
  
  return(genomic_bins)
}

#' List transposable element families.
#'
#' Retrieves all available families in the provided guideSet that match defined characteristics.
#'
#' @param guideSet guideSet object to query.
#' @param class Character. Only show families belonging to \code{class}.
#' @param low_complexity Logical. Should low complexity families be returned.
#' @return Character vector.
#' @export
listFamilies <- function(guideSet, 
                         repclass = NULL, 
                         low_complexity = FALSE)
{
  tes <- as_tibble(guideSet@tes) %>% select(repname, repclass) %>% distinct
  if (!is.null(repclass)) 
  {
    tes <- filter(tes, repclass %in% repclass)
  }
  
  if (!low_complexity)
  {
    tes <- filter(tes, repclass != 'Simple_repeat')
  }
  
  res <- tes %>%  pull(repname) %>% sort
  return(res)
}

clustGuides <- function(guideSet, 
                        min_Son = 0,
                        n_clust = 15)
{
  if(n_clust > 20) { stop('Maximal 20 clusters currently supported') }
  
  kmers <- as_tibble(guideSet@kmers) %>% select(-matches('kmer_clust|te_clust'))
  kmers_filt <- 
    kmers %>%
    filter(Son > min_Son & on_target == 1) %>%
    filter(valid)
    
  mat_full <- 
    kmers_filt %>%
    select(kmer_id, te_id, Son) %>%
    .tidyToSparse()
    
  #mat_full = log2(mat_full+1)
   
  # mat_slim <- kmers %>%
    # filter(on_target >= 0) %>%
    # mutate(on_target = on_target * Sbind) %>%
    # select(kmer_id, te_id, on_target) %>%
    # guideR::tidyToSparse()  

  print(paste0('Clustering ', nrow(mat_full), ' kmers into ', n_clust, ' groups'))
  kmer_cors <- as.matrix(qlcMatrix::cosSparse(t(mat_full)))
  #kmer_cors <- tgs_cor(as.matrix(t(mat_full)), spearman = TRUE)
  kmer_clusts <- tibble(kmer_id = as.numeric(rownames(mat_full)),
                        kmer_clust = as.numeric(cutree(fastcluster::hclust(tgs_dist(kmer_cors), 'ward.D2'), n_clust)))
                          
  print(paste0('Clustering ', ncol(mat_full), ' loci into ', n_clust, ' groups'))
  loci_cors <- as.matrix(qlcMatrix::cosSparse(mat_full))
  loci_clusts <- tibble(te_id = as.numeric(colnames(mat_full)),
                        te_clust = as.numeric(cutree(fastcluster::hclust(tgs_dist(loci_cors), 'ward.D2'), n_clust)))   

  kmers <- left_join(kmers, kmer_clusts, by = 'kmer_id') %>% left_join(., loci_clusts, by = 'te_id')                        
  guideSet@kmers$kmer_clust <- kmers$kmer_clust
  guideSet@kmers$te_clust <- kmers$te_clust
  
  return(guideSet)
}

.jellyfish <- function(guideSet, lower_count = 2)
{
  PAM <- guideSet@PAM
  kmer_length <- guideSet@guide_length + nchar(PAM)
  n_cores <- guideSet@.n_cores
  print ('Computing kmer universe')
  
  seq_file <- tempfile()
  kmers_file <- tempfile()
  
  # Write target seqs to file for jellyfish
  Biostrings::writeXStringSet(guideSet@targets$seq, filepath = seq_file, format = 'fasta')
  
  jellyfish_path <- system.file(package = 'Repguide', 'inst', 'bin', 'jellyfish-linux')
  
  cmd <- glue::glue("{jellyfish_path} count --mer-len {kmer_length} --size 100M --threads {n_cores} --out-counter-len 0 --lower-count {lower_count} --text {seq_file} --output {kmers_file} ")
  system(command = cmd)
  
  kmers_all <- fread(kmers_file, skip = 1)
  print('Select kmers with proper PAM')
  kmers_filt <-
    kmers_all %>%
    filter(grepl('[A,T,G,C]GG$', V1)) %>% # NGG PAM
    as_tibble %>%
    pull(V1)
    
  return(kmers_filt)
}   

makeGRangesFromDataFramePar <- function(df, keep.extra.columns = FALSE, n_cores = NULL)
{
  if (sum(colnames(df) %in% 'seqnames') == 0) { df <- df %>% rename(seqnames = chrom) }
  if (is.null(n_cores)) { n_cores <- max(parallel::detectCores(), length(unique(df$seqnames))) }
  
  doMC::registerDoMC(n_cores)
  
  gr <- 
    plyr::dlply(df, 
                .variable = 'seqnames', 
                .fun = makeGRangesFromDataFrame,
                keep.extra.columns = keep.extra.columns,
                .parallel = TRUE) %>%
    GRangesList %>%
    unlist(., recursive = TRUE, use.names = FALSE)
  
  return(gr)
}

.rmGaps = function(alignment, max_gap_freq = 0.8)
	{
	alignment_bin = ape::as.DNAbin(alignment)
	alignment_bin_wo_gaps = ape::del.colgapsonly(as.matrix(alignment_bin), max_gap_freq)
	
	res = as.character(as.list(alignment_bin_wo_gaps)) %>% 
		lapply(., paste0, collapse = '') %>% 
		unlist %>% 
		DNAStringSet
	return(res)
}

.tidyToSparse = function(df_tidy)
{
	colnames(df_tidy) = c('a', 'b', 'c')
	data = df_tidy %>%
		mutate_at(c('a', 'b'), funs(factor(.)))
		
	data_sparse = Matrix::sparseMatrix(as.integer(data$a), as.integer(data$b), x = data$c)
	colnames(data_sparse) = levels(data$b)
	rownames(data_sparse) = levels(data$a)
	
	return(data_sparse)
}

kmersToCov <- function(kmers, guide_length = 19) # input Kmer annotation df
{
  te_id_cov <-
    kmers %>% 
    select(te_id, pos_in_loc) %>% 
    mutate(pos = 1, rep = guide_length) %>% 
    uncount(rep) %>% 
    mutate(pos = rep(1:guide_length, times = nrow(.)/guide_length), 
           pos = pos + pos_in_loc - 1) %>% 
      select(-pos_in_loc) %>% 
    count(te_id, pos)
	
  return(te_id_cov)
}

.tibToM = function(tib) # takes as input a tibble with first column as rownames
{
	m = as.matrix(tib[,-1])
	rownames(m) = pull(tib, 1)
	return(m)
}
 
cols_long <- c("#984EA3", 
               "#771155", 
               "#AA4488", 
               "#CC99BB", 
               "#114477", 
               "#4477AA", 
               "#77AADD", 
               "#117777", 
               "#44AAAA", 
               "#77CCCC", 
               "#117744", 
               "#44AA77", 
               "#88CCAA",
               "#777711", 
               "#AAAA44", 
               "#DDDD77", 
               "#774411", 
               "#AA7744", 
               "#DDAA77", 
               "#771122", 
               "#AA4455", 
               "#DD7788", 
               "darkred", 
               "darkorange2", 
               "orange", 
               "gold", 
               "pink") 