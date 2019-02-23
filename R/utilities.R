# Consider using 'I4!=======44444+++++++' weighted phred from https://github.com/mhorlbeck/CRISPRiaDesign/blob/master/Library_design_walkthrough.md
# add guide mother seq performance could be improved using data.table
.bowtie <- function(guideSet,
                    n_mismatches = 0)
{
  genome <- guideSet@genome
  index_dir <- guideSet@refdir
  n_cores <- guideSet@.n_cores
  guide_length <- guideSet@guide_length
  kmers_file <- tempfile(tmpdir = guideSet@temp)
  align_file <- tempfile(tmpdir = guideSet@temp)
  
  kmers <- guideSet@kmers$guide_seq
  data.table::fwrite(list(kmers), kmers_file)
  
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

  # Throw binding sites with NGG mismatch
  base_pos <- paste(c(guide_length + 1, guide_length + 2) , collapse = '|')
  kmers_mapped <- kmers_mapped[!stringr::str_detect(kmers_mapped$mismatches, base_pos)] # 0-based!
  
  # delete 'N' from NGG mismatches
  base_pos <- paste0(guide_length, ':[ACGT]>[ACGT]')
  kmers_mapped$mismatches <- stringr::str_replace(kmers_mapped$mismatches, base_pos, '')
  
  # Add number of mismatches
  kmers_mapped$n_mismatches <- stringr::str_count(kmers_mapped$mismatches, stringr::fixed(':'))
  
  kmers_mapped$guide_seq <- # add guide mother seq
    left_join(as_tibble(kmers_mapped), 
              tibble(guide_seq = substring(kmers, 1, guide_length), kmer_id = 0:(length(kmers) -1)), 
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
  genome_fasta_fn <- tempfile(tmpdir = guideSet@temp)
    
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
    dplyr::rename(te_id = rowname) %>% 
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
    left_join(msa_long, ., by = c('te_id', 'pos')) %>%
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

.getTSS <- function(genes)
{
	transcripts = genes[genes$type == 'transcript']
	end(transcripts[strand(transcripts) == '+']) = start(transcripts[strand(transcripts) == '+'])
	#start(transcripts[strand(transcripts) == '+']) = start(transcripts[strand(transcripts) == '+']) - extend
	
	start(transcripts[strand(transcripts) == '-']) = end(transcripts[strand(transcripts) == '-'])
	#end(transcripts[strand(transcripts) == '-']) = end(transcripts[strand(transcripts) == '-']) + extend
	tss = unique(transcripts)
	return(tss)
}

.jellyfish <- function(guideSet, lower_count = 2)
{
  PAM <- guideSet@PAM
  kmer_length <- guideSet@guide_length + nchar(PAM)
  n_cores <- guideSet@.n_cores
  print ('Computing kmer universe')
  
  seq_file <- tempfile(tmpdir = guideSet@temp)
  kmers_file <- tempfile(tmpdir = guideSet@temp)
  
  # Write target seqs to file for jellyfish
  Biostrings::writeXStringSet(guideSet@targets$seq, filepath = seq_file, format = 'fasta')
  
  jellyfish_path <- system.file(package = 'Repguide', 'bin', ifelse(Sys.info()['sysname'] == 'Linux', 'jellyfish-linux', 'jellyfish-macosx'))
  
  cmd <- glue::glue("{jellyfish_path} count --mer-len {kmer_length} --size 100M --threads {n_cores} --out-counter-len 0 --lower-count {lower_count} {seq_file} --output {kmers_file} ")
  system(command = cmd)
  
  kmers_all <- data.table::fread(kmers_file, skip = 1)
  print('Select kmers with proper PAM')
  kmers_filt <-
    kmers_all %>%
    filter(grepl('[A,T,G,C]GG$', V1)) %>% # NGG PAM
    as_tibble %>%
    pull(V1)
    
  guideSet@kmers <- GRanges(seqnames = 1:length(kmers_filt), 
                            ranges = 1:length(kmers_filt), 
                            guide_seq = kmers_filt)
  return(guideSet)
}   

makeGRangesFromDataFramePar <- function(df, keep.extra.columns = FALSE, n_cores = NULL)
{
  if (sum(colnames(df) %in% 'seqnames') == 0) { df <- df %>% dplyr::rename(seqnames = chrom) }
  if (is.null(n_cores)) { n_cores <- max(parallel::detectCores(), length(unique(df$seqnames))) }
  
  #doMC::registerDoMC(n_cores)
  
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

.plotEmpty <- function(string = 'NA', size = 3, border = TRUE)
{
  p <- 
    ggplot() + 
      annotate('text', x = 4, y = 25, size = size, label = string) +
      theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
      
  if (!border)
  {
    p <- p + theme(panel.border = element_blank())
  }
  return(p)
}

.rmGaps <- function(alignment, max_gap_freq = 0.8)
	{
	alignment_bin = ape::as.DNAbin(alignment)
	alignment_bin_wo_gaps = ape::del.colgapsonly(as.matrix(alignment_bin), max_gap_freq)
	
	res = as.character(as.list(alignment_bin_wo_gaps)) %>% 
		lapply(., paste0, collapse = '') %>% 
		unlist %>% 
		DNAStringSet
	return(res)
}

.tidyToSparse <- function(df_tidy)
{
	colnames(df_tidy) = c('a', 'b', 'c')
	data = df_tidy %>%
		mutate_at(c('a', 'b'), funs(factor(.)))
		
	data_sparse = Matrix::sparseMatrix(as.integer(data$a), as.integer(data$b), x = data$c)
	colnames(data_sparse) = levels(data$b)
	rownames(data_sparse) = levels(data$a)
	
	return(data_sparse)
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