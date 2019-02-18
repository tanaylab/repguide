.plotCpG <- function(guideSet,
                     window_size = 250)
{
  if(length(guideSet@alignments) == 0)
  {
    p_cpg_dens <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = 'No alignment provided')
  } else {
    alignment <- guideSet@alignments
    
    k <- window_size
    
    liste <- list()
    for(repname_id in names(alignment))
    {
      liste[[repname_id]] <- 
        sapply(alignment[[repname_id]], function(x)
        {
          x <- as.character(x)
          stringr::str_count(substring(x, seq(1, nchar(x), 1), seq(1, nchar(x), 1) + k -1), 'CG') / k
        }) %>%
        as_tibble %>%
        mutate(pos = 1:nrow(.),
               repname = repname_id) %>%
        gather('te_id', 'CpG', -pos, -repname)
    }
    
    cpg_density <- 
      do.call(bind_rows, liste) %>%
      group_by(pos, repname) %>% 
      summarise(CpGmean = mean(CpG),
                CpGvar = sd(CpG)) %>%
      ungroup %>%
      mutate(CpGupper = CpGmean + CpGvar,
             CpGlower = CpGmean - CpGvar)
    
    .custom_breaks <- function(x) { c(0, max(x)) }
    .custom_labels <- function(x) { ifelse(x == 0, '0', scales::scientific(x)) }
    
    p_cpg_dens <-
      cpg_density %>% 
      ggplot(aes(pos, CpGmean, col = repname)) + 
        geom_line() +
        geom_ribbon(aes(ymin=CpGlower, ymax=CpGupper, fill = repname), linetype=2, alpha=0.05, colour = NA) +
        facet_wrap(~repname, ncol = 1, scales = 'free') +
        scale_y_continuous(breaks = .custom_breaks, labels = .custom_labels) +
        theme(legend.position = 'none',
              #plot.margin = unit(c(0, 0, 0, 0), "cm"),
              strip.background = element_blank(),
              strip.text.x = element_blank()) +
        ggplot2::guides(col = guide_legend(ncol = 1))
  }
  return(p_cpg_dens)
}

#' @export
plotMSA <- function(guideSet, 
                    # max_gap_perc = 0.8,
                    type = 'base')
{
  if(length(guideSet@alignments) == 0)
  {
    p_msa <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = 'No alignment provided')
  } else {
    msa <- guideSet@alignments
    if (class(msa) == 'DNAStringSetList')
    {
      msa <- 
        lapply(msa, function(x) { res <- msaToLong(x); return(res) }) %>%
        data.table::rbindlist(., idcol = 'repname')
    }
    
    # msa_filt <-
      # msa %>%
      # filter(gap_perc <= max_gap_perc) %>%
      # group_by(te_id) %>%
      # mutate(pos = 1:n()) %>%
      # filter(base != '-') %>%
      # ungroup
    
    if(type == 'CpG')
    {
      p_msa <-  
        msa %>%   
        ggplot(aes(pos, te_id, fill = CpG)) + 
          geom_raster() + 
          scale_fill_manual(values = c('lightgrey', 'black'))
    } 
    if(type == 'base')
    {
      p_msa <-  
        msa %>%   
        ggplot(aes(pos, te_id, fill = base)) + 
          geom_raster() + 
          scale_fill_manual(values = c('A' = "#5DA731", 'T' = "#D63317", 'G' = "#DFD93C", 'C' = "#2F4C9B", '-' = 'white', 'N' = 'lightgrey'))
    }
    
    p_msa <-
      p_msa +
      facet_wrap(~repname, ncol = 1, scales = 'free', strip.position = 'left') +
      theme(#plot.margin = unit(c(0, 0, 0, 0), "cm"),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(),
            legend.position = 'none')
  }     
  return(p_msa)
}

# plotTargetCov <- function(guideSet)
# {
  # kmers <- as_tibble(guideSet@kmers) %>% filter(on_target == 1)
  # combinations <- guideSet@combinations %>% filter(best) %>% unnest
  # families <- guideSet@families
  
  # ggdata <-
    # kmers %>% 
    # select(repname, n_mismatches, te_id, kmer_id) %>% 
    # distinct %>% 
    # left_join(., combinations %>% select(combi_id, kmer_id)) %>% 
    # complete(combi_id, n_mismatches, te_id, repname) %>% 
    # filter(!is.na(combi_id)) %>% 
    # mutate(kmer_id = !is.na(kmer_id)) %>% 
    # rename(covered = kmer_id) %>% 
    # count(n_mismatches, combi_id, repname, te_id, covered, sort = TRUE) %>% 
    # count(n_mismatches, combi_id, repname, covered, n) %>%
    # mutate(n = ifelse(!covered, 0, n),
           # n = ifelse(n > 5, '>5', n), 
           # n = factor(n, levels = c('>5', 5:0)))
  
  # p_target_cov <-
    # ggdata %>%
      # ggplot(aes(repname, nn, fill = n)) + 
        # geom_bar(stat = 'identity', position = 'fill') +
        # facet_grid(n_mismatches ~ combi_id) +
        #scale_fill_manual(values = structure(c(colorRampPalette(c('darkred', 'orange'))(6), 'lightgrey'), names = c('>5', '5', '4', '3', '2', '1', '0'))) +
        # theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  # return(p_target_cov)
# }

.plotClusts <- function(guideSet)
{
  kmers <- 
    as_tibble(guideSet@kmers) %>%
    filter(!is.na(kmer_clust) & !is.na(te_clust)) %>%
    filter(valid)
    
  if (nrow(kmers) > 0) 
  {
    #########################
    # Sample Kmers and Loci #
    #########################
    if (nrow(kmers) > 5e6)
    {
      set.seed(guideSet@.seed)
      kmers_slim <- kmers %>% select(kmer_id, te_id, kmer_clust, te_clust, best)
      
      kmers_subset <- 
        kmers_slim %>% 
        select(kmer_id, kmer_clust, best) %>%
        distinct %>%
        group_by(kmer_clust) %>% 
        summarise(kmer_id = list(c(sample(kmer_id, round(n()/10)), kmer_id[best]))) %>% 
        unnest %>% 
        pull(kmer_id)

      loci_subset <- 
        kmers_slim %>% 
        select(te_id, te_clust) %>%
        distinct %>%
        group_by(te_clust) %>% 
        summarise(te_id = list(c(sample(te_id, round(n()/10))))) %>% 
        unnest %>% 
        pull(te_id)  
        
      kmers <- kmers %>% filter(te_id %in% loci_subset & kmer_id %in% kmers_subset)
      xaxis_lab <-'Sampled '
      yaxis_lab <-'Sampled '
    } else {
      xaxis_lab <-''
      yaxis_lab <-''
    }
    
    kmers_dt <- as.data.table(kmers, sorted = FALSE)
    setkey(kmers_dt, 'kmer_id', 'te_id')
    kmer_hit_score <- kmers_dt[, .(Son = max(Son),
                                   kmer_clust = kmer_clust[1],
                                   te_clust = te_clust[1]),
                                   by = c('kmer_id', 'te_id')]
    
    # ggdata <-
      # kmer_hit_score[order(kmer_clust)
                   # ][, kmer_id := forcats::fct_inorder(as.character(kmer_id))
                   # ][order(te_clust)
                   # ][, te_id := forcats::fct_inorder(as.character(te_id))
                   # ][, Slocus := cut(Son, breaks = c(0, 0.25, 0.5, 1), include.lowest = TRUE)]
    
    
    ggdata <- # mildly slow
      kmer_hit_score %>%
      arrange(kmer_clust) %>% 
      mutate(kmer_id = forcats::fct_inorder(as.character(kmer_id))) %>% 
      arrange(te_clust) %>% 
      mutate(te_id = forcats::fct_inorder(as.character(te_id))) %>%
      mutate(Slocus = cut(Son, breaks = c(0, 0.25, 0.5, 1), include.lowest = TRUE))
      
    # Plotting
    row_clust_lines <- c(0, ggdata %>% count(kmer_id, kmer_clust) %>% count(kmer_clust) %>% pull(n) %>% cumsum)
    col_clust_lines <- c(0, ggdata %>% count(te_id, te_clust) %>% count(te_clust) %>% pull(n) %>% cumsum)
    
    n_kmers <- length(unique(ggdata$kmer_id))
    n_loci <- length(unique(ggdata$te_id))
    
    best_kmer_pos = which(levels(ggdata$kmer_id) %in% (kmers %>% filter(best) %>% pull(kmer_id) %>% unique))
    
    p_heatmap <- 
      ggdata %>% 
      ggplot(aes(te_id, kmer_id, fill = Slocus)) +
        geom_raster() + 
        scale_fill_manual(values = c('orange', 'darkred', 'purple'), drop = FALSE) +
        #scale_fill_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', midpoint = 0)
        annotate('segment', x = 0, xend = n_loci, y = row_clust_lines, yend = row_clust_lines, lwd = 0.1) +
        geom_vline(xintercept = col_clust_lines, lwd = 0.1) +
        xlab(paste0(xaxis_lab, 'target loci (n = ', n_loci, ')')) +
        ylab(paste0(yaxis_lab, 'guides (n = ', n_kmers, ')')) +
        coord_cartesian(xlim = c(0, n_loci + round(n_loci / 100)), clip="off") +
        annotate('segment', x = n_loci + round(n_loci / 100), xend = n_loci + 10, y = best_kmer_pos, yend = best_kmer_pos, 
                 lwd = 0.1, arrow = arrow(type = "open", length = unit(n_kmers/100000, "npc"))) +
        #coord_fixed() +
        theme(legend.position = 'bottom',
              legend.key.width = unit(.3,"cm"),
              legend.key.height = unit(.1,"cm"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank()) +
              #plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ggplot2::guides(fill = guide_legend(nrow = 1))
  } else {
    p_heatmap = ggplot() + annotate("text", x = 4, y = 25, size=8, label = 'No clustering provided')
  }
  
  return(p_heatmap)
}

.plotTSS <- function(guideSet) # TSS on consensus
{
  if(length(guideSet@alignments) == 0)
  {
    p_tss_cov <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = 'No alignment provided')
  } else {
    genome <- guideSet@genome
    cons <- guideSet@consensus
    families <- names(guideSet@consensus)
    targets <- guideSet@targets
    
    tss <- resize(guideSet@cis, 20, fix = 'center')
    tss_subs <- subsetByOverlaps(tss, targets)
    
    # Add corresponding repname and seq
    tss_subs$repname <- targets[findOverlaps(tss_subs, targets)@to]$repname
    tss_subs$seq <- getSeq(genome, tss_subs)
    if (length(tss_subs == 0)) { stop('No overlapping TSS found') }
    
    liste <- list()
    for (fam in families)
    {
      tss_seqs <- tss_subs[tss_subs$repname == fam]$seq
      con <- cons[fam]
     
      # Align TSS seqs against consensus
      tss_cov <- as.data.frame(coverage(pairwiseAlignment(tss_seqs, con, type = 'global-local'))) %>%
        rownames_to_column %>%
        as_tibble %>%
        dplyr::rename(pos = rowname, n = value) %>%
        mutate(repname = fam,
               pos = as.numeric(pos))
      liste[[fam]] = tss_cov
    }
    ggdata <- do.call(bind_rows, liste)
    
    .custom_breaks <- function(x) { c(0, floor(x)) }  
    p_tss_cov <- 
      ggdata %>% 
        ggplot(aes(pos, n, fill = repname)) + 
          geom_density(stat = 'identity', alpha = 0.5) +
          ggplot2::guides(fill = guide_legend(ncol = 1)) +
        # geom_bar(stat = 'identity') + 
          facet_wrap(~repname, ncol = 1, scales = 'free', strip.position = 'right') +
          scale_y_continuous(breaks = .custom_breaks) +
          theme(legend.position = 'none',
                #plot.margin = unit(c(0, 0, 0, 0), "cm"),
                strip.background = element_blank(),
                strip.text.y = element_blank())
  }
  return(p_tss_cov)
}

#' Generates targets QC plots
#' 
#' Computes multiple QC plots and arranges into a figure with subpanels:
#' \describe{
#'   \item{(a)}{ Genomic copy-number vs average size (in kbs) of selected (red) vs all other families in guideSet. }
#'   \item{b}{ Genomic copy-number of selected families on autosomes and sex chromosomes (other refers to alternative chromosome assemblies). }
#'   \item{c}{ Size distribution of selected families. }
#'   \item{d}{ Multiple sequence alignment of selected families. Rows represent individual loci and columns position on the alignment. }
#'   \item{e}{ CpG density along the consensus model. Positions are identical to (d). }
#'   \item{f}{ Position on the consensus model (as in e) of putative cis regulatory features. }
#' }
#' @param guideSet guideSet containing combinations.
#' @return list of ggplot objects.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13A')
#' gs <- addGuides(gs, guide_length = 16)
#' gs <- plotGuides(gs)  
#' }
#' @seealso [addTargets()], and [export()]
#' @export
plotTargets <- function(guideSet)
{
  n_fams <- length(guideSet@families)
  
  # Plot global overview of selected families
  p_global_overview <-
    guideSet@tes %>% 
    as_tibble %>% 
    group_by(repname) %>% 
      summarise(copy_number = n(), mean_kbs = mean(width) / 1000, repclass = repclass[1]) %>% 
    mutate(label = repname %in% guideSet@families) %>% 
    ggplot(aes(copy_number, mean_kbs, col = label)) + 
      geom_point(aes(size = label, alpha = label)) + 
      scale_x_log10() + 
      scale_y_log10() + 
      #ggrepel::geom_label_repel(point.padding = 0.5, max.iter = 50, show.legend = FALSE) +
      scale_color_manual(values = c('grey', 'red')) +
      scale_alpha_manual(values = c(0.1, 1)) +
      scale_size_manual(values = c(0.1, 2)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
            legend.position = 'none')
 
  # Plots number of loci per family (colored by auto/sex chromosome)
  p_cpn <- 
    guideSet@targets %>% 
    as_tibble %>%
    mutate(chrom = ifelse(grepl('X|Y', seqnames), as.character(seqnames), 'autosome'),
           #chrom = factor(chrom, levels = c('auto', 'chrX', 'chrY')),
           chrom = forcats::fct_lump(chrom, n = 3, other_level = 'other')) %>% 
    count(repname, chrom) %>% 
    ggplot(aes(repname, n, fill = chrom)) + 
      geom_bar(stat = 'identity') +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
            axis.text.x = element_text(angle=90, hjust=1),
            axis.title.x = element_blank(),
            legend.position = 'right',
            legend.key.width = unit(.4,"cm"),
            legend.key.height = unit(.2,"cm"),
            legend.title = element_blank()) +
      ggplot2::guides(fill = guide_legend(ncol = 1))
      
  # Plots the size distribution (density) per family    
  p_size_distr <-
    guideSet@targets %>% 
    as_tibble %>%
    ggplot(aes(x = repname, y = width, group = repname, fill = repname)) + 
      geom_violin() +
      geom_boxplot(width=.1, outlier.colour=NA) +
      scale_y_log10() +
      theme(legend.position = 'none',
            axis.text.x = element_text(angle=90, hjust=1),
            axis.title.x = element_blank())
            
  p_msa <- plotMSA(guideSet)
  
  p_tss_cov <- .plotTSS(guideSet)
  
  p_cpg <- .plotCpG(guideSet)
    
  plots <- list('a' = p_global_overview, 
                'b' = p_cpn, 
                'c' = p_size_distr, 
                'd' = p_msa, 
                'e' = p_cpg, 
                'f' = p_tss_cov)  
  
  plots_combined <- cowplot::plot_grid(plotlist = plots, labels = 'auto', ncol = 3, rel_heights = c(1, 2))
                                
  guideSet@plots$targets <- plots
  guideSet@plots$targets$figure <- plots_combined
  
  print(plots_combined)
  
  return(guideSet)
}

#' Generates guides QC plots
#' 
#' Computes multiple QC plots and arranges into a figure with subpanels:
#' \describe{
#'   \item{a}{ Number of guides passing blacklisting steps. }
#'   \item{b}{ Guide GC content distribution. Blacklisted guides are highlighted in red. }
#'   \item{c}{ Total on- and off-target score per valid (turquoise) and blacklisted (red) guides. }
#'   \item{d}{ Predicted binding score distribution per number of mismatch to genomic complement. }
#'   \item{e}{ Distance to the nearest cis regulatory feature (defined by \code{cis} in the \code{createGuideSet} function) and conversion into Scis score (red dotted line) }
#'   \item{f}{ On-target binding profile for valid guides (rows) along target loci (columns) (sampled for large matrices). Colors indicate predicted binding affinity. }
#'   \item{g}{ Selection of guide representative (large dot) per cluster (colors). Black color represents group of blacklisted guides. }
#' }
#' @param guideSet guideSet containing combinations.
#' @return list of ggplot objects.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13A')
#' gs <- addGuides(gs, guide_length = 16)
#' gs <- plotGuides(gs)  
#' }
#' @seealso [addGuides()], and [export()]
#' @export
plotGuides <- function(guideSet)
{
  kmers <- as_tibble(guideSet@kmers) 
  kmer_stats <- .kmerStats(kmers)
  
  # gs@kmers %>% as_tibble %>% filter(on_target < 0) %>% ggplot(aes(Scis, Sbind, col = Soff)) + geom_point(alpha = 0.1) + scale_colour_gradientn(colours = colorRampPalette(c('darkblue', 'orange', 'darkred'))(256))
  
  p_guide_filt_counts <-
    kmers %>% 
    select(kmer_id, valid_gc, valid_score, valid) %>% 
    distinct %>% 
    summarise(all = n(), 
              pass_gc = sum(valid_gc), 
              pass_score = sum(valid_score),
              pass_all = sum(valid)) %>% 
    gather() %>% 
    mutate(key = factor(key, levels = c('all', 'pass_gc', 'pass_score', 'pass_all')[order(value, decreasing = TRUE)])) %>% 
    ggplot(aes(key, value)) + 
      geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle=90, hjust=1),
            axis.title.x = element_blank())      #+
      #geom_text(aes(label = value), position = position_dodge(width=0.9), vjust=-0.25)
  
  p_gc_hist <-
    kmers %>%
      select(kmer_id, gc, valid_gc) %>%
      distinct %>%
      #mutate(repname = forcats::fct_lump(repname, 5)) %>%
      ggplot(aes(gc, fill = valid_gc)) +
        geom_histogram(binwidth = 0.05) +
        xlim(0,1) +
        theme(axis.text.x = element_text(angle=90, hjust=1),
              legend.position = 'none')
        
  # p_score_hist <-
    # kmers %>% 
    # .kmerStats(., full = TRUE) %>% 
    # select(kmer_id, Son_tot, Soff_tot, valid_score) %>% 
    # distinct %>% 
    # gather('type', 'score', Son_tot, Soff_tot) %>% 
    # ggplot(aes(score, fill = valid_score)) + 
      # geom_histogram(position = 'identity', alpha = 0.5, bins = 50) + 
      # scale_x_log10() + 
      # facet_wrap(~type, nrow = 1, scales = 'free_x') + 
      # theme(legend.position = 'none')
  
  p_valid_score_gc <- # slow, improve!
    kmer_stats %>% 
    mutate(valid = kmer_id %in% 
                  (kmers %>% filter(valid) %>% pull(kmer_id))) %>%
    ggplot() +
      geom_point(aes(Soff_tot, Son_tot, col = valid), alpha = 0.05) +
      #geom_density_2d(aes(Soff_tot, Son_tot)) +
      scale_x_log10() +
      scale_y_log10() +
      #scale_color_manual(values = c('black', 'red')) +
      #scale_shape_manual(values = c(4, 19)) +
      #facet_grid(~valid_gc) +
      theme(legend.position = 'none')
      
  p_sbind_boxplot <- # samples to plot distribution
    kmers %>% 
      filter(valid) %>%
      filter(n_mismatches != 0) %>%
      sample_n(min(nrow(.), 10000), replace = FALSE) %>%
      ggplot(aes(x = as.factor(n_mismatches), y = Sbind, fill = as.factor(n_mismatches), group = n_mismatches)) + 
        #geom_violin() +
        geom_boxplot(width=.1, outlier.colour=NA) +
        theme(legend.position = 'none')
  
  # Plot cis decay - shows all cis hits regardless of Sbind or n_mismatches
  scaling_fact <- kmers %>% filter(on_target < 0 & cis_dist <= 100) %>% count %>% pull
  p_cis_decay <-
    kmers %>% 
      filter(on_target < 0) %>% #only off-targets
      filter(cis_dist <= 5000) %>%
      ggplot() + 
        geom_histogram(aes(x = cis_dist, stat(count / scaling_fact)), binwidth = 150) +
        geom_line(aes(x = cis_dist, y = Scis), col = 'red', lty = 2) + # cis decay function
        scale_y_continuous(limits = c(0, 5000), name = 'Scis', sec.axis = sec_axis(~.*scaling_fact, name = 'Count')) + 
        theme(legend.position = 'none')
  
  p_guide_hit_heatmap <- .plotClusts(guideSet)
  
  p_best_per_clust <- #slow 
    left_join(kmer_stats, 
              kmers %>% select(kmer_id, kmer_clust, best) %>% distinct) %>%
    ggplot(aes(Soff_tot, Son_tot, col = as.factor(kmer_clust))) +
      geom_point(aes(size = best, alpha = best)) +
      facet_wrap(~kmer_clust, ncol = 2) +
      #scale_color_manual(values = c('black', 'red')) +
      scale_x_log10() +
      scale_y_log10() +
      #scale_shape_manual(values = c(4, 19)) +
      scale_alpha_manual(values = c(0.1, 1)) +
      scale_size_manual(values = c(0.1, 2)) +
      ggplot2::guides(col = guide_legend(title = 'Cluster', ncol=2),
             shape = FALSE,
             size = FALSE) +
      theme(axis.text.x = element_text(angle=90, hjust=1),
            #axis.ticks.x = element_blank(),
            #axis.text.y = element_blank(),
            #axis.ticks.y = element_blank(),
            legend.position = 'none',
            panel.spacing = unit(0, 'lines'),
            strip.background = element_blank(),
            strip.text.x = element_blank())
    
  plots <- list('a' = p_guide_filt_counts, 
                'b' = p_gc_hist, 
                'c' = p_valid_score_gc, 
                'd' = p_sbind_boxplot, 
                'e' = p_cis_decay, 
                'f' = p_guide_hit_heatmap,
                'g' = p_best_per_clust) 
  
  p1 <- cowplot::plot_grid(p_guide_filt_counts,
                           p_gc_hist, 
                           p_valid_score_gc,
                           labels = 'auto', nrow = 1, rel_widths = c(1, 2, 2))
  # p2 <- cowplot::plot_grid(p_score_hist,
                           # p_valid_score_gc,
                           # labels = c('c', 'd'), ncol = 1)
  p2 <- cowplot::plot_grid(p_sbind_boxplot,
                           p_cis_decay,
                           labels = c('d', 'e'), nrow = 1, rel_widths = c(1, 2))
  
  p3 <- cowplot::plot_grid(p_guide_hit_heatmap,
                           p_best_per_clust,
                           labels = c('f', 'g'), ncol = 2, rel_widths = c(2, 1))
  
  plots_combined <- cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1, 1, 2.5))
  
  guideSet@plots$guides <- plots
  guideSet@plots$guides$figure <- plots_combined
  
  print(plots_combined)
  return(guideSet)
}

#' Generates combinations QC plots
#' 
#' Computes multiple QC plots and arranges into a figure with subpanels:
#' \describe{
#'   \item{a}{ Table showing total on/off-target scores (Son_tot and Soff_tot) for the best combination per number of guides (columns). }
#'   \item{b}{ Total on- and off-target score per number of greedy search iterations. Colors indicate the number of guides (from yellow (at least 2) to black (defined by \code{max_guides} in the \code{addCombinations} function). }
#'   \item{c}{ Total on- and off-target score per number of guides. }
#'   \item{d}{ Total on- and off-target score for all computated combinations per number of guides. }
#'   \item{e}{ Off-target coverage for the best guide combination per number of guides and off-target type (e.g. transposon family). Colors indicate predicted binding affinity to the off-target. }
#'   \item{f}{ Distance to the nearest cis regulatory feature (defined by \code{cis} in the \code{createGuideSet} function) for the best combination number of guides. }
#'   \item{g}{ On-target coverage for the best guide combination per number of guides and on-target transposon family. Colors indicate predicted binding affinity to the target. }
#'   \item{h}{ Coverage for the best guide combination on the consensus model of transposon families. Colors indicate the number of guides (from yellow (1) to black (\code{max_guides})) }
#' }
#' @param guideSet guideSet containing combinations.
#' @return list of ggplot objects.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13A')
#' gs <- plotTargets(gs)
#' }
#' @seealso [addCombinations()], and [export()]
#' @export
plotCombinations <- function(guideSet)
{ 
  if (length(guideSet@combinations) == 0) { stop ('Call addCombinations on guideSet before calling plotCombinations') }
  combinations <- guideSet@combinations 
  kmers <- 
    guideSet@kmers %>%
    as_tibble %>%
    right_join(., 
              combinations %>% filter(best) %>% unnest, by= 'kmer_id', 
              suffix = c('_kmer', '_combi')) %>% 
    mutate(Sbind_bin = cut(Sbind, breaks = c(0, 0.25, 0.5, 1), include.lowest = TRUE))
  
  n_fams <- length(guideSet@families)
  n_guides <- max(kmers$n_guides)
  n_guides_col <- colorRampPalette(c('gold', 'orange', 'darkred', 'black'))(n_guides)  

  p_table <-
    combinations %>%
    filter(best) %>%
    select(n_guides, Son_tot, Soff_tot, on_tot, off_tot) %>%
    mutate_all(round, digit = 1) %>%
    as.matrix %>%
    t() %>% 
    gridExtra::tableGrob(., theme = gridExtra::ttheme_default(base_size = 8))
  
  .custom_breaks <- function(x) { c(0, floor(max(x))) }
  p_on_off_bar <- 
    combinations %>% 
    filter(best) %>% 
    gather('type', 'counts', Son_tot, Soff_tot) %>% 
    ggplot(aes(n_guides, counts, fill = as.factor(n_guides))) + 
      geom_bar(stat = 'identity') +
      facet_wrap(~type, scales = 'free') +
      scale_fill_manual(values = n_guides_col) +
      scale_y_continuous(breaks = .custom_breaks) +
      theme(legend.position = 'none')
  
  p_score_scatter <-
    combinations %>%
    ggplot(aes(Soff_tot, Son_tot, col = best)) +
      geom_point(aes(size = best, alpha = best)) + 
      facet_wrap(~n_guides, nrow = 1) +
      scale_color_manual(values = c('black', 'red')) +
      scale_alpha_manual(values = c(0.1, 1)) +
      scale_size_manual(values = c(0.1, 2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = 'none',
            panel.spacing = unit(0.25, 'lines'))
            
  if (sum(grepl('greedy', combinations$combi_id)))
  {
    .custom_breaks_x <- function(x) { round(x) }
    .custom_breaks_y <- function(x) { c(ceiling(min(x)), floor(max(x))) }
    p_greedy_iterations_on <-
      combinations %>% 
      filter(!is.na(iterations)) %>% 
      ggplot(aes(iterations, Son_tot, col = as.factor(n_guides))) + 
        geom_line() + 
        geom_point(size = 0.5) +
        scale_x_continuous(breaks = .custom_breaks_x) +
        #scale_y_continuous(breaks = .custom_breaks_y) +
        scale_color_manual(values = n_guides_col) +
        theme(legend.position = 'none')
    p_greedy_iterations_off <-
      combinations %>% 
      filter(!is.na(iterations)) %>% 
      ggplot(aes(iterations, Soff_tot, col = as.factor(n_guides))) + 
        geom_line() + 
        geom_point(size = 0.5) +
        scale_x_continuous(breaks = .custom_breaks_x) +
        #scale_y_continuous(breaks = .custom_breaks_y) +
        scale_color_manual(values = n_guides_col) +
        theme(legend.position = 'none')
  } else { 
    p_greedy_iterations_on <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = 'No greedy optimization') 
    p_greedy_iterations_off <- p_greedy_iterations_on
  }
  p_greedy_iterations <- cowplot::plot_grid(p_greedy_iterations_on, p_greedy_iterations_off)
  
  .custom_breaks <- function(x) { c(0, floor(max(x))) }
  p_offtarget_cov <- 
    kmers %>% 
    filter(on_target < 0) %>% 
    group_by(n_guides, unique_id) %>% 
      summarise(Sbind = sum(Sbind), repname = unique(repname)) %>% 
    ungroup %>% 
    mutate(Sbind_bin = cut(Sbind, breaks = c(0, 0.5, 1, 2, 4, Inf), include.lowest = TRUE)) %>% 
    count(n_guides, repname, Sbind_bin) %>% 
    mutate(repname = forcats::fct_lump(repname, n = 4, ties.method = 'first')) %>%
    ggplot(aes(n_guides, n, fill = Sbind_bin)) + 
      geom_bar(stat = 'identity') + 
      scale_y_continuous(breaks = .custom_breaks) +
      facet_wrap(~repname, nrow = 1, scales = 'free_y') +
      theme(#axis.text.x = element_text(angle = 90),
            #panel.spacing = unit(0.25, 'lines'))
            legend.position = 'right',
            legend.key.width = unit(.3,"cm"),
            legend.key.height = unit(.3,"cm")) +
      ggplot2::guides(fill = guide_legend(title = 'Sbind', ncol = 1))
  
  # Plot offtarget dist histogram
  ggdata <- 
    kmers %>% 
    filter(on_target < 0) %>% 
    filter(cis_dist < 5000)
  if (nrow(ggdata) > 0)
  {
    p_offtarget_dist <-
      ggdata %>%
      ggplot(aes(cis_dist, fill = as.factor(n_guides))) + 
        geom_histogram(binwidth = 250) + 
        facet_wrap(~n_guides, nrow = 1) + 
        xlim(0, 5000) +
        scale_fill_manual(values = n_guides_col) +
        theme(legend.position = 'none',
              axis.text.x = element_text(angle = 90))
  } else {
    p_offtarget_dist <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = 'No cis proximal off targets')   
  }
  
  p_target_cov <- 
    kmers %>% 
    filter(on_target > 0) %>% 
    group_by(n_guides, te_id) %>% 
      summarise(Son = sum(Son), repname = unique(repname)) %>% 
    ungroup %>% 
    mutate(Son_bin = cut(Son, breaks = c(0, 0.5, 1, 2, 4, Inf), include.lowest = TRUE)) %>% 
    count(n_guides, repname, Son_bin) %>% 
    ggplot(aes(n_guides, n, fill = Son_bin)) + 
      geom_bar(stat = 'identity') + 
      facet_wrap(~repname, nrow = 1, scales = 'free_y') +
      theme(panel.spacing = unit(0.25, 'lines'),
            legend.position = 'right',
            legend.key.width = unit(.3,"cm"),
            legend.key.height = unit(.3,"cm")) +
      ggplot2::guides(fill = guide_legend(title = 'Sbind', ncol = 1))

  # p_cis_hist <-   
    # kmers %>%
    # filter(on_target < 0) %>%
    # filter(cis_dist <= 5000) %>%
    # group_by(n_guides, unique_id) %>% 
      # summarise(Sbind = sum(Sbind), cis_dist = mean(cis_dist)) %>% 
    # ungroup %>%
    # mutate(Sbind_bin = cut(Sbind, breaks = c(0, 0.5, 1, 2, Inf), include.lowest = TRUE)) %>%
    # ggplot(aes(cis_dist, fill = Sbind_bin)) +
      # geom_histogram(binwidth = 500) +
      # facet_wrap(~n_guides, nrow = 1) +
      # theme(axis.text.x = element_text(angle = 90, hjust = 1),
            # panel.spacing = unit(0.25, 'lines'),
            # strip.background = element_blank(),
            # strip.text.x = element_blank())
    
  if (length(guideSet@alignments) > 0) 
  {
    .custom_breaks <- function(x) { c(0, floor(max(x))) }
    p_cons_cov <-
      kmers %>% 
      filter(on_target == 1) %>% 
      mutate(con_bin = cut(con_pos, breaks = seq(0, max(con_pos) + 250, 250), include.lowest = TRUE, dig.lab = 5)) %>%
      count(repname, n_guides, con_bin) %>%
      ggplot(aes(con_bin, n, fill = as.factor(n_guides))) + 
        geom_bar(stat = 'identity', position = 'dodge') + 
        facet_wrap(~repname, scales = 'free_y', ncol = 2, strip.position = 'left') +
        #scale_y_continuous(breaks = .custom_breaks) +
        scale_fill_manual(values = n_guides_col) +
        scale_y_continuous(breaks = .custom_breaks) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
              axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = 'none')
              # strip.background = element_blank(),
              # strip.text.x = element_blank())
  } else { 
    p_cons_cov <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = 'No alignment provided') 
  }
  
  # p_heatmap <- # add clustering!
    # kmers %>% 
    # group_by(n_guides, kmer_id, te_id) %>% 
      # summarise(Son = sum(Son)) %>% 
    # ungroup %>% 
    # ggplot(aes(as.factor(kmer_id), as.factor(te_id), fill = Son)) + 
      # ggrastr::geom_tile_rast() + 
      # facet_wrap(~n_guides, nrow = 1) +
      # theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
            # axis.text.x = element_text(angle = 90, hjust = 1),
            # panel.grid.major = element_blank())
        
  
  plots <- list('a' = p_table, 
                'b' = p_greedy_iterations,
                'c' = p_on_off_bar, 
                'd' = p_score_scatter, 
                'e' = p_offtarget_cov, 
                'f' = p_offtarget_dist, 
                'g' = p_target_cov,
                'h' = p_cons_cov) 
                
  p1 <- cowplot::plot_grid(p_table, labels = 'a')
  p2 <- cowplot::plot_grid(p_greedy_iterations,
                           p_on_off_bar,
                           labels = c('b', 'c'), nrow = 1)
  p3 <- cowplot::plot_grid(p_score_scatter,  
                           p_offtarget_cov,
                           p_offtarget_dist,
                           p_target_cov,
                           labels = c('d', 'e', 'f', 'g'), ncol = 1)
  p4 <- cowplot::plot_grid(p_cons_cov,
                           labels = 'h')
                           
  plots_combined <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1, rel_heights = c(1, 1, 4, 2))
  guideSet@plots$combinations <- plots
  guideSet@plots$combinations$figure <- plots_combined
  
  print(plots_combined)
  return(guideSet)
}