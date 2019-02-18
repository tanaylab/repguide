#' Essential human gene promoter coordinates
#'
#' Genomic coordinates (hg38) of essential human gene promoters. Transcription start sites of essential genes (Hart et al., 2015) were determined from Gencode annotation (v28) and extended 1 kb up- and downstream.
#'
#' @format bed file with 13944 genomic intervals (each 2 kb in size).
#' @name essential_genes
#'
#' @source [Gencode v28 annotation](https://www.gencodegenes.org/human/), and essential genes [Hart et al., 2015](http://tko.ccbr.utoronto.ca/Data/core-essential-genes-sym_HGNCID)
#' @section hg38_essential_gene_tss_coords.bed:
#'
#' This data is used to exemplify \code{blacklist} of the \code{createGuideSet} function.
NULL