% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute.R
\name{clustGuides}
\alias{clustGuides}
\title{Cluster guideRNAs}
\usage{
clustGuides(guideSet, min_Son = 0, n_clust = 15, alpha = 10)
}
\arguments{
\item{guideSet}{guideSet containing guide mappings}

\item{min_Son}{Numeric from 0 through 1. Only considers genomic target binding sites above \code{min_Son} score.}

\item{n_clust}{Integer from 1 to 20. Number of clusters to group guides into. Passed to \code{cutree()} function.}
}
\value{
guideSet object with clustered guides.
}
\description{
Clusters guides based on their target binding profile
}
\examples{
\dontrun{
gs <- createGuideSet(Hsapiens, tes = te_annotation_df)
gs <- addTargets(gs, targets = 'LTR13A')
gs <- addGuides(gs, guide_length = 16)
gs <- compClusts(gs, min_Son = 0.25, n_clust = 10)
gs <- plotGuides(gs)
}
}
\seealso{
\code{\link[=addGuides]{addGuides()}}, and \code{\link[=plotGuides]{plotGuides()}}
}
