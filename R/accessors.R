#' @include guideSet.R
NULL

#' Returns genomic mappings from guideSet
#'
#' Retrieves data.frame with guide's genomic mappings and annotation.
#' @export
setGeneric("mappings", function(guideSet) standardGeneric("mappings"))
setMethod("mappings", signature("guideSet"), function(guideSet) {
  out <- guideSet@kmers
  return(out)
})

#' Returns transposable element families from guideSet
#'
#' Retrieves all available families in the provided guideSet that match defined characteristics.
#'
#' @param guideSet guideSet object to query.
#' @param pattern Regular expression. Returns only families that match \code{pattern}.
#' @param fixed Logical. If \code{TRUE}, \code{pattern} is matched as is.
#' @return Character vector.
#' @examples
#' \dontrun{
#' gs <- createGuideSet(genome = Hsapiens, tes = te_annotation_df)
#' families <- repnames(gs, pattern = 'LTR12') # returns all families containing 'LTR12' in their name.
#'}
#' @seealso [createGuideSet()]
#' @export
repnames <- function(guideSet, 
                     pattern = NULL,
                     fixed = FALSE)
{
  if (length(pattern) > 1) { stop ('Pattern must be of length 1') }
  families <- unique(guideSet@tes$repname)
  if (!is.null(pattern))
  {
    families <- grep(pattern, families, fixed = fixed, value = TRUE)
  }
  
  if (length(families) == 0)
  {
    families <- unique(guideSet@tes$repname)
    fam_dists <- structure(stringdist::stringdist(pattern, unique(guideSet@tes$repname)), names = families)
    suggestions <- paste(names(head(sort(fam_dists), 5)), collapse = ', ')
    stop (glue::glue("{pattern} not found in guideSet. Did you mean any of {suggestions}?"))
  }
  
  return(families)
}

#' Returns blacklisted regions from guideSet
#' @export
setGeneric("blacklisted", function(guideSet) standardGeneric("blacklisted"))
setMethod("blacklisted", signature("guideSet"), function(guideSet) {
  out <- guideSet@blacklist
  return(out)
})

#' Returns blacklisted regions from guideSet
#' @export
setGeneric("combinations", function(guideSet) standardGeneric("combinations"))
setMethod("combinations", signature("guideSet"), function(guideSet) {
  out <- guideSet@combinations
  return(out)
})