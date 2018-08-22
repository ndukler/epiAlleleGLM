methods::setGeneric("getLabelIndices", function(obj,query) {
  standardGeneric("getLabelIndices")
})

#' Get all site indices for a label 
#'
#' Gets indicies of all sites with a given label
#' @param obj a rateModel object
#' @param query a label
#' @name getLabelIndices
#' @return a numeric vector of indices
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("getLabelIndices", signature(obj = "rateModel"), function(obj,query) {
  return(obj@siteLabels[.(query)]$index)
})