#' Get all site indices for a label 
#'
#' Gets indicies of all sites with a given label
#' @param obj a rateModel object
#' @param label a site label
#' @name getLabelIndices
#' @return a numeric vector of indices
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setGeneric("getLabelIndices", function(obj,label) {
  standardGeneric("getLabelIndices")
})

#' @name getLabelIndices
#' @rdname getLabelIndices
methods::setMethod("getLabelIndices", signature(obj = "rateModel"), function(obj,label) {
  resp=obj@siteLabels[.(label),on="siteLabel",nomatch=0]$index
  if(length(resp)==0){
    warning("No matching label to query")
  }
  return(resp)
})