#' Get all site labels for an index 
#'
#' Gets label of all sites with a given index
#' @param obj a rateModel object
#' @param i a vector of site indices
#' @name getIndexLabels
#' @return a factor vector of labels
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setGeneric("getIndexLabels", function(obj,i) {
  standardGeneric("getIndexLabels")
})

#' @name getIndexLabels
#' @rdname getIndexLabels
methods::setMethod("getIndexLabels", signature(obj = "rateModel"), function(obj,i) {
  resp=obj@siteLabels[i,nomatch=0]$siteLabel
  if(length(resp)!=length(i)){
    warning("At least one out of range query index")
  }
  return(resp)
})