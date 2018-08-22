methods::setGeneric("getParams", function(obj) {
  standardGeneric("getParams")
})

#' Returns alleleData
#'
#' Gets the alleleData from a rate model object
#' @param obj rateModel 
#' @name getParams
#' @return alleleData object
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("getParams", signature(obj = "rateModel"), function(obj) {
  return(obj@params)
})

