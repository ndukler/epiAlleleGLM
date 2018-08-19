methods::setGeneric("getSiteInfo", function(obj) {
  standardGeneric("getSiteInfo")
})

#' Extract site information table
#'
#' Extracts site information table
#' @param x alleleData or rateModel object
#' @name getEdgeTable
#' @return data.table containing edge ids
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("getSiteInfo", signature(obj = "alleleData"), function(obj) {
  return(obj@siteInfo)
})
methods::setMethod("getSiteInfo", signature(obj = "rateModel"), function(obj) {
  return(obj@data@siteInfo)
})