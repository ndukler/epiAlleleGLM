methods::setGeneric("getSiteInfo", function(obj) {
  standardGeneric("getSiteInfo")
})

#' Extract site information table
#'
#' Extracts site information table
#' @param obj alleleData or rateModel object
#' @name getSiteInfo
#' @return data.table containing site info
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("getSiteInfo", signature(obj = "alleleData"), function(obj) {
  return(obj@siteInfo)
})
methods::setMethod("getSiteInfo", signature(obj = "rateModel"), function(obj) {
  return(getSiteInfo(getAlleleData(obj)))
})