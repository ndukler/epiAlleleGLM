#' Extract site information table
#'
#' Extracts site information table
#' @param obj alleleData or rateModel object
#' @name nSites
#' @return number of sites in data
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setGeneric("nSites", function(obj) {
  standardGeneric("nSites")
})

#' @name nSites
#' @rdname nSites
methods::setMethod("nSites", signature(obj = "alleleData"), function(obj) {
  return(obj@nSites)
})
methods::setMethod("nSites", signature(obj = "rateModel"), function(obj) {
  return(nSites(getAlleleData(obj)))
})