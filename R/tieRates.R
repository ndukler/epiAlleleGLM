methods::setGeneric("tieRates", function(obj,siteLabelA,edgeGroupA,siteLabelB,edgeGroupB) {
  standardGeneric("tieRates")
})

#' Get all site indices for a label 
#'
#' Gets indicies of all sites with a given label
#' @param obj a rateModel object
#' @param siteLabelA a label
#' @param edgeGroupA foo
#' @param siteLabelB a label
#' @param edgeGroupB foo
#' @name tieRates
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("tieRates", signature(obj = "rateModel"), function(obj,siteLabelA,edgeGroupA,siteLabelB,edgeGroupB) {
  return(NULL)
})