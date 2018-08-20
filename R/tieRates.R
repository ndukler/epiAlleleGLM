methods::setGeneric("tieRates", function(rateMod,siteLabelA,edgeGroupA,siteLabelB,edgeGroupB) {
  standardGeneric("tieRates")
})

#' Get all site indices for a label 
#'
#' Gets indicies of all sites with a given label
#' @param rateMod a rateModel object
#' @param siteLabelA 
#' @param edgeGroupA foo
#' @param siteLabelB a label
#' @param edgeGroupB foo
#' @name tieRates
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("tieRates", signature(rateMod = "rateModel"), function(rateMod,siteLabelA,edgeGroupA,siteLabelB,edgeGroupB) {
  
})