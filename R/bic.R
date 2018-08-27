#' Compute BIC for rate models
#'
#' Computes the bayesoan information criterion (BIC) for a given rateModel.
#' Note that when comparing two models, the one with the smaller BIC is preferred. 
#' @param obj rateModel object
#' @name bic
#' @include rateModel-class.R
#' @rdname bic
#' @return numeric BIC value
#' @examples
#' 
#' @export
methods::setGeneric("bic", function(obj,i) {
  standardGeneric("bic")
})

#' @name bic
#' @rdname bic
methods::setMethod("bic", signature(obj = "rateModel"), function(obj) {
  ll=logLikelihood(obj=obj)
  k=sum(obj@fixed)
  return(log(nSites(obj))*k-2*ll)
})