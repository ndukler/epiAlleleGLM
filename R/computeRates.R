#' Compute rates
#'
#' Compute per site rates from regression model
#' Uses rate = exp(X %*% B) where X is the design matrix and B are the coefficients
#' @param obj rateModel
#' @param eg the edgeGroup to compute rates for
#' @include rateModel-class.R
#' @rdname computeRates
#' @examples
#' 
#' @export
methods::setGeneric("computeRates", function(obj,...) {
  standardGeneric("computeRates")
})

#' @name computeRates
#' @rdname computeRates
methods::setMethod("computeRates", signature(obj = "rateModel"), function(obj,eg) {
  return(as.matrix(xp(getRateDM(obj) %*% getParams(obj)[getRateIndex(obj = obj,edgeGroup = eg)],ncol=1)))
})