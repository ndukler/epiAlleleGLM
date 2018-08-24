#' Compute scaled negative log-likelihood
#'
#' Compute scaled negative log-likelihood for given rate model object
#' @param x a set of parameters to compute the log-likelihood for
#' @param obj rateModel object
#' @param stickParams
#' @return A rateModel object with updated parameter values
#' @name scaledLL
#' @include rateModel-class.R
#' @rdname scaledLL
methods::setGeneric("scaledLL", function(x,obj,scale,...) {
  standardGeneric("scaledLL")
})

#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="missing",obj = "rateModel"), function(x,obj,scale) {
  ll=logLikelihood(obj = obj)
  return(ll)
})

#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="numeric",obj = "rateModel"), function(x,obj,scale,stickParams=TRUE) {
  ll=logLikelihood(x = x,obj = obj,stickParams = stickParams)
  return(ll*scale)
})