#' Compute log-likelihood
#'
#' Compute log-likelihood for given rate model object
#' @param obj rateModel object
#' @return A rateModel object with updated parameter values
#' @name logLikelihood
#' @include rateModel-class.R
#' @rdname logLikelihood
#' @examples
#' 
#' @export
methods::setGeneric("logLikelihood", function(obj,threads) {
  standardGeneric("logLikelihood")
})

#' @name logLikelihood
#' @rdname logLikelihood
methods::setMethod("logLikelihood", signature(obj = "rateModel"), function(obj,threads=1) {
  if(!is.numeric(threads) || length(threads) > 1){
    stop("Threads must be a single numeric value")
  }
  
  
})