#' Compute pi
#'
#' Compute per site stationary distribution per allele
#' 
#' @param obj rateModel
#' @include rateModel-class.R
#' @rdname computePi
#' @examples
#' 
#' @export
methods::setGeneric("computePi", function(obj) {
  standardGeneric("computePi")
})

#' @name computePi
#' @rdname computePi
methods::setMethod("computePi", signature(obj = "rateModel"), function(obj) {
  co=list()
  ## Get the coefficients for all alleles
  for(a in 2:getAlleleData(obj)@nAlleles){
    co[[as.character(a)]]=getParams(obj)[getPiIndex(obj,a)]
  }
  ## Place coefficients in coefficient x allele matrix
  coMat=do.call("cbind",co)
  ## Compute probabilities for each allele
  potential=as.matrix(cbind(1,exp(getPiDM(obj) %*% coMat)))
  return(potential/rowSums(potential))
})
