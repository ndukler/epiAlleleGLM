#' Set parameter values
#'
#' Get parameter values given a set of indices
#' @param obj rateModel object
#' @param i a numeric vector of indices to set to \'value\' argument
#' @param value a numeric vector of values to set 
#' @return A rateModel object with updated parameter values
#' @name setParamValue
#' @include rateModel-class.R
#' @rdname setParamValue
#' @examples
methods::setGeneric("setParamValue", function(obj,i,value) {
  standardGeneric("setParamValue")
})

#' @name setParamValue
#' @rdname setParamValue
methods::setMethod("setParamValue", signature(obj = "rateModel"), function(obj,i,value) {
  if(!is.numeric(i)){
    stop("Indices must be numeric")
  }
  if(length(i)!=length(value)){
    stop("There must be the same number of values as indices")
  }
  setValues(obj@params,i,value) ## use rcpp function to update object without causing duplication
})