#' Get parameter values
#'
#' Get parameter values given a set of indices
#' @param obj rateModel object
#' @param i a numeric vector of indices
#' @name getParamValue
#' @include rateModel-class.R
#' @rdname getParamValue
#' @examples
#' 
#' @export
methods::setGeneric("getParamValue", function(obj,i) {
  standardGeneric("getParamValue")
})

#' @name getParamValue
#' @rdname getParamValue
methods::setMethod("getParamValue", signature(obj = "rateModel"), function(obj,i) {
  if(!is.numeric(i)){
    stop("Indices must be numeric")
  }
  return(obj@params[i])
})