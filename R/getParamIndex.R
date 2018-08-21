#' Get parameter index
#'
#' Returns parameter index as data.table 
#' @param obj rateModel object
#' @name getParamIndex
#' @include rateModel-class.R
#' @rdname getParmIndex
#' @return A data.table of parameter indices, keyed by edgeGroups and siteLabels
#' @examples
#' 
#' @export
methods::setGeneric("getParamIndex", function(obj) {
  standardGeneric("getParamIndex")
})

#' @name getParamIndex
#' @rdname getParamIndex
methods::setMethod("getParamIndex", signature(obj = "rateModel"), function(obj) {
  return(obj@paramEnviron$paramIndex)
})