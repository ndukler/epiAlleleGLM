#' Get rate design matrix
#'
#' Get rate design matrix
#' @param obj rateModel
#' @include rateModel-class.R
#' @rdname getRateDM
#' @examples
#' 
#' @export
methods::setGeneric("getRateDM", function(obj,...) {
  standardGeneric("getRateDM")
})

#' @name getRateDM
#' @rdname getRateDM
methods::setMethod("getRateDM", signature(obj = "rateModel"), function(obj) {
  return(obj@rateDM)
})