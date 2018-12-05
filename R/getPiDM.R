#' Get pi design matrix
#'
#' Get pi design matrix
#' @param obj rateModel
#' @include rateModel-class.R
#' @rdname getPiDM
#' @examples
#' 
#' @export
methods::setGeneric("getPiDM", function(obj,...) {
  standardGeneric("getPiDM")
})

#' @name getPiDM
#' @rdname getPiDM
methods::setMethod("getPiDM", signature(obj = "rateModel"), function(obj) {
  return(obj@piDM)
})