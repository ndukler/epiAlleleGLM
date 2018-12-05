#' Get pi parameter index
#'
#' Returns indecies of the stationary distribution parameters associated with a given siteLabel 
#' @param obj rateModel object
#' @param allele The column index of the corresponding allele
#' @name getPiIndex
#' @include rateModel-class.R
#' @rdname getPiIndex
#' @return A numeric vector of parameter indices
#' @examples
#' 
#' @export
methods::setGeneric("getPiIndex", function(obj,allele) {
  standardGeneric("getPiIndex")
})

#' @name getPiIndex
#' @rdname getPiIndex
methods::setMethod("getPiIndex", signature(obj = "rateModel"), function(obj,allele) {
  if(allele<=1){
    stop("Allele values must be >1")
  }
  al=allele
  ind=obj@piIndex[.(al)]$index
  return(ind)
})