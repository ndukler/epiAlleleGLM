#' Get pi parameter index
#'
#' Returns indecies of the stationary distribution parameters associated with a given siteLabel 
#' @param obj rateModel object
#' @param siteLabel A character vector of length one for one of the site labels
#' @name getPiIndex
#' @include rateModel-class.R
#' @rdname getPiIndex
#' @return A numeric vector of parameter indices
#' @examples
#' 
#' @export
methods::setGeneric("getPiIndex", function(obj,siteLabel) {
  standardGeneric("getPiIndex")
})

#' @name getPiIndex
#' @rdname getPiIndex
methods::setMethod("getPiIndex", signature(obj = "rateModel"), function(obj,siteLabel) {
  sl=siteLabel
  ind=getParamIndex(obj)[.(sl),on="siteLabel"]$piIndex[1]
  return(ind:(ind+getAlleleData(obj)@nAlleles-1))
})