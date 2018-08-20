methods::setGeneric("clone", function(object) {
  standardGeneric("clone")
})

#' Copy object
#'
#' Creates a copy of a rateModel object
#' @param object rateModel 
#' @name clone
#' @return rateModel object
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("clone", signature(object = "rateModel"), function(object) {
  ## Create deep copy of environment
  paramEnviron = new.env()
  paramEnviron$paramIndex=object@paramEnviron$paramIndex
  paramEnviron$params=object@paramEnviron$params
  paramEnviron$fixed=object@paramEnviron$fixed
  ## Return copy of object
  return(methods::new("rateModel",alleleData=object@alleleData,edgeGroups=object@edgeGroups,siteLabelCriteria=object@siteLabelCriteria,
               siteLabels=object@siteLabels,paramEnviron=paramEnviron))
})

