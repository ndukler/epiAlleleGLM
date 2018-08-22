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
  adEnviron = new.env()
  adEnviron$alleleData=getAlleleData(object)
  ## Return copy of object
  return(methods::new("rateModel",alleleData=adEnviron,edgeGroups=object@edgeGroups,siteLabelCriteria=object@siteLabelCriteria,
               siteLabels=object@siteLabels,params=object@params,paramIndex=object@paramIndex,fixed=object@fixed))
})

