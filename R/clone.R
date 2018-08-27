methods::setGeneric("clone", function(object,deep=TRUE) {
  standardGeneric("clone")
})

#' Copy object
#'
#' Creates a copy of a rateModel object
#' @param object rateModel
#' @param deep Logical value as to whether to make a deep copy of the data (default = TRUE)
#' @name clone
#' @return rateModel object
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("clone", signature(object = "rateModel",deep="logical"), function(object,deep=TRUE) {
  if(deep){
    ## Create deep copy of environment
    adEnviron = new.env()
    adEnviron$alleleData=getAlleleData(object)
  } else {
    adEnviron=object@alleleData
  }
  ## Deep vector copies
  p=numeric(length(object@params))
  p[1:length(p)]=object@params[1:length(p)]
  f=logical(length(object@fixed))
  f[1:length(f)]=object@fixed[1:length(f)]
  ## Return copy of object
  return(methods::new("rateModel",alleleData=adEnviron,edgeGroups=data.table::copy(object@edgeGroups),
                      siteLabelCriteria=data.table::copy(object@siteLabelCriteria),siteLabels=data.table::copy(object@siteLabels),
                      params=p,paramIndex=data.table::copy(object@paramIndex),fixed=f))
})

