#' Compute log-likelihood
#'
#' Compute log-likelihood for given rate model object
#' @param x a set of parameters to compute the log-likelihood for
#' @param obj rateModel object
#' @return A rateModel object with updated parameter values
#' @name logLikelihood
#' @include rateModel-class.R
#' @rdname logLikelihood
#' @examples
#' 
#' @export
methods::setGeneric("logLikelihood", function(x,obj) {
  standardGeneric("logLikelihood")
})

# obj=rateMod

#' @name logLikelihood
#' @rdname logLikelihood
methods::setMethod("logLikelihood", signature(x="missing",obj = "rateModel"), function(x,obj) {
  siteTypes=levels(obj@siteLabels$siteLabel)
  nTips=length(getTree(obj)$tip.label)
  # l = siteTypes[1]
  siteGroupLik=foreach::foreach(l=siteTypes) %dopar% {
    ## Extract subset of data for site
    data=getAlleleData(obj)@data[getLabelIndices(obj,l),,drop=FALSE]
    ## Create traversal table and rates
    tt=data.table::data.table(getTree(obj)$edge)
    data.table::setnames(x = tt,old=colnames(tt),c("parent","child"))
    rates=getParams(obj)[getRateIndex(obj,edges = tt,siteLabel = l)]
    pi=getParams(obj)[getPiIndex(obj,siteLabel = l)]
    ## Compute transition matricies
    logTransMat=epiAllele:::branchRateMatrix(rate = rates,branch.length =  getTree(obj)$edge.length,pi = pi)
    siteLik=treeLL(data=data,tMat=logTransMat,traversal=as.matrix(tt-1),nTips=nTips,logPi=log(pi))
    return(sum(siteLik))
  }
  return(sum(unlist(siteGroupLik)))
})

#' @name logLikelihood
#' @rdname logLikelihood
methods::setMethod("logLikelihood", signature(x="numeric",obj = "rateModel"), function(x,obj) {
  if(length(x) != length(getParams(obj))){
    stop("Length of x must be equal to the number of parameters in the model.")
  }
  ## Set parameter values without object duplication
  setValues(x,1:length(x),obj)
  ## Pass on to logLikelihood method that doesn't specify parameters
  logLikelihood(obj)
})