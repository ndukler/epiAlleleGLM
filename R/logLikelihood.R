#' Compute log-likelihood
#'
#' Compute log-likelihood for given rate model object
#' @param x a set of parameters to compute the log-likelihood for
#' @param obj rateModel object
#' @param stickParams
#' @return A rateModel object with updated parameter values
#' @name logLikelihood
#' @include rateModel-class.R
#' @rdname logLikelihood
#' @examples
#' 
#' @export
methods::setGeneric("logLikelihood", function(x,obj,...) {
  standardGeneric("logLikelihood")
})

#' @name logLikelihood
#' @rdname logLikelihood
methods::setMethod("logLikelihood", signature(x="missing",obj = "rateModel"), function(x,obj) {
  siteTypes=levels(obj@siteLabels$siteLabel)
  nTips=length(getTree(obj)$tip.label)
  # l = siteTypes[1]
  `%myPar%` <- ifelse(foreach::getDoParRegistered(), yes = foreach::`%dopar%`, no = foreach::`%do%`)
  siteGroupLik=foreach::foreach(l=siteTypes) %myPar% {
    ## Extract subset of data for site
    data=getAlleleData(obj)@data[getLabelIndices(obj,l),,drop=FALSE]
    ## Create traversal table and rates
    tt=data.table::data.table(getTree(obj)$edge)
    data.table::setnames(x = tt,old=colnames(tt),c("parent","child"))
    rates=getParams(obj)[getRateIndex(obj,edges = tt,siteLabel = l)]
    pi=getParams(obj)[getPiIndex(obj,siteLabel = l)]
    ## Compute transition matricies
    logTransMat=branchRateMatrix(rate = rates,branch.length =  getTree(obj)$edge.length,pi = pi)
    siteLik=treeLL(data=data,tMat=logTransMat,traversal=as.matrix(tt-1),nTips=nTips,logPi=log(pi))
    return(sum(siteLik))
  }
  ll=sum(unlist(siteGroupLik))
  return(ll)
})

#' @name logLikelihood
#' @rdname logLikelihood
methods::setMethod("logLikelihood", signature(x="numeric",obj = "rateModel"), function(x,obj,stickParams=TRUE) {
  if(!stickParams){
    if(length(x) != length(getParams(obj))){
      stop("Length of x must be equal to the number of parameters in the model if stickParams is not TRUE.")
    }
    y=x
  } else{
    xPiLen=length(x)-max(getParamIndex(obj)$rateIndex)
    expPiLen = length(levels(obj@siteLabels$siteLabel))*(getAlleleData(obj)@nAlleles-1)
    ## Check that the vector of stick breaking parameters are of the right length and 
    ## if so, convert them to probabilities
    if(xPiLen==expPiLen){
      probPi=epiAllele:::multiStickToProb(x[(max(getParamIndex(obj)$rateIndex)+1):length(x)],
                                          getAlleleData(obj)@nAlleles-1)
      y=c(x[1:max(getParamIndex(obj)$rateIndex)],probPi)
    } else{
      stop("Incorrect number of parameters for stick breaking process sent to logLikelihood function")
    }
  }
  ## Set parameter values without object duplication
  setParamValue(obj=obj,i = 1:length(y),val = y)
  ## Pass on to logLikelihood method that doesn't specify parameters
  return(logLikelihood(obj=obj))
})