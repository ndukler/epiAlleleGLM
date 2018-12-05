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
  nTips=length(getTree(obj)$tip.label)
  nSites=getAlleleData(obj)@nSites
  ## Create traversal table and rates
  tt=data.table::data.table(getTree(obj)$edge)
  data.table::setnames(x = tt,old=colnames(tt),c("parent","child"))
  ## Map edges to edge groups
  ttAug=getEdgeGroupTable(obj)[with(tt,paste0(parent,"-",child))]
  ## Iterate over unique edgeGroups and compute rates for them
  rateMatrix=Matrix::Matrix(1,nrow=nSites,ncol=nrow(tt))
  for(e in unique(ttAug$edgeGroup)){
      rateMatrix[,e==ttAug$edgeGroup]=epiAllele:::computeRates(obj,e)
  }
  pi=epiAllele:::computePi(obj)
  `%myPar%` <- ifelse(foreach::getDoParRegistered(), yes = foreach::`%dopar%`, no = foreach::`%do%`)
  siteLik=foreach::foreach(i=1:nSites) %myPar% {
    ## Compute transition matricies
    ltm=branchRateMatrix(rate = rateMatrix[i,],branch.length =  getTree(obj)$edge.length,pi = pi[i,])
    ## Re-sort logTransMat so that the matricies are ordered in the list by the node number of the child
    logTransMat=list()
    logTransMat[tt$child]=ltm
    logTransMat[[tt$parent[nrow(tt)]]]=matrix(0,length(pi),length(pi)) ## placeholder matrix to avoid error when passing to Rcpp
    ## Compute log-likelihood
    siteLik=treeLL(data=data,tMat=logTransMat,traversal=as.matrix(tt-1),nTips=nTips,logPi=log(pi))
    return(siteLik)
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