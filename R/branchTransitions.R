#' Compute transitions per branch
#'
#' Compute the number of expected transitions per branch between states based on the marginal joint probabilities
#' @param obj rateModel object
#' @param subset The site indices for which should be summed over to compute the total number of expected transitions (NOT IMPLEMENTED)
#' @param ncores Number of cores to use
#' @return A 3-dimensional vectors where the first index is the branch, the second is the state of the parent, and the third 
#' is the state of the child. The entries are the exprected number of transitions on that branch for the corresponding parent
#' and child states.
#' @name branchTransitions
#' @include rateModel-class.R
#' @rdname branchTransitions
#' @examples
#' 
#' @export
methods::setGeneric("branchTransitions", function(obj,subset=NULL,ncores=1) {
  standardGeneric("branchTransitions")
})

#' @name branchTransitions
#' @rdname branchTransitions
methods::setMethod("branchTransitions", signature(obj = "rateModel"), function(obj,subset=NULL,ncores=1) {
  siteTypes=levels(obj@siteLabels$siteLabel)
  nTips=length(getTree(obj)$tip.label)
  ## Create a list of all siblings for each node 
  siblingList=list()
  for(n in sort(unique(as.numeric(getTree(obj)$edge)))){
    siblingList[[n]]=getSiblings(getTree(obj),n)-1
  }
  #  l = siteTypes[1]
  ## `%myPar%` <- ifelse(foreach::getDoParRegistered(), yes = foreach::`%dopar%`, no = foreach::`%do%`)
  groupTrans=foreach::foreach(l=siteTypes,.inorder = TRUE, .final = function(x) setNames(x, siteTypes)) %do% {
    ## Extract subset of data for site
    data=getAlleleData(obj)@data[getLabelIndices(obj,l),,drop=FALSE]
    ## Create traversal table and rates
    tt=data.table::data.table(getTree(obj)$edge)
    data.table::setnames(x = tt,old=colnames(tt),c("parent","child"))
    rates=getParams(obj)[getRateIndex(obj,edges = tt,siteLabel = l)]
    pi=getParams(obj)[getPiIndex(obj,siteLabel = l)]
    ## Compute transition matricies
    ltm=branchRateMatrix(rate = rates,branch.length =  getTree(obj)$edge.length,pi = pi)
    ## Re-sort logTransMat so that the matricies are ordered in the list by the node number of the child
    logTransMat=list()
    logTransMat[tt$child]=ltm
    logTransMat[[tt$parent[nrow(tt)]]]=matrix(0,length(pi),length(pi)) ## placeholder matrix to avoid error when passing to Rcpp
    ## Compute the marginal transitions
    margTrans=epiAllele:::marginalTransitionsCpp(data=data,tMat=logTransMat,traversal=as.matrix(tt-1),nTips=nTips,
                                                 logPi=log(pi),siblings=siblingList,ncores=ncores)
    out=list()
    for(k in 1:nrow(tt)){
      out[[paste0(tt$parent[k],"-",tt$child[k])]]=margTrans[k,,]
    }
    return(out)
  }
  return(groupTrans)
})