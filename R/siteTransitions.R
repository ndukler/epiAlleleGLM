#' Compute transitions per branch
#'
#' Compute the number of expected transitions per site between states based on the marginal joint probabilities
#' @param obj rateModel object
#' @param subset The site indices for which should be summed over to compute the total number of expected transitions (NOT IMPLEMENTED)
#' @param ncores Number of cores to use
#' @return A data.table where each is labeled by site index, parent allele, and child allele, and the number transitions over the whole tree
#' @name siteTransitions
#' @include rateModel-class.R
#' @rdname siteTransitions
#' @examples
#' 
#' @export
methods::setGeneric("siteTransitions", function(obj,subset=NULL,ncores=1) {
  standardGeneric("siteTransitions")
})

#' @name siteTransitions
#' @rdname siteTransitions
methods::setMethod("siteTransitions", signature(obj = "rateModel"), function(obj,subset=NULL,ncores=1) {
  siteTypes=levels(obj@siteLabels$siteLabel)
  nTips=length(getTree(obj)$tip.label)
  ## Create a list of all siblings for each node 
  siblingList=list()
  for(n in sort(unique(as.numeric(getTree(obj)$edge)))){
    siblingList[[n]]=getSiblings(getTree(obj),n)-1
  }
  #  l = siteTypes[1]
  ## `%myPar%` <- ifelse(foreach::getDoParRegistered(), yes = foreach::`%dopar%`, no = foreach::`%do%`)
  groupTrans=foreach::foreach(l=siteTypes,.inorder = TRUE, .final = function(x) data.table::rbindlist(x)) %do% {
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
    margTrans=siteGainLossCpp(data=data,tMat=logTransMat,traversal=as.matrix(tt-1),nTips=nTips,
                                                 logPi=log(pi),siblings=siblingList,ncores=ncores)
    ## Place cube in table
    ind=getLabelIndices(obj,l)
    out=data.table::as.data.table(expand.grid(site=getLabelIndices(obj,l),parent=1:length(pi),child=1:length(pi)),transitions=0)
    data.table::setkeyv(out,cols = c("site","parent","child"))
    for(i in 1:nrow(data)){
      mm=reshape2::melt(margTrans[i,,])
      out[.(ind[i],mm$Var1,mm$Var2),transitions:=mm$value]
    }
    return(out)
  }
  return(groupTrans[order(site)])
})