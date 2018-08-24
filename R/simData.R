#' Simulate data
#'
#' Simulates data without an error model
#' @param nSites number of sites to simulate
#' @param tr tree
#' @param rate a vector, either of length one, or the same as the number of branches
#' @param pi stationary distribution of allele frequencies
#' @name simData
#' @return a list
#' @export 
simData <- function(nSites,tr,rate,pi){
  ## Unroot tree if it's rooted
  if(ape::is.rooted(tr)){
    write("Unrooting tree...")
    tr=ape::unroot(tr)
  } 
  ## Check that rate is either of length one or has the same length as there are tree edges
  if(length(rate) !=1 && length(rate) != length(tr$edge.length)){
    stop("Rate must either be one or the same as the number of edges in the UNROOTED tree.")
  }
  if(length(rate)==1){
    rate=rep(rate,length(tr$edge.length))
  }
  ## Rescale tree edges
  trRescale=tr
  trRescale$edge.length=tr$edge.length*rate
  ## Normalize the stationary frequency
  pi=pi/sum(pi)
  nAlleles=length(pi)
  ## Normalize the rate
  temp=matrix(1,ncol=nAlleles,nrow = nAlleles)
  diag(temp)=0
  Q=temp %*% diag(pi) ## constuction that guarentees detailed balance: pi_i*q_ij = pi_j*q_ji
  diag(Q)=-rowSums(Q)
  normRate=rate/sum(-diag(Q)*pi)
  ## Create matrix to hold simulated data
  simDat=matrix(nrow = nSites,ncol=length(tr$tip.label))
  colnames(simDat)=tr$tip.label
  for(i in 1:nSites){
    simDat[i,] <- ape::rTraitDisc(phy = trRescale,rate=1,k = length(pi),freq=pi,ancestor = FALSE,
                                  root.value = sample(x=1:length(pi),size = 1,prob = pi))
  }
  aData=lapply(split(t(simDat), f =colnames(simDat)),function(x){
    z=matrix(0,nrow = length(x),ncol=nAlleles)
    for(k in 1:nrow(z)){z[k,x[k]]=1}
    return(z)
  })
  ## Create edgeGroup table
  eTab=data.table::data.table(parent=trRescale$edge[,1],child=trRescale$edge[,2],edgeID=paste(trRescale$edge[,1],trRescale$edge[,2],sep="-"),
                         edgeGroup=paste0("e",as.numeric(factor(rate))))
  ## Map of edge groups to rates
  rateMap=unique(data.table::data.table(edgeGroup=paste0("e",as.numeric(factor(rate))),rate))
  return(list(data=aData,tr=tr,edgeTable=eTab,rateMap=rateMap,pi=pi))
}