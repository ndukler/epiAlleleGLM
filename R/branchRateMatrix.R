#' Compute exponentiated rate matrix for all branches
#'
#' @param rate vector of mutation rates, one per branch
#' @param branch.length vector of branch lengths
#' @param pi stationary distribution of allele frequencies
#' @param log if TRUE returns log transition matricies
#' @name branchRateMatrix
#' @return a list of transition matricies, one entry for each specified branch length. The rows are the parental allele,
#' the columns are the child allele.
branchRateMatrix <- function(rate,branch.length,pi,log=TRUE){
  ## Check that there are the same number of rates as there are branch lengths
  if(length(rate) != length(branch.length)){
    stop("There must be the same number of rates as branch lengths")
  }
  nAlleles=length(pi)
  ## Construct rate matrix
  temp=matrix(1,ncol=nAlleles,nrow = nAlleles)
  diag(temp)=0
  Q=temp %*% diag(pi) ## constuction that guarentees detailed balance: pi_i*q_ij = pi_j*q_ji
  diag(Q)=-rowSums(Q)
  ## Standardize rate matrix
  norm=1/sum(-diag(Q)*pi)
  qN=Q*norm
  ## Compute per branch scaling factor
  sf = rate * branch.length
  ## Exponentiate rate matrix to transition matrix
  if(log==TRUE){
    qList=lapply(as.list(sf),function(x) log(as.matrix(Matrix::expm(qN*x))))
  } else {
    qList=lapply(as.list(sf),function(x) as.matrix(Matrix::expm(qN*x)))
  }
  return(qList)
}