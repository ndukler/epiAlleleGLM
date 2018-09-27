#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"
#include "msgPassing.h"

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Parallelized C++ equivalent
// [[Rcpp::export]]
arma::cube marginalTransitionsCpp(const NumericMatrix& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                                  const NumericVector& logPi,const NumVecList& siblings, int ncores = 1) {
  unsigned int nAlleles=logPi.size();
  // Create cube to hold the expected number of transitions over all sites
  arma::cube expectedTransitions(traversal.nrow(),nAlleles,nAlleles,arma::fill::zeros);
  unsigned int nNode=Rcpp::max(traversal(_,0))+1; // The total number of nodes on the tree
  unsigned int root=traversal(traversal.nrow()-1,0); // Get the index of the root node
  // Iterate over all sites
  #pragma omp parallel for num_threads(ncores)
  for(int i=0;i<data.nrow();i++){
    // First do the forward message passing up to the root
    arma::mat alpha=postorderMessagePassing((Rcpp::NumericVector) data(i,_),tMat,traversal,nTips,logPi,nNode);
    arma::mat beta=preorderMessagePassing((Rcpp::NumericVector) data(i,_),tMat,traversal,nTips,logPi,siblings,nNode,root);
  }
  return(expectedTransitions);
}
