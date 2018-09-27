#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"
#include "msgPassing.h"
using namespace Rcpp;

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
    // Do the message passing
    NumericMatrix alpha=postorderMessagePassing((Rcpp::NumericVector) data(i,_),tMat,traversal,nTips,logPi,nNode);
    NumericMatrix beta=preorderMessagePassing((Rcpp::NumericVector) data(i,_),tMat,traversal,nTips,logPi,alpha,siblings,nNode,root);
    for(int e=0; e<traversal.nrow();e++){
      int parentInd=traversal(e,0);
      int childInd=traversal(e,1);
      for(int a=0;a<nAlleles;a++){ // iterate over parent alleles
        for(int b=0;b<nAlleles;b++){ // iterate over child alleles
          expectedTransitions(e,a,b) = beta(parentInd,a)+ tMat[childInd](a,b) + alpha(childInd,b);
        }
      }
      // Compute partition function
      double Z = logSumExpArma(arma::vectorise((arma::mat) expectedTransitions(arma::span(e),arma::span::all,arma::span::all)));
      expectedTransitions(arma::span(e),arma::span::all,arma::span::all)=expectedTransitions(arma::span(e),arma::span::all,arma::span::all)-Z;
    }
  }
  return(expectedTransitions);
}
