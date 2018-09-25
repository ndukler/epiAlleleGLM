#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Parallelized C++ equivalent
// [[Rcpp::export]]
arma::cube postorderMessagePassing(const NumericMatrix& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                     const NumericVector& logPi,int ncores = 1) {
  unsigned int sites=data.nrow();
  unsigned int nAlleles=logPi.size();
  unsigned int nNode=Rcpp::max(traversal(_,0))+1; // The total number of nodes on the tree
  arma::cube poTab(sites,nNode,nAlleles,arma::fill::zeros);  
    
  // Added an omp pragma directive to parallelize the loop with ncores
  #pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0;i<sites;i++){
    // Initialize values for the tips of the tree
    for(unsigned int n=0;n<nTips;n++){
      for(unsigned int a=0;a<nAlleles;a++){
        poTab(i,n,a)=data(i,(n*nAlleles)+a);
      }
    }
    // Now compute the probability for the interior nodes
    for(unsigned int n=0;n<traversal.nrow();n++){
      unsigned int parentInd=traversal(n,0);
      unsigned int childInd=traversal(n,1);
      arma::mat tMatA = as<arma::mat>(tMat[n]);
      for(unsigned int a=0;a<nAlleles;a++){ // iterate over all parental alleles
        poTab(i,parentInd,a) = poTab(i,parentInd,a) + logSumExpArma((arma::vec) poTab(arma::span(i), arma::span(childInd), arma::span::all)+ tMatA.row(a).t());
      }
    }
  }
  return(poTab);
}
