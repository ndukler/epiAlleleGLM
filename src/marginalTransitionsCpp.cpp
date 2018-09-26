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
arma::cube marginalTransitionsCpp(const NumericMatrix& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                                  const NumericVector& logPi,const NumVecList& siblings, int ncores = 1) {
  unsigned int nAlleles=logPi.size();
  // Create cube to hold the expected number of transitions over all sites
  arma::cube expectedTransitions(traversal.nrow(),nAlleles,nAlleles,arma::fill::zeros);  
  
  
}
