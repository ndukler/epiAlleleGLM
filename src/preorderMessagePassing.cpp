#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"

// [[Rcpp::export]]
arma::mat preorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                                   const NumericVector& logPi,const NumVecList& siblings, int nNode, int root, int ncores = 1) {
  unsigned int nAlleles=logPi.size();
  arma::mat poTab(nNode,nAlleles,arma::fill::zeros);  
  
  // Initialize the root
  for(unsigned int a=0;a<nAlleles;a++){
      poTab(root,a)=logPi(a);
  }
  // Now compute the probability for the interior nodes
  for(int n=traversal.nrow()-1;n>=0;n--){
    unsigned int parentInd=traversal(n,0);
    unsigned int childInd=traversal(n,1);
    arma::mat tMatA = as<arma::mat>(tMat[n]); 
    for(unsigned int a=0;a<nAlleles;a++){ // iterate over all parental alleles
      poTab(childInd,a) = logSumExpArma((arma::vec) poTab(arma::span(parentInd), arma::span::all) + tMatA.col(a));
    }
  }
  return(poTab);
}
