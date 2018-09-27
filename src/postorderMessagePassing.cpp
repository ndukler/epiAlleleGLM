#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"

// Forward message passing for a single site
// [[Rcpp::export]]
arma::mat postorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                                  const NumericVector& logPi,unsigned int nNode,int ncores = 1) {
  unsigned int nAlleles=logPi.size();
  arma::mat poTab(nNode,nAlleles,arma::fill::zeros);  
  
  // Initialize values for the tips of the tree
  for(unsigned int n=0;n<nTips;n++){
    for(unsigned int a=0;a<nAlleles;a++){
      poTab(n,a)=data((n*nAlleles)+a);
    }
  }
  // Now compute the probability for the interior nodes
  for(unsigned int n=0;n<traversal.nrow();n++){
    unsigned int parentInd=traversal(n,0);
    unsigned int childInd=traversal(n,1);
    arma::mat tMatA = as<arma::mat>(tMat[n]);
    for(unsigned int a=0;a<nAlleles;a++){ // iterate over all parental alleles
      poTab(parentInd,a) = poTab(parentInd,a) + logSumExpArma((arma::vec) poTab(arma::span(childInd), arma::span::all)+ tMatA.row(a).t());
    }
  }
  return(poTab);
}
