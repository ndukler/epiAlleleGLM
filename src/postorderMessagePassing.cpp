#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"

// Forward message passing for a single site
// [[Rcpp::export]]
NumericMatrix postorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const int nTips, 
                                  const NumericVector& logPi,unsigned int nNode) {
  unsigned int nAlleles=logPi.size();
  NumericMatrix poTab(nNode,nAlleles);  
  
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
    for(unsigned int a=0;a<nAlleles;a++){ // iterate over all parental alleles
      poTab(parentInd,a) = poTab(parentInd,a) + logSumExp(poTab(childInd,_) + tMat[childInd](a,_));
    }
  }
  return(poTab);
}
