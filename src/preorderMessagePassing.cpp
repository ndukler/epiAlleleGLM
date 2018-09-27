#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"

// [[Rcpp::export]]
arma::mat preorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                                   const NumericVector& logPi,const NumericMatrix& alpha,const NumVecList& siblings, int nNode, int root, int ncores = 1) {
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
    NumericVector msgHolder(nAlleles*nAlleles);
    int iter=0;
    for(unsigned int a=0;a<nAlleles;a++){ // iterate over all parental alleles
      for(unsigned int b=0;b<nAlleles;b++){ // iterate over all sibbling alleles
        // Parental contibution to message
        
        // Sibbing contribution to message
        iter++;
      }
    }
  }
  return(poTab);
}
