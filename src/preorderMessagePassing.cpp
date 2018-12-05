#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAlleleGLM_types.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix preorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const int nTips, 
                                   const NumericVector& logPi,const NumericMatrix& alpha,const NumVecList& siblings, int nNode, int root) {
  int nAlleles=logPi.size();
  NumericMatrix poTab(nNode,nAlleles);  
  // Initialize the root
  for(int a=0;a<nAlleles;a++){
      poTab(root,a)=logPi(a);
  }
  // Now compute the probability for the interior nodes
  for(int n=traversal.nrow()-1;n>=0;n=n-1){
    int parentInd=traversal(n,0);
    int childInd=traversal(n,1);
    for(unsigned int a=0;a<nAlleles;a++){ // iterate over all alleles of focal node
      NumericVector msgHolder(nAlleles);
      for(unsigned int b=0;b<nAlleles;b++){ // iterate over all states of parental node summing over all sibblings
        double parentContrib = poTab(parentInd,b) + tMat[childInd](b,a); // calculate the parent contribution
        double sibContrib = 0; // holds the sibling contribution
        for(int s=0;s<siblings[childInd].size();s++){ // iterate over siblings
          sibContrib=sibContrib+logSumExp(alpha(siblings[childInd](s),_)+tMat[siblings[childInd](s)](b,_));
        }
        msgHolder(b)=parentContrib+sibContrib;
      }
      poTab(childInd,a)=logSumExp(msgHolder);
    }
  }
  return(poTab);
}