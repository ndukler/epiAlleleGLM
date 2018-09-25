#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logSumExp.h"
#include "epiAllele_types.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector treeLL(const NumericMatrix& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                     const NumericVector& logPi) {
  int sites=data.nrow();
  int nAlleles=logPi.size();
  NumericVector siteLik(sites); // Numeric vector of the for the logProbability of each site
  int nNode=Rcpp::max(traversal(_,0))+1; // The total number of nodes on the tree
  int root=traversal(traversal.nrow()-1,0); // Get the index of the root node

  // loop over sites
  for(int i=0;i<sites;i++){
    NumericMatrix nodeLogProb(nNode,nAlleles);
    // Initialize values for the tips of the tree
    for(int n=0;n<nTips;n++){
      for(int a=0;a<nAlleles;a++){
        nodeLogProb(n,a)=data(i,(n*nAlleles)+a);
      }
    }
    // Rcout << nodeLogProb << std::endl;
    // Now compute the probability for the interior nodes
    for(int n=0;n<traversal.nrow();n++){
      int parentInd=traversal(n,0);
      int childInd=traversal(n,1);
      for(int a=0;a<nAlleles;a++){ // iterate over all parental alleles
        nodeLogProb(parentInd,a) = nodeLogProb(parentInd,a) + logSumExp(nodeLogProb(childInd,_)+tMat[n](a,_));
      }
    }
    siteLik(i)=logSumExp(nodeLogProb(root,_)+logPi);
  }
  return(siteLik);
}
