#ifndef MSGPASS_H
#define MSGPASS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "epiAlleleGLM_types.h"
using namespace Rcpp;

NumericMatrix postorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const int nTips, 
                                  const NumericVector& logPi,unsigned int nNode);
NumericMatrix preorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const int nTips, 
                                     const NumericVector& logPi,const NumericMatrix& alpha,const NumVecList& siblings, int nNode, int root);
#endif