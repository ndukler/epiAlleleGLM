#ifndef MSGPASS_H
#define MSGPASS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "epiAllele_types.h"
using namespace Rcpp;

arma::mat postorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                                  const NumericVector& logPi, unsigned int nNode,int ncores = 1);
arma::mat preorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, 
                                 const NumericVector& logPi,const NumVecList& siblings, int nNode, int root, int ncores = 1);
#endif