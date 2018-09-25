#ifndef STICKFUN_H
#define STICKFUN_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

NumericVector probToStick(NumericVector x);
NumericVector stickToProb(NumericVector x);

#endif