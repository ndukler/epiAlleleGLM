#ifndef STICKFUN_H
#define STICKFUN_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector probToStick(NumericVector x);
NumericVector stickToProb(NumericVector x);

#endif