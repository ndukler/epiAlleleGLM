#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' probToStick 
//'
//' Convert a set of probabilities that sum to one to a set of stick breaking parameters
//' @param x parameter vector of probabilities
//' @name probToStick
//' @return a vector of parameters for a stick breaking parameters
// [[Rcpp::export]]
NumericVector probToStick(NumericVector x) {
  NumericVector stick(x.size()-1);
  double remain=1;
  for(int i=0; i<x.size()-1;i++){
    stick(i)=x(i)/remain;
    remain=remain-x(i);
  }
  return(stick);
}