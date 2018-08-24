#include <Rcpp.h>
using namespace Rcpp;

//' stickToProb 
//'
//' Convert a set of stick breaking parameters to a set of probabilities that sum to one
//' @param x parameter vector of stick breaking process parameters
//' @name stickToProb
//' @return a vector of probabilities that sum to one
// [[Rcpp::export]]
NumericVector stickToProb(NumericVector x) {
  NumericVector prob(x.size()+1);
  double remain=1;
  for(int i=0; i<x.size();i++){
    prob(i)=remain*x(i);
    remain=remain-prob(i);
  }
  prob(x.size())=remain;
  return(prob);
}