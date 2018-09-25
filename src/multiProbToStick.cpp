#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "stickFunctions.h"
using namespace Rcpp;

//' multiProbToStick 
//'
//' Convert a multiple sets of stick breaking parameters to a set of probabilities that sum to one
//' @param x the parameter of stick breaking parameters
//' @param width the width of the stick breakink process
//' @name multiProbToStick
//' @return a vector of probabilities that sum to one
// [[Rcpp::export]]
NumericVector multiProbToStick(const NumericVector& x, int width) {
  if(x.size() % width != 0){
    throw std::range_error("Number of probabilities must be a integer multiple of width");
  }
  int numProc=x.size()/width;
  NumericVector y((width-1)*numProc);
  int xPos=0;
  int yPos=0;
  for(int i=0;i<numProc;i++){
    y[Range(yPos,yPos+width-2)]=probToStick(x[Range(xPos,xPos+width-1)]);
    xPos=xPos+width;
    yPos=yPos+width-1;
  }
  return(y);
}