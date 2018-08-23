#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void setValues(NumericVector& x, NumericVector& ind, NumericVector& val) {
  for(int i=0; i<ind.size();i++){
    x(ind(i)-1)=val(i);
  }
}