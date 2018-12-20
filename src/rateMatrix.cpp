#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rateMatrix(const arma::vec & pi,double rate, double branchLength) {
  int nAlleles = pi.size();
  // initalize temp matrix with all ones
  arma::mat temp(nAlleles,nAlleles,arma::fill::ones);
  // Set diagonal elements to zero
  temp.diag().zeros();
  // constuction that guarentees detailed balance: pi_i*q_ij = pi_j*q_ji
  arma::mat piM = arma::diagmat(pi);
  arma::mat Q=temp*piM; 
  Q.diag()= -arma::sum(Q,1); 
  // Standardize rate matrix and scale by branch length and rate
  double norm = 1/sum(-Q.diag()%pi);
  Q*=norm*branchLength*rate;
  // exponentiate rate matrix and return
  return(arma::expmat(Q));
}