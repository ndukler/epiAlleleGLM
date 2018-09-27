// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "epiAllele_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logSumExp
double logSumExp(NumericVector x);
RcppExport SEXP _epiAllele_logSumExp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logSumExp(x));
    return rcpp_result_gen;
END_RCPP
}
// marginalTransitionsCpp
arma::cube marginalTransitionsCpp(const NumericMatrix& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, const NumericVector& logPi, const NumVecList& siblings, int ncores);
RcppExport SEXP _epiAllele_marginalTransitionsCpp(SEXP dataSEXP, SEXP tMatSEXP, SEXP traversalSEXP, SEXP nTipsSEXP, SEXP logPiSEXP, SEXP siblingsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumMatList& >::type tMat(tMatSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traversal(traversalSEXP);
    Rcpp::traits::input_parameter< const double >::type nTips(nTipsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const NumVecList& >::type siblings(siblingsSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(marginalTransitionsCpp(data, tMat, traversal, nTips, logPi, siblings, ncores));
    return rcpp_result_gen;
END_RCPP
}
// multiProbToStick
NumericVector multiProbToStick(const NumericVector& x, int width);
RcppExport SEXP _epiAllele_multiProbToStick(SEXP xSEXP, SEXP widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    rcpp_result_gen = Rcpp::wrap(multiProbToStick(x, width));
    return rcpp_result_gen;
END_RCPP
}
// multiStickToProb
NumericVector multiStickToProb(const NumericVector& x, int width);
RcppExport SEXP _epiAllele_multiStickToProb(SEXP xSEXP, SEXP widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    rcpp_result_gen = Rcpp::wrap(multiStickToProb(x, width));
    return rcpp_result_gen;
END_RCPP
}
// postorderMessagePassing
NumericMatrix postorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const int nTips, const NumericVector& logPi, unsigned int nNode);
RcppExport SEXP _epiAllele_postorderMessagePassing(SEXP dataSEXP, SEXP tMatSEXP, SEXP traversalSEXP, SEXP nTipsSEXP, SEXP logPiSEXP, SEXP nNodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumMatList& >::type tMat(tMatSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traversal(traversalSEXP);
    Rcpp::traits::input_parameter< const int >::type nTips(nTipsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nNode(nNodeSEXP);
    rcpp_result_gen = Rcpp::wrap(postorderMessagePassing(data, tMat, traversal, nTips, logPi, nNode));
    return rcpp_result_gen;
END_RCPP
}
// preorderMessagePassing
NumericMatrix preorderMessagePassing(const NumericVector& data, const NumMatList& tMat, const NumericMatrix& traversal, const int nTips, const NumericVector& logPi, const NumericMatrix& alpha, const NumVecList& siblings, int nNode, int root);
RcppExport SEXP _epiAllele_preorderMessagePassing(SEXP dataSEXP, SEXP tMatSEXP, SEXP traversalSEXP, SEXP nTipsSEXP, SEXP logPiSEXP, SEXP alphaSEXP, SEXP siblingsSEXP, SEXP nNodeSEXP, SEXP rootSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumMatList& >::type tMat(tMatSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traversal(traversalSEXP);
    Rcpp::traits::input_parameter< const int >::type nTips(nTipsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumVecList& >::type siblings(siblingsSEXP);
    Rcpp::traits::input_parameter< int >::type nNode(nNodeSEXP);
    Rcpp::traits::input_parameter< int >::type root(rootSEXP);
    rcpp_result_gen = Rcpp::wrap(preorderMessagePassing(data, tMat, traversal, nTips, logPi, alpha, siblings, nNode, root));
    return rcpp_result_gen;
END_RCPP
}
// probToStick
NumericVector probToStick(NumericVector x);
RcppExport SEXP _epiAllele_probToStick(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(probToStick(x));
    return rcpp_result_gen;
END_RCPP
}
// setValues
void setValues(NumericVector& x, NumericVector& ind, NumericVector& val);
RcppExport SEXP _epiAllele_setValues(SEXP xSEXP, SEXP indSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type ind(indSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type val(valSEXP);
    setValues(x, ind, val);
    return R_NilValue;
END_RCPP
}
// stickToProb
NumericVector stickToProb(NumericVector x);
RcppExport SEXP _epiAllele_stickToProb(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(stickToProb(x));
    return rcpp_result_gen;
END_RCPP
}
// treeLL
NumericVector treeLL(const NumericMatrix& data, const NumMatList& tMat, const NumericMatrix& traversal, const double nTips, const NumericVector& logPi);
RcppExport SEXP _epiAllele_treeLL(SEXP dataSEXP, SEXP tMatSEXP, SEXP traversalSEXP, SEXP nTipsSEXP, SEXP logPiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumMatList& >::type tMat(tMatSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traversal(traversalSEXP);
    Rcpp::traits::input_parameter< const double >::type nTips(nTipsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logPi(logPiSEXP);
    rcpp_result_gen = Rcpp::wrap(treeLL(data, tMat, traversal, nTips, logPi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epiAllele_logSumExp", (DL_FUNC) &_epiAllele_logSumExp, 1},
    {"_epiAllele_marginalTransitionsCpp", (DL_FUNC) &_epiAllele_marginalTransitionsCpp, 7},
    {"_epiAllele_multiProbToStick", (DL_FUNC) &_epiAllele_multiProbToStick, 2},
    {"_epiAllele_multiStickToProb", (DL_FUNC) &_epiAllele_multiStickToProb, 2},
    {"_epiAllele_postorderMessagePassing", (DL_FUNC) &_epiAllele_postorderMessagePassing, 6},
    {"_epiAllele_preorderMessagePassing", (DL_FUNC) &_epiAllele_preorderMessagePassing, 9},
    {"_epiAllele_probToStick", (DL_FUNC) &_epiAllele_probToStick, 1},
    {"_epiAllele_setValues", (DL_FUNC) &_epiAllele_setValues, 3},
    {"_epiAllele_stickToProb", (DL_FUNC) &_epiAllele_stickToProb, 1},
    {"_epiAllele_treeLL", (DL_FUNC) &_epiAllele_treeLL, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_epiAllele(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
