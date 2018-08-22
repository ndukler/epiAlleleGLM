// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

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

static const R_CallMethodDef CallEntries[] = {
    {"_epiAllele_logSumExp", (DL_FUNC) &_epiAllele_logSumExp, 1},
    {"_epiAllele_setValues", (DL_FUNC) &_epiAllele_setValues, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_epiAllele(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
