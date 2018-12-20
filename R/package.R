#' epiAllele
#'
#' Description of your package
#'
#' @docType package
#' @author Noah Dukler <ndukler@cshl.edu>
#' @importFrom data.table :=
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom ggtree %<+%
#' @import Rcpp
#' @useDynLib epiAlleleGLM
#' @name epiAlleleGLM
NULL

Rcpp::loadModule(module = "paramIndex", what = TRUE)