#' Get edgeGroupTable
#'
#' Returns edgeGroups table as data.table 
#' @param obj rateModel object
#' @name getEdgeGroupTable
#' @include rateModel-class.R
#' @rdname getEdgeGroupTable
#' @examples
#' 
#' @export
methods::setGeneric("getEdgeGroupTable", function(obj) {
  standardGeneric("getEdgeGroupTable")
})

#' @name getEdgeGroupTable
#' @rdname getEdgeGroupTable
methods::setMethod("getEdgeGroupTable", signature(obj = "rateModel"), function(obj) {
  return(obj@edgeGroups)
})