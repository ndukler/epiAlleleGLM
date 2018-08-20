#' Get edgeGroup
#'
#' Returns edgeGroups table as data.table 
#' @param obj rateModel object
#' @param parent node number(s) for parents
#' @param child node number(s) for children
#' @name getEdgeGroup
#' @include rateModel-class.R
#' @rdname getEdgeGroup
#' @examples
#' 
#' @export
methods::setGeneric("getEdgeGroup", function(obj) {
  standardGeneric("getEdgeGroup")
})

#' @name getEdgeGroup
#' @rdname getEdgeGroup
methods::setMethod("getEdgeGroup", signature(obj = "rateModel"), function(obj) {
  
})