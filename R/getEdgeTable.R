methods::setGeneric("getEdgeTable", function(obj) {
  standardGeneric("getEdgeTable")
})

#' Get edge table
#'
#' Returns edge table. Useful for construction of rateModel class.
#' @param x alleleData
#' @name getEdgeTable
#' @return data.table containing edge ids
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("getEdgeTable", signature(obj = "alleleData"), function(obj) {
  data.table::data.table(parent=obj@tree$edge[,1],child=obj@tree$edge[,2])
})
  
