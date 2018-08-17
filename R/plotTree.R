methods::setGeneric("plotTree", function(obj,...) {
  standardGeneric("plotTree")
})

#' Tree visualization
#'
#' Plots tree with various helpful highlightings
#' @param x alleleData
#' @name plotTree
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("plotTree", signature(obj = "alleleData"), function(obj,edgeIds=FALSE) {
  ape::plot.phylo(obj@tree, edge.width = 2, label.offset = 0.1)
  ape::nodelabels()
  ape::tiplabels()
})
