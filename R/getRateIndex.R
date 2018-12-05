#' Get rate parameter index
#'
#' Returns index of a rate parameter associated with a given edge/edgeGroup and siteLabel 
#' @param obj rateModel object
#' @param edges Can either be a data.frame/data.table with two columns, \'parent\' and \'child\' or a character vector of edgeIDs
#' @param edgeGroup Instead of the edges argument, a vector of edge groups can be supplied. Only one can be specified.
#' @name getRateIndex
#' @include rateModel-class.R
#' @rdname getRateIndex
#' @return A numeric vector of parameter indices
#' @examples
#' 
#' @export
methods::setGeneric("getRateIndex", function(obj,edges,edgeGroup) {
  standardGeneric("getRateIndex")
})

#' @name getRateIndex
#' @rdname getRateIndex
methods::setMethod("getRateIndex", signature(obj = "rateModel",edges="data.frame",edgeGroup="missing"), function(obj,edges=NULL,edgeGroup) {
    if(!setequal(colnames(edges),c("parent","child"))){
      stop("Table must contain columns \'parent\' and \'child\'")  
    } 
    edgeIDs=paste(edges$parent,edges$child,sep = "-")
    getRateIndex(obj,edgeIDs)
})

#' @name getRateIndex
#' @rdname getRateIndex
methods::setMethod("getRateIndex", signature(obj = "rateModel",edges="character",edgeGroup="missing"), function(obj,edges,edgeGroup) {
  eg=getEdgeGroupTable(obj)[.(edges)]$edgeGroup
  ind=obj@rateIndex[.(edgeGroup=eg)]$rateIndex
  return(ind)
})

#' @name getRateIndex
#' @rdname getRateIndex
methods::setMethod("getRateIndex", signature(obj = "rateModel",edges="missing",edgeGroup="character"), function(obj,edges,edgeGroup) {
  eg=edgeGroup
  ind=obj@rateIndex[.(eg)]$index
  return(ind)
})
