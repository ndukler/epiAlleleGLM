#' Get parameter index
#'
#' Returns parameter index as data.table 
#' @param obj rateModel object
#' @param edges Can either be a data.frame/data.table with two columns, \'parent\' and \'child\' or a character vector of edgeIDs
#' @param siteLabel A character vector of length one for one of the site labels
#' @name getRateIndex
#' @include rateModel-class.R
#' @rdname getRateIndex
#' @return A numeric vector of parameter indices
#' @examples
#' 
#' @export
methods::setGeneric("getRateIndex", function(obj,edges,siteLabel) {
  standardGeneric("getRateIndex")
})

#' @name getRateIndex
#' @rdname getRateIndex
methods::setMethod("getRateIndex", signature(obj = "rateModel",edges="data.frame"), function(obj,edges,siteLabel) {
    if(!setequal(colnames(edges),c("parent","child"))){
      stop("Table must contain columns \'parent\' and \'child\'")  
    } 
    edgeIDs=paste(edges$parent,edges$child,sep = "-")
    getRateIndex(obj,edgeIDs,siteLabel)
})

#' @name getRateIndex
#' @rdname getRateIndex
methods::setMethod("getRateIndex", signature(obj = "rateModel",edges="character"), function(obj,edges,siteLabel) {
  if(length(siteLabel)>1){
    stop("Only one site label can be specified at a time")
  }
  sl=siteLabel
  ind=obj@paramEnviron$paramIndex[.(getEdgeGroupTable(obj)[edges]$edgeGroup,sl)]$rateIndex
   return(ind)
})
