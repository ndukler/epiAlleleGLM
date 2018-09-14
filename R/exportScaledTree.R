#' Exports scaled tree
#'
#' Exports list of scaled trees, one per entry per site label
#' @param obj rateModel
#' @name exportScaledTree
#' @include rateModel-class.R
#' @rdname exportScaledTree
#' @export
methods::setGeneric("exportScaledTree", function(obj,...) {
  standardGeneric("exportScaledTree")
})

#' @name exportScaledTree
#' @rdname exportScaledTree
methods::setMethod("exportScaledTree", signature(obj = "rateModel"), function(obj) {
  ## Compute edge coloring
  tr=getTree(obj)
  temp=list()
  trList=list()
  for(l in levels(obj@siteLabels$siteLabel)){
    temp[[l]]=data.table::data.table(node = getEdgeGroupTable(obj)$child,index=getRateIndex(obj,edges = getEdgeGroupTable(obj)[,.(parent,child)],siteLabel = l))
  }
  temp=data.table::rbindlist(temp,idcol=TRUE)
  ## Add rate values to temp data.table
  temp[,value:=getParamValue(obj,index)]
  temp[,index:=factor(index)]
  
  ## Create list of trees
  for(l in levels(obj@siteLabels$siteLabel)){
    trList[[l]]=tr
    trList[[l]]$edge.length=tr$edge.length*temp[.id==l][order(index)]$value
  }
  class(trList) <- "multiPhylo"
  return(trList)
})
