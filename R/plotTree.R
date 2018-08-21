#' Tree visualization
#'
#' Plots tree with various helpful highlightings
#' @param x alleleData
#' @name plotTree
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @rdname plotTree
#' @examples
#' 
#' @export
methods::setGeneric("plotTree", function(obj,...) {
  standardGeneric("plotTree")
})

#' @name plotTree
#' @rdname plotTree
methods::setMethod("plotTree", signature(obj = "alleleData"), function(obj) {
  g <- ggtree::ggtree(getTree(obj))+
    ggplot2::geom_label(ggplot2::aes(label=node), hjust=0.5)+
    ggtree::theme_tree2()
  return(g)
})

#' @name plotTree
#' @rdname plotTree
methods::setMethod("plotTree", signature(obj = "rateModel"), function(obj,colorByRate=c("index","value")) {
  ## Checks and default setting
  if(length(colorByRate)>1) colorByRate="index"
  if(!colorByRate %in% c("index","value")) stop("colorByRate must be either \'index'\ or \'value\'")
  
  ## Compute edge coloring
  tr=getTree(obj)
  temp=list()
  trList=list()
  for(l in levels(obj@siteLabels$siteLabel)){
    temp[[l]]=data.table::data.table(node = getEdgeGroupTable(obj)$child,index=getRateIndex(obj,edges = getEdgeGroupTable(obj)[,.(parent,child)],siteLabel = l))
    trList[[l]]=tr
  }
  temp=data.table::rbindlist(temp,idcol="siteLabel")
  temp[,color:=viridisLite::viridis(max(index))[index]]
  class(trList) <- "multiPhylo"
  

  
  ## plot building
  g <- ggtree::ggtree(trList)+
    ggplot2::geom_label(ggplot2::aes(label=node), hjust=0.5)+
    ggtree::theme_tree2()+
    ggplot2::facet_wrap(~siteLabel, scale="free")
  
  g %<+% temp[index==1] + aes(color=color) 
  return(g)

  
})
