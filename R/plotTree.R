#' Tree visualization
#'
#' Plots tree with various helpful highlightings
#' @param obj alleleData
#' @param colorByRate Can be set to either color by the grouping of rates or their actual values. Only usable with rateModel object.
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
  temp=data.table::rbindlist(temp,idcol=TRUE)
  temp[,index:=factor(index)]
  temp[,value:=obj@paramEnviron$params[index]]
  class(trList) <- "multiPhylo"
  

  ## plot building
  g <- ggtree::ggtree(trList)+
    ggplot2::geom_label(ggplot2::aes(label=node),color="black", hjust=0.5)+
    ggtree::theme_tree2()+
    ggplot2::facet_wrap(~.id, scale="free") + 
    ggplot2::theme(legend.position = "bottom")+
    ggplot2::theme(legend.position="right")
  ## Merge in new data
  g$data=merge(g$data,temp,by = c(".id","node"),all.x = TRUE)   
  ## Color by indicated option
  if(colorByRate=="index"){
    g=g+ggplot2::aes(color=index)+
      ggplot2::scale_color_viridis_d(breaks=levels(temp$index),end=0.7)
  } else if(colorByRate=="value"){
    g=g+ggplot2::aes(color=value)+
      ggplot2::scale_color_viridis_c(end=0.7)
  }
  
  return(g)
})
