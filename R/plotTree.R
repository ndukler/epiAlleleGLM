#' Tree visualization
#'
#' Plots tree with various helpful highlightings
#' @param obj alleleData
#' @param colorByRate Can be set to either color by the grouping of rates or their actual values. Only usable with rateModel object.
#' @param scaleBranches If TRUE, branches are scaled by rate (default = FALSE)
#' @param offset How much to offset tip labels by
#' @param xmax Maximum value of x axis in tree plots (Note: if you want different x axis for different tree 
#' plots use \link[ggtree]{xlim_expand} on the returned object)
#' @param nodeLabels if TRUE label nodes with numbers
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
methods::setMethod("plotTree", signature(obj = "alleleData"), function(obj,offset=0.1,xmax=NULL,nodeLabels=TRUE) {
  tempTree=getTree(obj)
  g <- ggtree::ggtree(tempTree)+
    ggtree::theme_tree2()+
    ggtree::geom_tiplab(show.legend=FALSE,color="black",offset=offset)
  if(nodeLabels){g <- g+ggplot2::geom_label(ggplot2::aes(label=node), hjust=0.5)}
    if(!is.null(xmax)) {g=g+ggtree::xlim_tree(xlim = xmax)}
  return(g)
})

#' @name plotTree
#' @rdname plotTree
methods::setMethod("plotTree", signature(obj = "rateModel"), function(obj,colorByRate=c("index","value"),scaleBranches=FALSE,
                                                                      offset=0.1,xmax=NULL,nodeLabels=TRUE) {
  ## Checks and default setting
  if(length(colorByRate)>1) colorByRate="index"
  if(!colorByRate %in% c("index","value")) stop("colorByRate must be either \'index'\ or \'value\'")
  
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
    ## Scale branches in tree if option set
    if(scaleBranches){
      trList[[l]]$edge.length=tr$edge.length*temp[.id==l][order(index)]$value
    }
  }
  class(trList) <- "multiPhylo"

  ## plot building
  g <- ggtree::ggtree(trList)+
    ggtree::theme_tree2()+
    ggplot2::facet_wrap(~.id, scale="free") + 
    ggplot2::theme(legend.position = "bottom")+
    ggplot2::theme(legend.position="right")+
    ggtree::geom_tiplab(show.legend=FALSE,color="black",offset=offset)
  if(nodeLabels){g <- g+ggplot2::geom_label(ggplot2::aes(label=node), hjust=0.5)}
  if(!is.null(xmax)) {g=g+ggtree::xlim_tree(xlim = xmax)}
  ## Merge in new data
  g$data=merge(g$data,temp,by = c(".id","node"),all.x = TRUE)   
  ## Color by indicated option
  if(colorByRate=="index"){
    if (requireNamespace("randomcoloR", quietly = TRUE)){
      g=g+ggplot2::aes(color=index)+
        ggplot2::scale_color_manual(values=randomcoloR::distinctColorPalette(length(levels(g$data$index))),breaks=levels(temp$index))
    } else {
      g=g+ggplot2::aes(color=index)
    }
  } else if(colorByRate=="value"){
    g=g+ggplot2::aes(color=value)+
      ggplot2::scale_color_viridis_c(end=0.7)
  }
  
  return(g)
})
