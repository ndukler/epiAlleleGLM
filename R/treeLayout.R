#' Returns a data table defining the line segments of a phylogenetic tree.
#'
#' This function takes a \code{\link{phylo}} or \code{\link{phyloseq-class}} object
#' and returns a list of two \code{\link{data.table}}s suitable for plotting
#' a phylogenetic tree with \code{\link[ggplot2]{ggplot}}2.
#' 
#' @param phy (Required). The \code{\link{phylo}} or \code{\link{phyloseq-class}}
#'  object (which must contain a \code{\link{phylo}}genetic tree)
#'  that you want to converted to \code{\link{data.table}}s
#'  suitable for plotting with \code{\link[ggplot2]{ggplot}}2.
#'
#' @param ladderize (Optional). Boolean or character string (either
#'  \code{FALSE}, \code{TRUE}, or \code{"left"}).
#'  Default is \code{FALSE} (no ladderization).
#'  This parameter specifies whether or not to \code{\link[ape]{ladderize}} the tree 
#'  (i.e., reorder nodes according to the depth of their enclosed
#'  subtrees) prior to plotting.
#'  This tends to make trees more aesthetically pleasing and legible in
#'  a graphical display.
#'  When \code{TRUE} or \code{"right"}, ``right'' ladderization is used.
#'  When set to \code{FALSE}, no ladderization is applied.
#'  When set to \code{"left"}, the reverse direction
#'  (``left'' ladderization) is applied.
#'  
#' @return
#'  A list of two \code{\link{data.table}}s, containing respectively 
#'  a \code{data.table} of edge segment coordinates, named \code{edgeDT},
#'  and a \code{data.table} of vertical connecting segments, named \code{vertDT}.
#'  See \code{example} below for a simple demonstration.
#' 
#' @seealso
#' An early example of this functionality was borrowed directly, with permission,
#' from the package called \code{ggphylo}, 
#' released on GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' by its author Gregory Jordan \email{gjuggler@@gmail.com}.
#' That original phyloseq internal function, \code{tree.layout}, has been
#' completely replaced by this smaller and much faster user-accessible 
#' function that utilizes performance enhancements from standard 
#' \code{\link{data.table}} magic as well as \code{\link{ape-package}}
#' internal C code.
#' 
#' @export
#' @examples
#' library("ggplot2")
#' data("esophagus")
#' phy = phy_tree(esophagus)
#' phy <- ape::root(phy, "65_2_5", resolve.root=TRUE)
#' treeSegs0 = tree_layout(phy)
#' treeSegs1 = tree_layout(esophagus)
#' edgeMap = aes(x=xleft, xend=xright, y=y, yend=y)
#' vertMap = aes(x=x, xend=x, y=vmin, yend=vmax)
#' p0 = ggplot(treeSegs0$edgeDT, edgeMap) + geom_segment() + geom_segment(vertMap, data=treeSegs0$vertDT)
#' p1 = ggplot(treeSegs1$edgeDT, edgeMap) + geom_segment() + geom_segment(vertMap, data=treeSegs1$vertDT)
#' print(p0)
#' print(p1)
#' plot_tree(esophagus, "treeonly")
#' plot_tree(esophagus, "treeonly", ladderize="left")
treeLayout = function(phy, ladderize=FALSE){
  if(inherits(phy, "phyloseq")){
    phy = phy_tree(phy)
  }
  if(!inherits(phy, "phylo")){
    stop("tree missing or invalid. Please check `phy` argument and try again.")
  }
  if(is.null(phy$edge.length)){
    # If no edge lengths, set them all to value of 1 (dendrogram).
    phy$edge.length <- rep(1L, times=nrow(phy$edge))
  }
  # Perform ladderizing, if requested
  if(ladderize != FALSE){
    if(ladderize == "left"){
      phy <- ladderize(phy, FALSE)
    } else if(ladderize==TRUE | ladderize=="right"){
      phy <- ape::ladderize(phy, TRUE)
    } else {
      stop("You did not specify a supported option for argument `ladderize`.")
    }
  }
  # 'z' is the tree in postorder order used in calls to .C
  # Descending order of left-hand side of edge (the ancestor to the node)
  z = ape::reorder.phylo(phy, order="postorder")
  # Initialize some characteristics of the tree.
  Ntip = length(phy$tip.label)
  ROOT = Ntip + 1
  nodelabels = phy$node.label
  # Horizontal positions
  xx = ape::node.depth.edgelength(phy)
  # vertical positions
  yy = ape::node.height(phy = phy, clado.style = FALSE)
  # Initialize an edge data.table 
  # Don't set key, order matters
  edgeDT = data.table::data.table(phy$edge, edge.length=phy$edge.length, OTU=NA_character_)
  # Add tip.labels if present
  if(!is.null(phy$tip.label)){
    # Initialize OTU, set node (V2) as key, assign taxa_names as OTU label
    edgeDT[, OTU:=NA_character_]
    data.table::setkey(edgeDT, V2)
    edgeDT[V2 <= Ntip, OTU:=phy$tip.label]
  }
  # Add the mapping for each edge defined in `xx` and `yy` 
  edgeDT[, xleft:=xx[V1]]
  edgeDT[, xright:=xx[V2]]
  edgeDT[, y:=yy[V2]]
  # Next define vertical segments
  vertDT = edgeDT[, list(x=xleft[1], vmin=min(y), vmax=max(y)), by=V1, mult="last"]
  if(!is.null(phy$node.label)){
    # Add non-root node labels to edgeDT
    edgeDT[V2 > ROOT, x:=xright]
    edgeDT[V2 > ROOT, label:=phy$node.label[-1]]
    # Add root label (first node label) to vertDT
    setkey(vertDT, V1)
    vertDT[J(ROOT), y:=mean(c(vmin, vmax))]
    vertDT[J(ROOT), label:=phy$node.label[1]]
  }
  return(list(edgeDT=edgeDT, vertDT=vertDT))
}