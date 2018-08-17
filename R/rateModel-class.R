#' Class rateModel
#'
#' Class \code{rateModel} holds the allelic probabilities and corresponding phylogenetic tree
#'
#' @name rateModel-class
#' @rdname rateModel-class
#' @exportClass rateModel
methods::setClass("rateModel", representation(alleleData = "alleleData",rateGroups="data.table",siteLabelCriteria="character",
                                              siteLabels="rle"))

#' rateModel
#'
#' Contructs an object that holds the allele data and the phylogenetic tree
#' @param data An alleleData object
#' @param siteLabelCriteria A charecter vector of the columns in siteData contained in alleleData label genomic regions
#' @param lineageTable A table with columns parent,child,edgeID, and rateGroup, where rateGroup is used to specify how the rates are tied between the branches
#' @name rateModel
#' @return an alleleData object
#' @examples
#'
#' @export
rateModel <- function(data,siteLabelCriteria=NULL,lineageTable=NULL){
  ## Check that data is of class alleleData
  if(class(data)!="alleleData"){
    stop("data must be of class alleleData")
  }
  ## Check that all siteLabelCriteria elements are columns in alleleData@siteInfo
  if(is.null(siteLabelCriteria)){
    siteLabels=rle(rep("1",nrow(data@data)))   
  } else(!all(siteLabelCriteria %in% colnames(data@siteInfo))){
    stop("Not all siteLabelCriteria are columns in data@siteInfo")
  }
  ## Check if lineage table contains all edges 
}