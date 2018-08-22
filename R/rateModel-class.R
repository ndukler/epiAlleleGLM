#' Class rateModel
#'
#' Class \code{rateModel} holds a rate model for a set of alleles
#'
#' @slot alleleData an allele data object
#' @slot edgeGroups a data.table with four columns: parent, child, edgeID, edgeGroup. EdgeIDs must match those of alleleData object.
#' @slot siteLabelCriteria a character vector containing the names of the columns in alleleData@siteInfo to use to label sites
#' @slot paramEnviron A locked environment that contains three elements params, paramIndex, and fixed
#'
#' @name rateModel-class
#' @rdname rateModel-class 
#' @include rateModelValidityCheck.R
#' @importClassesFrom data.table data.table
#' @exportClass rateModel
methods::setClass("rateModel", slots=c(alleleData = "environment",edgeGroups="data.table",
                                              siteLabelCriteria="character",siteLabels="data.table",
                                              params="numeric",paramIndex="data.table",fixed="logical"),
                  validity = rateModelValidityCheck)

#' rateModel
#'
#' Contructs an object that holds a rate model for a set of alleles
#' @param data An alleleData object
#' @param siteLabelCriteria A charecter vector of the columns in siteData contained in alleleData label genomic regions
#' @param lineageTable A table with columns parent,child,edgeID, and rateGroup, where rateGroup is used to specify how the rates are tied between the branches
#' @param rate Starting value for rate parameter. (Default: 0.1)
#' @param pi Starting frequencies for allele stationary distribution. Must be the same length as the number of alleles. Default is a uniform probability.
#' @name rateModel
#' @return an rateModel object
#' @examples
#'
#' @export
rateModel <- function(data,siteLabelCriteria=NULL,lineageTable=NULL,rate=NULL,pi=NULL){
  ## ** Validity checks and default setting** ##
  ## Check that data is of class alleleData
  if(class(data)!="alleleData"){
    stop("data must be of class alleleData")
  }
  ## Check that all siteLabelCriteria elements are columns in alleleData@siteInfo
  if(is.null(siteLabelCriteria)){
    siteLabels=factor(rep("1",nrow(data@data)))
  } else if(!all(siteLabelCriteria %in% colnames(data@siteInfo))){
    diff=setdiff(siteLabelCriteria,colnames(data@siteInfo))
    diff=paste0("\'",diff,"\'")
    stop(paste("Labels",diff,"are not columns in data@siteInfo",collapse = ", "))
  }
  ## Checks for lineage table validity
  if(!is.null(lineageTable)){
    if(!is.data.frame(lineageTable)){
      stop("Lineage table must be a data.frame or data.table")
    } else if(!all(c("edgeID","edgeGroup") %in% colnames(lineageTable))){
      stop("Lineage table must contain the columns \'edgeID\' and \'edgeGroup\'")      
    } else if(!setequal(getEdgeTable(data)$edgeID,lineageTable$edgeID)){
      stop("EdgeIDs in the alleleData object and the supplied lineageTable do not match. Run getEdgeTable(data) to view alleleData edgeIds.")  
    } else if(any(table(lineageTable$edgeID)>1)){
      stop("Duplicated edgeIDs in lineageTable")
    }
  } else {
    ## Create default lineage table
    lineageTable=getEdgeTable(data)
    lineageTable[,edgeGroup:="e1"]
  }
  ## If pi is not NULL, check that it is the same length as the number of alleles, otherwise set to default
  if(!is.null(pi)){
    if(length(pi) != data@nAlleles){
      stop("The length of the pi vector must be the same as the number of alleles")
    }
  } else {
    pi=rep(1,data@nAlleles)
  }
  ## If rate is not null, check that it is either length 1 or the same length as the number of edge groups, 
  ## otherwise set to default
  if(!is.null(rate)){
    if(length(rate) != 1){
      stop("Length of the rate vector must be one")
    } else if(any(rate<=0)){
      stop("All rates must be greater than zero")
    } 
  } else {
    rate=0.1    
  }
  
  ## ** Intermediate reformating and computation ** ##
  ## Create environment to hold parameter values and associated indices, etc.
  adEnviron=new.env()
  adEnviron$alleleData=data
  
  ## Renormalize pi
  pi=pi/sum(pi)
  
  ## Standardize format for lineageTable
  lineageTable[,c("parent","child"):=as.list(as.integer(unlist(strsplit(edgeID,split = "-")))),by="edgeID"]
  data.table::setcolorder(lineageTable,c("parent","child","edgeID","edgeGroup"))
  data.table::setkeyv(x = lineageTable,cols = c("edgeID"))
  
  ## Create labels for each site
  siteLabels=getSiteInfo(data)[,..siteLabelCriteria][,.(siteLabel=do.call(paste, c(.SD, sep = "_")))]
  siteLabels[,siteLabel:=factor(siteLabel)]
  siteLabels[,index:=1:nrow(siteLabels)]
  data.table::setkeyv(siteLabels,c("index","siteLabel"))
  
  ## Create parameter index
  paramIndex=data.table::data.table(expand.grid(edgeGroup=as.character(unique(lineageTable$edgeGroup)),
                                                   siteLabel=as.character(levels(siteLabels$siteLabel)),stringsAsFactors = FALSE))
  paramIndex=paramIndex[order(siteLabel,edgeGroup)] ## Order this way so that parameter indicies are grouped by site
  paramIndex[,rateIndex:=1:nrow(paramIndex)]
  paramIndex[,piIndex:=(max(rateIndex)+1)+((as.numeric(as.factor(siteLabel))-1)*3)]
  data.table::setkeyv(x = paramIndex,cols = c("edgeGroup","siteLabel"))
  
  ## Build the parameter vector
  params=numeric(max(paramIndex$piIndex)+data@nAlleles-1)
  params[1:max(paramIndex$rateIndex)]=rate
  params[(max(paramIndex$rateIndex)+1):(max(paramIndex$piIndex)+data@nAlleles-1)]=rep(pi,length(unique(paramIndex$siteLabel)))
  
  ## Set fixed vector to default to FALSE
  fixed=logical(length(params))
  
  ## ** Object construction ** ##
  methods::new("rateModel",alleleData=adEnviron,edgeGroups=lineageTable,siteLabelCriteria=siteLabelCriteria,
               siteLabels=siteLabels,params=params,paramIndex=paramIndex,fixed=fixed)
}
