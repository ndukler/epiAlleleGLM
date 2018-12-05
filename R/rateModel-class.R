#' Class rateModel
#'
#' Class \code{rateModel} holds a rate model for a set of alleles
#'
#' @slot alleleData an allele data object in a locked environment
#' @slot edgeGroups a data.table with four columns: parent, child, edgeID, edgeGroup. EdgeIDs must match those of alleleData object.
#' @slot rateFormula A formula that uses the variables in the allele data object to compute turnover rate
#' @slot piFormula A formula that uses the variables in the allele data object to compute pi
#' @slot rateDM the design matrix for the rate
#' @slot piDM the design matrix for pi
#' @slot params a vector of parameter values
#' @slot rateIndex a parameter index for the rate coefficients
#' @slot piIndex a parameter index for the pi coefficients
#' @slot fixed a logical vector indicating which variables are fixed
#'
#' @name rateModel-class
#' @rdname rateModel-class 
#' @include rateModelValidityCheck.R
#' @importClassesFrom data.table data.table
#' @importClassesFrom Matrix Matrix
#' @exportClass rateModel
methods::setClass("rateModel", slots=c(alleleData = "environment",edgeGroups="data.table",
                                       rateFormula="formula",piFormula="formula",rateDM="Matrix",
                                       piDM="Matrix",params="numeric",rateIndex="data.table",piIndex="data.table",
                                       fixed="logical"),
                  validity = rateModelValidityCheck)

#' rateModel
#'
#' Contructs an object that holds a rate model for a set of alleles
#' @param data An alleleData object
#' @param rateFormula A formula that uses the variables in the allele data object to compute turnover rate
#' @param piFormula A formula that uses the variables in the allele data object to compute pi (if NULL uses the same as rateFormula)
#' @param lineageTable A table with columns parent,child,edgeID, and rateGroup, where rateGroup is used to specify how the rates are tied between the branches
#' @name rateModel
#' @return an rateModel object
#' @examples
#'
#' @export
rateModel <- function(data,rateFormula,piFormula=NULL,lineageTable=NULL){
  ## ** Validity checks and default setting** ##
  ## Check that data is of class alleleData
  if(class(data)!="alleleData"){
    stop("data must be of class alleleData")
  }
  ## Check that a rate formula is specified as a formula
  if(class(rateFormula)!="formula"){
    stop("rateFormula must be a formula")
  }
  ## If piFormula is NULL set equal to rateFormula
  if(is.null(piFormula)){
    write("piFormula is not specified, using same formula as rateFormula...")
    piFormula=rateFormula
  }
  ## Check that all coavariates specified in rateFormula are contained in siteInfo
  if(!all(all.vars(rateFormula) %in% colnames(getSiteInfo(data)))){
    stop("Some of the covariates in the rateFormula are not present in the siteInfo of the alleleData object.")
  }
  ## Check that all coavariates specified in piFormula are contained in siteInfo
  if(!all(all.vars(piFormula) %in% colnames(getSiteInfo(data)))){
    stop("Some of the covariates in the piFormula are not present in the siteInfo of the alleleData object.")
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

  ## ** Intermediate reformating and computation ** ##
  ## Create environment to hold parameter values and associated indices, etc.
  adEnviron=new.env()
  adEnviron$alleleData=data
  
  ## Create the design matrix for the rate variable
  rateDM=Matrix::Matrix(model.matrix(rateFormula,getSiteInfo(data)))
  ## Create pi design matrix
  if(isTRUE(all.equal.formula(rateFormula,piFormula))){
    piDM=rateDM
  } else {
    piDM=Matrix::Matrix(model.matrix(piFormula,getSiteInfo(data)))
  }
  
  ## Standardize format for lineageTable
  lineageTable[,c("parent","child"):=as.list(as.integer(unlist(strsplit(edgeID,split = "-")))),by="edgeID"]
  data.table::setcolorder(lineageTable,c("parent","child","edgeID","edgeGroup"))
  data.table::setkeyv(x = lineageTable,cols = c("edgeID"))
  
  ## Create parameter index
  rateIndex=data.table::data.table(expand.grid(edgeGroup=as.character(unique(lineageTable$edgeGroup))),
                                 column=1:ncol(rateDM),stringsAsFactors = FALSE)
  piIndex=data.table::data.table(expand.grid(allele=2:data@nAlleles,
                               column=1:ncol(piDM),stringsAsFactors = FALSE))
  ## Add covariate names
  rateIndex[,name:=colnames(rateDM)[column]]
  piIndex[,name:=colnames(piDM)[column]]
  ## Order in index so they match the column order of the design matricies
  rateIndex[,index:=1:nrow(rateIndex)]
  piIndex[,index:=(nrow(rateIndex)+1):(nrow(rateIndex)+nrow(piIndex))]
  ## Set index jeys  
  data.table::setkeyv(x = rateIndex,cols = c("edgeGroup","column"))
  data.table::setkeyv(x = piIndex,cols = c("allele","column"))
  
  ## Build the parameter vector
  params=rep(1,nrow(rateIndex)+nrow(piIndex))
  
  ## Set fixed vector to default to FALSE
  fixed=logical(length(params))
  
  ## ** Object construction ** ##
  methods::new("rateModel",alleleData=adEnviron,edgeGroups=lineageTable,rateFormula=rateFormula,
              piFormula=piFormula,rateDM=rateDM,piDM=piDM,params=params,rateIndex=rateIndex,piIndex=piIndex,
              fixed=fixed)
}
