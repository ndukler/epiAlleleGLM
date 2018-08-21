methods::setGeneric("tieRates", function(obj,siteLabelA,edgeGroupA,siteLabelB,edgeGroupB) {
  standardGeneric("tieRates")
})

#' Tie rates across siteLabels or edgeGroups 
#'
#' Sets the rate of the set of branches specified by edgeGroupB and siteLabelB equal to that of A
#' @param obj a rateModel object
#' @param siteLabelA site label for rate A
#' @param edgeGroupA edge group for rate A
#' @param siteLabelB site label for rate B
#' @param edgeGroupB edge group for rate B
#' @name tieRates
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("tieRates", signature(obj = "rateModel"), function(obj,siteLabelA,edgeGroupA,siteLabelB,edgeGroupB) {
  if(!is.vector(siteLabelA) || length(siteLabelA)!=1){
    stop("siteLabelA must be a vector of length one")
  }
  if(!is.vector(siteLabelB) || length(siteLabelB)!=1){
    stop("siteLabelB must be a vector of length one")
  }
  if(!is.vector(edgeGroupA) || length(edgeGroupA)!=1){
    stop("edgeGroupA must be a vector of length one")
  }
  if(!is.vector(edgeGroupB) || length(edgeGroupB)!=1){
    stop("edgeGroupB must be a vector of length one")
  }
  ## Get parameter index
  pI=getParamIndex(obj)
  ## Check that queries are valid
  if(nrow(pI[.(edgeGroupA,siteLabelB),,nomatch=0])==0){
    stop("Invalid query B")
  }
  if(nrow(pI[.(edgeGroupA,siteLabelA),,nomatch=0])==0){
    stop("Invalid query A")
  }
  ## tie parameters, then recompute indices
  pI[,indexLevel:=rateIndex]
  pI[.(edgeGroupB,siteLabelB),indexLevel:=pI[.(edgeGroupA,siteLabelA),.(rateIndex)]] ## Set indexLevel of B so it points to A
  pI[,indexLevel:=as.numeric(as.factor(indexLevel))]
  ## Update parameter vector
  indB=getRateIndex(obj,edgeGroup = as.character(edgeGroupB),siteLabel = as.character(siteLabelB))
  unlockBinding("params", obj@paramEnviron)
  obj@paramEnviron$params[pI$rateIndex]=obj@paramEnviron$params[pI$indexLevel] ## Set parameters so that tied parameters are equal
  obj@paramEnviron$params=obj@paramEnviron$params[-indB] ## Remove B index param
  lockBinding("params", obj@paramEnviron)
  ## Update paramIndex and store
  pI[,rateIndex:=indexLevel]
  pI[,piIndex:=piIndex-1]
  pI[,indexLevel:=NULL]
  unlockBinding("paramIndex", obj@paramEnviron)
  obj@paramEnviron$paramIndex=pI
  lockBinding("paramIndex", obj@paramEnviron)
})