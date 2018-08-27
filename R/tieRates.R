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
#' @return A semi-deep copy of the rateModel object with tied parameters (points to same environment that held the data)
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
  pI=data.table::copy(getParamIndex(obj))
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
  ## Create updated parameter vector
  indB=getRateIndex(obj,edgeGroup = as.character(edgeGroupB),siteLabel = as.character(siteLabelB))
  params=obj@params
  params[pI$rateIndex]=params[pI$indexLevel] ## Set parameters so that tied parameters are equal
  params=params[-indB] ## Remove B index param
  ## Update paramIndex and store
  pI[,rateIndex:=indexLevel]
  pI[,piIndex:=piIndex-1]
  pI[,indexLevel:=NULL]
  ## Return rebuilt object
  return(methods::new("rateModel",alleleData=obj@alleleData,edgeGroups=obj@edgeGroups,siteLabelCriteria=obj@siteLabelCriteria,
               siteLabels=obj@siteLabels,params=params,paramIndex=pI,fixed=obj@fixed[-indB]))
})