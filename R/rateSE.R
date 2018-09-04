methods::setGeneric("rateSE", function(obj,alpha=0.95) {
  standardGeneric("rateSE")
})

#' Estimates the standard error of the rate
#'
#' Estimates the standard error for all rates in the model using the hessian
#' @param obj rateModel
#' @name rateSE
#' @return a data.table with the parameter values and standard errors
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("rateSE", signature(obj = "rateModel"), function(obj,alpha=NULL) {
  ## Get the parameter index
  parInd=getParamIndex(obj)
  ## Get parameter values
  paramVals=getParams(obj)[1:max(parInd$rateIndex)]
  ## Compute the hessian
  hess=numDeriv::hessian(func=logLikelihood,x=getParams(obj),obj=obj,stickParams=FALSE)
  ## Subset hessian to only include the rates
  hess=hess[1:max(parInd$rateIndex),1:max(parInd$rateIndex)]
  ## Invert the hessian
  fisherInfo<-solve(-hess)
  stdErr<-sqrt(diag(fisherInfo))
  ## Construct the out table to return
  outTab=parInd[order(rateIndex),.(edgeGroup,siteLabel,rateIndex)]
  outTab[,value:=paramVals]
  outTab[,se:=stdErr]
  return(outTab)
})