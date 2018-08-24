methods::setGeneric("fit", function(obj,scale) {
  standardGeneric("fit")
})

#' Fits rate model
#'
#' Fits rate model object and returns fitted model
#' @param obj rateModel
#' @param scale a scale factor to apply to log-likelihood, defaults to -1/nsites
#' @name fit
#' @return a list incuding a rateModel object and information about the optimization
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("fit", signature(obj = "rateModel"), function(obj,scale=NULL) {
  ## scale defaults to -1/nsites
  if(is.null(scale)){
    sca=-nrow(getAlleleData(obj)@data)
  }
  pI=getParamIndex(obj)
  x=getParams(obj)
  rateEnd=max(pI$rateIndex)
  ## Convert parameters to stick breaking parameterization
  y=c(x[1:rateEnd],multiProbToStick(x[(rateEnd+1):length(x)],getAlleleData(obj)@nAlleles))
  ## Set default box constraints, between 0 and 1
  ub=c(rep(Inf,rateEnd),rep(1,length(y)-rateEnd))
  lb=numeric(length(y))
  ## NEED TO FIX SO THAT RATE PARAMETERS CAN BE FIXED INDIVIDUALLY BUT PI IS ALL OR NOTHING PER SITE!!!!!
  ## Where the parameter values are fixed, set lb=ub=value
  ## ub[obj@fixed]=x[obj@fixed]
  ## lb[obj@fixed]=x[obj@fixed]
  optMod=optim(y,fn = logLikelihood,obj=obj,stickParams=TRUE,lower = lb,upper = ub,method="L-BFGS-B",
               control = list(fnscale=sca))
  probs=multiStickToProb(optMod$par[(rateEnd+1):length(optMod$par)],width=getAlleleData(obj)@nAlleles-1)
  setParamValue(obj,i = 1:length(obj@params),value = c(optMod$par[1:rateEnd],probs))
  return(with(optMod,list(model=obj,value=value,counts=counts,convergence=convergence,message=message)))
})

