methods::setGeneric("fit", function(obj,scale) {
  standardGeneric("fit")
})

#' Fits rate model
#'
#' Fits rate model object and returns fitted model
#' @param obj rateModel
#' @param scale a scale factor to apply to log-likelihood, defaults to -1/nsites
#' @name fit
#' @return rateModel object
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("fit", signature(obj = "rateModel"), function(obj,scale=NULL) {
  ## scale defaults to -1/nsites
  if(is.null(scale)){
    scale=-1/nrow(getAlleleData(obj)@data)
  }
  pI=getParamIndex(obj)
  x=getParams(obj)
  ## Set default box constraints, between 0 and 1
  ub=c(rep(Inf,max(pI$rateIndex)),rep(1,length(x)-max(pI$rateIndex)))
  lb=numeric(length(x))
  ## Where the parameter values are fixed, set lb=ub=value
  ub[obj@fixed]=x[obj@fixed]
  lb[obj@fixed]=x[obj@fixed]
  optMod=optim(x,fn =logLikelihood,obj=obj,method="L-BFGS-B",control = list(fnscale=scale))
  
  
})

