#' @description  Validity checker for alleleData-class
rateModelValidityCheck <- function(object){
  errors=c()
  if(class(object@alleleData) != "alleleData"){
    errors= c(errors,"alleleData slot must be of class alleleData")
  }
  
  ## Check that there are the same number of rows in the siteInfo data.frame as there are rows in the data
  if (length(errors) == 0) TRUE else errors
}