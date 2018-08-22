#' @description  Validity checker for alleleData-class
rateModelValidityCheck <- function(object){
  errors=c()
  if(!setequal(names(object@alleleData),c("alleleData"))){
    errors=c(errors,"object@alleleData must have one binding, \'alleleData\'") 
  } else if(class(object@alleleData$alleleData)!="alleleData"){
    errors=c(errors,"object@alleleData$alleleData is not of class allele Data")
  }
  if(!setequal(colnames(object@edgeGroups),c("parent","child","edgeID","edgeGroup"))){
    errors=c(errors,"Colnames of edgeGroups must be \'parent\', \'child\', \'edgeID\',\'edgeGroup\'")
  }
  if(!setequal(data.table::key(object@edgeGroups),c("edgeID"))){
    errors=c(errors,"Key of edgeGroups must be \'edgeID\'")
  }
  if(!setequal(colnames(object@siteLabels),c("siteLabel","index"))){
    errors=c(errors,"object@siteLabels must have two columns, \'siteLabel\' and \'index\'")  
  }

  ## Checks for paramIndex 
  if(!setequal(colnames(object@paramIndex),c("edgeGroup","siteLabel","rateIndex","piIndex"))){
    errors=c(errors,"object@paramIndex must have columns \'edgeGroup\', \'siteLabel\', \'rateIndex\', \'piIndex\'")
  } else if(!setequal(data.table::key(object@paramIndex),c("edgeGroup","siteLabel"))){
    errors=c(errors,"object@paramIndex must have keys \'edgeGroup\', \'siteLabel\'")
  }
  ## Check for the params
  expLen= max(object@paramIndex$piIndex)+object@alleleData$alleleData@nAlleles-1 ## expected length of parameter vec
  if(!is.numeric(object@params) || length(object@params) != expLen){
    errors=c(errors,paste("object@params must be a numeric vector of length",expLen)) 
  }
  ## Checks for the fixed vector
  if(!is.logical(object@fixed) || length(object@fixed) != expLen){
    errors=c(errors,paste("object@fixed must be a logical vector of length",expLen)) 
  }

  ## Lock the environment and all the bindings in it if all tests passed
  if (length(errors) == 0){
    lockEnvironment(env = object@alleleData,bindings = TRUE)
  }
  ## Return errors if there are any
  if (length(errors) == 0) TRUE else errors
}