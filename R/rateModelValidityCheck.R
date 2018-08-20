#' @description  Validity checker for alleleData-class
rateModelValidityCheck <- function(object){
  errors=c()
  if(!setequal(colnames(object@edgeGroups),c("parent","child","edgeID","edgeGroup"))){
    errors=c(errors,"Colnames of edgeGroups must be \'parent\', \'child\', \'edgeID\',\'edgeGroup\'")
  }
  if(!setequal(data.table::key(object@edgeGroups),c("parent","child"))){
    errors=c(errors,"Keys of edgeGroups must be \'parent\', \'child\'")
  }
  if(!setequal(colnames(object@siteLabels),"siteLabel")){
    errors=c(errors,"object@siteLabels must have only one column, \'siteLabel\'")  
  }
  if(!setequal(names(object@paramEnviron),c("params","paramIndex","fixed"))){
    errors=c(errors,"object@paramEnviron must have three slots, \'params\', \'paramIndex\', \'fixed\'") 
  } else {
    ## Checks for paramIndex in environment
    if(!data.table::is.data.table(object@paramEnviron$paramIndex)){
      errors=c(errors,"object@paramEnviron$paramIndex must be a data.table") 
    } else if(!setequal(colnames(object@paramEnviron$paramIndex),c("edgeGroup","siteLabel","rateIndex","piIndex"))){
      errors=c(errors,"object@paramEnviron$paramIndex must have columns \'edgeGroup\', \'siteLabel\', \'rateIndex\', \'piIndex\'")
    } else if(!setequal(data.table::key(object@paramEnviron$paramIndex),c("edgeGroup","siteLabel"))){
      errors=c(errors,"object@paramEnviron$paramIndex must have keys \'edgeGroup\', \'siteLabel\'")
    }
    ## Check for the params
    expLen= max(object@paramEnviron$paramIndex$piIndex)+object@alleleData@nAlleles-1 ## expected length of parameter vec
    if(!is.numeric(object@paramEnviron$params) || length(object@paramEnviron$params) != expLen){
      errors=c(errors,paste("object@paramEnviron$params must be a numeric vector of length",expLen)) 
    }
    ## Checks for the fixed vector
    if(!is.logical(object@paramEnviron$fixed) || length(object@paramEnviron$fixed) != expLen){
      errors=c(errors,paste("object@paramEnviron$fixed must be a logical vector of length",expLen)) 
    }
  }
  ## Lock the environment and all the bindings in it if all tests passed
  if (length(errors) == 0){
    lockEnvironment(env = object@paramEnviron,bindings = TRUE)
  }
  ## Return errors if there are any
  if (length(errors) == 0) TRUE else errors
}