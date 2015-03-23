
setGeneric("getNetType", function(NetM) standardGeneric("getNetType"))

#' Generic Function -- Extracts the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType = function(NetM) { 
  NULL 
}

#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkModel = function(NetM) {
  "none" 
}


#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkStruct = function(NetM) {
  "none" 
}



#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkModelPair = function(NetM) {
  return(c(getNetType(NetM@m1), getNetType(NetM@m2)))
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelSBM objects, this returns "block"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'block'
#' 
#' @export
#' 
getNetType.NetworkModelSBM = function(NetM) { "block" }


# Defines the NetworkModelLSM class

#' Returns the type of network model
#' 
#' Specifically for NetworkModelLSM objects, this returns "latent"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'latent'
#' 
#' @export
#' 
getNetType.NetworkModelLSM = function(NetM) { "latent" }


#' Extract the type of network model
#' 
#' @param NetM Object of class NetworkStructList
#' 
#' @return character vector of types
#' 
#' @export
#' 
getNetType.NetworkStructList = function(NetM) {
  return(sapply(NetM@models, getNetType))
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelHRG objects, this returns "tree"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'tree'
#' 
#' @export
#' 
getNetType.NetworkModelHRG = function(NetM) { "tree" }



#' Returns the type of network model
#' 
#' Specifically for NetworkModelRND objects, this returns "random"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'random'
#' 
#' @export
#' 
getNetType.NetworkModelRND = function(NetM) { "random" }



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructRND)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructRND = function(NetM) {
  return("random")
}

# Generic Function Definitions --------------------------------------------

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructHRG)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructHRG = function(NetM) {
  return("tree")
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructSBM)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructSBM = function(NetM) {
  return("block")
}


# setMethod ---------------------------------------------------------------

setMethod("getNetType", signature(NetM = "NetworkModel"), getNetType.NetworkModel)
setMethod("getNetType", signature(NetM = "NetworkModelHRG"), getNetType.NetworkModelHRG)
setMethod("getNetType", signature(NetM = "NetworkModelLSM"), getNetType.NetworkModelLSM)
setMethod("getNetType", signature(NetM = "NetworkModelPair"), getNetType.NetworkModelPair)
setMethod("getNetType", signature(NetM = "NetworkModelSBM"), getNetType.NetworkModelSBM)
setMethod("getNetType", signature(NetM = "NetworkModelRND"), getNetType.NetworkModelRND)
setMethod("getNetType", signature(NetM = "NetworkStruct"), getNetType.NetworkStruct)
setMethod("getNetType", signature(NetM = "NetworkStructList"), getNetType.NetworkStructList)
setMethod("getNetType", signature(NetM = "NetworkStructSBM"), getNetType.NetworkStructSBM)
setMethod("getNetType", signature(NetM = "NetworkStructRND"), getNetType.NetworkStructRND)
setMethod("getNetType", signature(NetM = "NetworkStructHRG"), getNetType.NetworkStructHRG)

