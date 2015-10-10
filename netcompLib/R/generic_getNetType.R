##@S Function to extract the type of network model

setGeneric("getNetType", function(Net) standardGeneric("getNetType"))

#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net [\code{\link{NetworkModel}} OR \code{\link{NetworkStruct}}] :: Network object
#' 
#' @return [char] :: type of network object
#' 
#' @export
#' 
getNetType = function(Net) { 
  NULL 
}


getNetType.NetworkModel = function(Net) {
  "none" 
}


getNetType.NetworkStruct = function(Net) {
  "none" 
}


getNetType.NetworkModelPair = function(Net) {
  return(c(getNetType(Net@m1), getNetType(Net@m2)))
}


getNetType.NetworkModelSBM = function(Net) { "block" }


getNetType.NetworkModelLSM = function(Net) { "latent" }


getNetType.NetworkStructList = function(Net) {
  return(sapply(Net@models, getNetType))
}


getNetType.NetworkModelHRG = function(Net) {
  "tree" 
}


getNetType.NetworkModelRND = function(Net) {
  "random" 
}


getNetType.NetworkStructRND = function(Net) {
  return("random")
}


getNetType.NetworkStructHRG = function(Net) {
  return("tree")
}

 
getNetType.NetworkStructSBM = function(Net) {
  return("block")
}


# setMethod ---------------------------------------------------------------
setMethod("getNetType", signature(Net = "NetworkModel"), getNetType.NetworkModel)
setMethod("getNetType", signature(Net = "NetworkModelHRG"), getNetType.NetworkModelHRG)
setMethod("getNetType", signature(Net = "NetworkModelLSM"), getNetType.NetworkModelLSM)
setMethod("getNetType", signature(Net = "NetworkModelPair"), getNetType.NetworkModelPair)
setMethod("getNetType", signature(Net = "NetworkModelSBM"), getNetType.NetworkModelSBM)
setMethod("getNetType", signature(Net = "NetworkModelRND"), getNetType.NetworkModelRND)
setMethod("getNetType", signature(Net = "NetworkStruct"), getNetType.NetworkStruct)
setMethod("getNetType", signature(Net = "NetworkStructList"), getNetType.NetworkStructList)
setMethod("getNetType", signature(Net = "NetworkStructSBM"), getNetType.NetworkStructSBM)
setMethod("getNetType", signature(Net = "NetworkStructRND"), getNetType.NetworkStructRND)
setMethod("getNetType", signature(Net = "NetworkStructHRG"), getNetType.NetworkStructHRG)
