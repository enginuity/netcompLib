##@S Function to extract the type of network model

## TODO: [Fully Documented] (remove this marking eventually)

setGeneric("getNetType", function(NetM) standardGeneric("getNetType"))

#' Extract the type of network model
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



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModel = function(NetM) {
  "none" 
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStruct = function(NetM) {
  "none" 
}



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelPair = function(NetM) {
  return(c(getNetType(NetM@m1), getNetType(NetM@m2)))
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelSBM = function(NetM) { "block" }


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelLSM = function(NetM) { "latent" }


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStructList = function(NetM) {
  return(sapply(NetM@models, getNetType))
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelHRG = function(NetM) { "tree" }



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelRND = function(NetM) { "random" }



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStructRND = function(NetM) {
  return("random")
}



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStructHRG = function(NetM) {
  return("tree")
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
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

