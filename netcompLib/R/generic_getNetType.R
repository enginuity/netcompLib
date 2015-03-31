##@S Function to extract the type of network model

## TODO: [Fully Documented] (remove this marking eventually)
## TODO: [Fix Documentation] Documentation for parameter 'Net', 'NetM', 'NetS' should be updated (and made consistent throughout the entire project)

setGeneric("getNetType", function(Net) standardGeneric("getNetType"))

#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType = function(Net) { 
  NULL 
}



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModel = function(Net) {
  "none" 
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStruct = function(Net) {
  "none" 
}



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelPair = function(Net) {
  return(c(getNetType(Net@m1), getNetType(Net@m2)))
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelSBM = function(Net) { "block" }


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelLSM = function(Net) { "latent" }


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStructList = function(Net) {
  return(sapply(Net@models, getNetType))
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelHRG = function(Net) { "tree" }



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkModelRND = function(Net) { "random" }



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStructRND = function(Net) {
  return("random")
}



#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType.NetworkStructHRG = function(Net) {
  return("tree")
}


#' Extract the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param Net Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
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

