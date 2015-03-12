# This file defines the NetworkModel class.

## TODO: Add validity checks?



# Defines the NetworkModel Class and its subclasses -----------------------
setClass("NetworkModel", representation(Nnodes = "numeric"))
setClass("NetworkModelLSM", contains = "NetworkModel")
setClass("NetworkModelHRG", contains = "NetworkModel")
setClass("NetworkModelRND", contains = "NetworkModel")

# 
# setMethod("getNetType", signature(NetM = "NetworkModelLSM"), function(NetM) { "latent" })
# getNetType(new("NetworkModelLSM"))


#' Instantiates an object of class NetworkModel
#' 
#' @return Object of class NetworkModel
#' 
#' @export
#' 
NetworkModel = function() {
  # creates a default networkmodel object -- this is the default constructor, but probably should never be used. the specific ones for a specific model should be used. 
  # maybe want to make this eventually call the network generation methods
  
  ## currently this does nothing but return a lame object, that satisfies the generic methods (although the generic methods would not do much?)
  
  NetM = new("NetworkModel", Nnodes = 10)
}


# Define Generic Methods --------------------------------------------------


### New method to extract the number of nodes from a network
setGeneric("getNnodes", function(NetM) standardGeneric("getNnodes"))

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNnodes)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNnodes = function(NetM) { NULL }

#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkModel
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkModel = function(NetM) {
  # NetM should be object of type NetworkModel
  return(NetM@Nnodes)
}
setMethod("getNnodes", signature("NetworkModel"), getNnodes.NetworkModel)








# Define Default print/summary methods ------------------------------------





## Methods for network-model specific information
## These are the functions that need to be defined for each specific network model. 


# Define generic functions that are model-specific ------------------------


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
getNetType = function(NetM) { NULL }




setGeneric("getEdgeProbMat", function(NetM) standardGeneric("getEdgeProbMat"))

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getEdgeProbMat)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getEdgeProbMat = function(NetM) { return(NULL) }




setGeneric("sampleNetwork", function(NetM, Nobs = 1, ...) standardGeneric("sampleNetwork"))

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sampleNetwork)
#' <What does this function do>
#' 
#' @param NetM temp
#' @param Nobs temp
#' @param ... temp
#' 
#' @return temp
#' 
#' @export
#' 
sampleNetwork = function(NetM, Nobs = 1, ...) { return(NULL) }


# Define generic null-model values so this doesn't crash ------------------

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkModel)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkModel)
#' <What does this function do>
#' 
#' @param netM temp
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkModel = function(NetM) { "none" }
setMethod("getNetType", signature(NetM = "NetworkModel"), getNetType.NetworkModel)

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getEdgeProbMat.NetworkModel)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getEdgeProbMat.NetworkModel = function(NetM) { NULL }
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModel"), getEdgeProbMat.NetworkModel)

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sampleNetwork.NetworkModel)
#' <What does this function do>
#' 
#' @param NetM temp
#' @param Nobs temp
#' @param ... temp
#' 
#' @return temp
#' 
#' @export
#' 
sampleNetwork.NetworkModel = function(NetM, Nobs = 1, ...) { NULL }
setMethod("sampleNetwork", signature = (NetM = "NetworkModel"), sampleNetwork.NetworkModel)


