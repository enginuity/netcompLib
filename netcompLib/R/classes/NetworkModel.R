# This file defines the NetworkModel class.

## TODO: Add validity checks?


# Defines the NetworkModel Class and its subclasses -----------------------
setClass("NetworkModel", representation(Nnodes = "numeric", Type = "character"))
setClass("SBMNetworkModel", contains = "NetworkModel")
setClass("LSMNetworkModel", contains = "NetworkModel")
setClass("HRGNetworkModel", contains = "NetworkModel")


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
  
  NetM = new("NetworkModel", Nnodes = 10, Type = NA)
}

# Define Generic Methods --------------------------------------------------


### New method to extract the number of nodes from a network
setGeneric("getNnodes", function(NetM) standardGeneric("getNnodes"))

getNnodes.NetworkModel = function(NetM) {
  # NetM should be object of type NetworkModel
  return(NetM@Nnodes)
}

setMethod("getNnodes", signature("NetworkModel"), getNnodes.NetworkModel)




### New method to extract the type of network model
setGeneric("getNetType", function(NetM) standardGeneric("getNetType"))

getNetType.NetworkModel = function(NetM) {
  return(NetM@Type)
}

setMethod("getNetType", signature("NetworkModel"), getNetType.NetworkModel)
















# 
# ## test 
# 
# 
# test = new("NetworkModel", Nnodes = 20)
# test = new("SBMNetworkModel", Nnodes = 3)
# test
# getNnodes(test)
# 
# showMethods(classes = "NetworkModel")
# 
# getNodes
# 
# 
# 
