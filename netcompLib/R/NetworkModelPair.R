# This file defines the NetworkModelPair class.

## TODO: Add validity checks?

# Defines the NetworkModelPair Class and its subclasses -----------------------
setClass("NetworkModelPair", representation(m1 = "NetworkModel", m2 = "NetworkModel", is_null = "logical"), contains("NetworkModel"))

#' Instantiates an object of class NetworkModelPair
#' 
#' This is done by providing both network models. If m2 is not given and is_null is FALSE, there is an error. Otherwise, m2 will be ignored if is_null is TRUE. (ie null hypothesis means both models are the same, so m1)
#' There is also no check that m1 and m2 are on the same network size, but they should be. -- add this in...
#' 
#' @return Object of class NetworkModelPair
#' 
#' @export
#' 
NetworkModelPair = function(m1, m2 = NULL, is_null = FALSE) {
  if (is_null) {
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m1, is_null = TRUE)
  } else {
    if (is.null(m2)) { stop("---You must provide the second model if the null hypothesis is FALSE.---") }
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m2, is_null = FALSE)
  }
  return(netmp)
}


# Define Generic Methods --------------------------------------------------


#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkModelPair
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkModelPair = function(NetM) {
  # NetM should be object of type NetworkModelPair
  return(NetM@Nnodes)
}
setMethod("getNnodes", signature("NetworkModelPair"), getNnodes.NetworkModelPair)


#' Samples any number of observed networks from a given network model
#' 
#' Note -- this implementation requires the edge probability matrix. It's not really efficient (since it re-computes the edge probability matrix each time)
#' 
#' @param NetM Network Model object
#' @param Nobs Number of observed networks to sample
#' @param Nsim Number of simulations (number of times to do Nobs samples)
#' 
#' @return List of two lists/arrays of network models. 
#' 
#' @export
#' 
sampleNetwork.NetworkModelPair = function(NetM, Nobs = 1, Nsim = 1) { 
  return(list(sampleNetwork(NetM@m1), sampleNetwork(NetM@m2)))
}
setMethod("sampleNetwork", signature = (NetM = "NetworkModelPair"), sampleNetwork.NetworkModelPair)




# Define Default print/summary methods ------------------------------------





## Methods for network-model specific information
## These are the functions that need to be defined for each specific network model. 



# Define generic null-model values so this doesn't crash ------------------

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
setMethod("getNetType", signature(NetM = "NetworkModelPair"), getNetType.NetworkModelPair)

#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM Network Model object
#' 
#' @return NULL - no edge probability matrix for generic network model object
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelPair = function(NetM) {
  return(list(getEdgeProbMat(NetM@m1),getEdgeProbMat(NetM@m2)))
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelPair"), getEdgeProbMat.NetworkModelPair)

#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelPair = function(NetM) {
  netsl = new("NetworkStructList", Nnodes = getNnodes(NetM), 
              models = list(extractStruct(NetM@m1), extractStruct(NetM@m2)))
  return(netsl)
}
setMethod("extractStruct", signature = (NetM = "NetworkModelPair"), extractStruct.NetworkModelPair)
