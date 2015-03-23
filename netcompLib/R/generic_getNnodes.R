
setGeneric("getNnodes", function(NetM) standardGeneric("getNnodes"))

#' Generic Function -- Extracts the number of nodes in the network model
#' 
#' @param NetM NetworkModel object
#' 
#' @return Numeric -- Number of nodes in network model
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


#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkStruct
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkStruct = function(NetM) {
  # NetM should be object of type NetworkStruct
  return(NetM@Nnodes)
}


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


# setMethod ---------------------------------------------------------------

setMethod("getNnodes", signature("NetworkModel"), getNnodes.NetworkModel)
setMethod("getNnodes", signature("NetworkModelPair"), getNnodes.NetworkModelPair)
setMethod("getNnodes", signature("NetworkStruct"), getNnodes.NetworkStruct)

