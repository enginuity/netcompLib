##@S Function to extract the number of nodes in the network

setGeneric("getNnodes", function(Net) standardGeneric("getNnodes"))

#' Extracts the number of nodes in a network model
#' 
#' @param Net [\code{\link{NetworkModel}} OR \code{\link{NetworkStruct}}] :: Network object
#' 
#' @return [int] :: Number nodes in network
#' 
#' @export
#' 
getNnodes = function(Net) { NULL }


getNnodes.NetworkModel = function(Net) {
  return(Net@Nnodes)
}


getNnodes.NetworkStruct = function(Net) {
  return(Net@Nnodes)
}


getNnodes.NetworkModelPair = function(Net) {
  return(Net@Nnodes)
}


# setMethod ---------------------------------------------------------------
setMethod("getNnodes", signature("NetworkModel"), getNnodes.NetworkModel)
setMethod("getNnodes", signature("NetworkModelPair"), getNnodes.NetworkModelPair)
setMethod("getNnodes", signature("NetworkStruct"), getNnodes.NetworkStruct)

