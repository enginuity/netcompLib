##@S Function to extract the number of nodes in the network

## TODO: [Fully Documented] (remove this marking eventually)

setGeneric("getNnodes", function(NetM) standardGeneric("getNnodes"))

#' Extracts the number of nodes in a network model
#' 
#' @param NetM NetworkModel object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes = function(NetM) { NULL }


#' Extracts the number of nodes in a network model
#' 
#' @param NetM NetworkModel object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes.NetworkModel = function(NetM) {
  # NetM should be object of type NetworkModel
  return(NetM@Nnodes)
}


#' Extracts the number of nodes in a network model
#' 
#' @param NetM NetworkStruct object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes.NetworkStruct = function(NetM) {
  # NetM should be object of type NetworkStruct
  return(NetM@Nnodes)
}



#' Extracts the number of nodes in a network model
#' 
#' @param NetM NetworkModel object
#' 
#' @return Numeric -- Number of nodes in network model
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

