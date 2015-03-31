##@S Function to extract the number of nodes in the network

## TODO: [Fully Documented] (remove this marking eventually)

setGeneric("getNnodes", function(Net) standardGeneric("getNnodes"))

#' Extracts the number of nodes in a network model
#' 
#' @param Net NetworkModel or NetworkStruct object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes = function(Net) { NULL }


#' Extracts the number of nodes in a network model
#' 
#' @param Net NetworkModel or NetworkStruct object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes.NetworkModel = function(Net) {
  # Net should be object of type NetworkModel
  return(Net@Nnodes)
}


#' Extracts the number of nodes in a network model
#' 
#' @param Net NetworkModel or NetworkStruct object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes.NetworkStruct = function(Net) {
  # Net should be object of type NetworkStruct
  return(Net@Nnodes)
}



#' Extracts the number of nodes in a network model
#' 
#' @param Net NetworkModel or NetworkStruct object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes.NetworkModelPair = function(Net) {
  # Net should be object of type NetworkModelPair
  return(Net@Nnodes)
}


# setMethod ---------------------------------------------------------------

setMethod("getNnodes", signature("NetworkModel"), getNnodes.NetworkModel)
setMethod("getNnodes", signature("NetworkModelPair"), getNnodes.NetworkModelPair)
setMethod("getNnodes", signature("NetworkStruct"), getNnodes.NetworkStruct)

