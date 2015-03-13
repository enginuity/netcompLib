# Defines the NetworkModelRND class
setClass("NetworkModelRND", representation(assign = "numeric", probmat = "matrix"), contains = "NetworkModel")

#' Constructor for RND network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkModelRND = function(Nnodes = 10, model_param = set_model_param()) {
  ## TODO: [Implement] - fill this in eventually 
  stop("Not implemented yet")
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelRND objects, this returns "random"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'random'
#' 
#' @export
#' 
getNetType.NetworkModelRND = function(NetM) { "random" }
setMethod("getNetType", signature(NetM = "NetworkModelRND"), getNetType.NetworkModelRND)


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelRND = function(NetM) {
  ## TODO: [Implement]
  
  stop("Not implemented yet")
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelRND"), getEdgeProbMat.NetworkModelRND)

