# Defines the NetworkModelLSM class
setClass("NetworkModelLSM", representation(assign = "numeric", probmat = "matrix"), contains = "NetworkModel")

#' Constructor for LSM network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelLSM object
#' 
#' @export
#' 
NetworkModelLSM = function(Nnodes = 10, model_param = set_model_param()) {
  ## TODO: - fill this in eventually 
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelLSM objects, this returns "latent"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'latent'
#' 
#' @export
#' 
getNetType.NetworkModelLSM = function(NetM) { "latent" }
setMethod("getNetType", signature(NetM = "NetworkModelLSM"), getNetType.NetworkModelLSM)


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelLSM = function(NetM) {
  ## TODO - fix this
  
  return(res)
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelLSM"), getEdgeProbMat.NetworkModelLSM)

