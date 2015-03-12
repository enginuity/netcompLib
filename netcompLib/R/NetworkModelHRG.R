# Defines the NetworkModelHRG class
setClass("NetworkModelHRG", representation(assign = "numeric", probmat = "matrix"), contains = "NetworkModel")

#' Constructor for HRG network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelHRG object
#' 
#' @export
#' 
NetworkModelHRG = function(Nnodes = 10, model_param = set_model_param()) {
  ## TODO: - fill this in eventually 
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelHRG objects, this returns "tree"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'tree'
#' 
#' @export
#' 
getNetType.NetworkModelHRG = function(NetM) { "tree" }
setMethod("getNetType", signature(NetM = "NetworkModelHRG"), getNetType.NetworkModelHRG)


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelHRG = function(NetM) {
  ## TODO - fix this
  
  return(res)
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelHRG"), getEdgeProbMat.NetworkModelHRG)

