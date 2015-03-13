# Defines the NetworkStructSBM class
setClass("NetworkStructSBM", representation(groups = "numeric", counts = "numeric", expand = "list"), contains = "NetworkStruct")

#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelSBM object
#' 
#' @export
#' 
NetworkStructSBM = function(Nnodes = 10, model_param = set_model_param()) {
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelSBM(Nnodes = Nnodes, model_param = model_param)
  return(extractStruct(NetM))
}


