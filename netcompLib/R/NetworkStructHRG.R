# Defines the NetworkStructHRG class
setClass("NetworkStructHRG", representation(tree_list = "list", expand = "list", counts = "numeric"), contains = "NetworkStruct")

#' Constructor for RND network structure
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
NetworkStructHRG = function(Nnodes = 10, model_param = set_model_param()) {
  
  stop("Not written")
  
  netm = new("NetworkModelRND", Nnodes = Nnodes, counts = sizes, ids = temp, 
             prob = runif(n = rnd_Ngroups, min = model_param$pmin, max = model_param$pmax))
  return(netm)
  
}


