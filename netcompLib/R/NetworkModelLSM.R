# Defines the NetworkModelLSM class
setClass("NetworkModelLSM", representation(locs = "matrix", alpha = "numeric"), contains = "NetworkModel")

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
  
  K = model_param$latent_nclass
  D = model_param$latent_dim
  gen_norm = model_param$latent_isgennorm
  sd_center = model_param$latent_sdcenter
  ## TODO: [TEMP] REMOVE these variables. 
  
  # K,gen_norm = TRUE, D = 2, sd_center = 5) {
  
  group_assign = sample(1:K, size = Nnodes, replace = TRUE)
  centers = list()
  
  for(j in 1:K) {
    if (gen_norm) {
      centers[[j]] = rnorm(n = D, mean = 0, sd = sd_center)
    } else {
      centers[[j]] = runif(n = D, min = -2 * sd_center, max = 2 * sd_center)
    }
  }
  
  locs = matrix(nrow= Nnodes, ncol = D)
  for(j in 1:Nnodes) {
    locs[j,] = centers[[group_assign[j]]] + rnorm(n = D, mean = 0, sd = 1)
  }
  
  netm = new("NetworkModelLSM", Nnodes = Nnodes, locs = locs, alpha = runif(1,0,10))
  return(netm)
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
  Nnodes = getNnodes(NetM)
  odds = exp(NetM@alpha - dist(NetM@locs))
  prob_mat = matrix(0, nrow = Nnodes, ncol = Nnodes)
  prob_mat[lower.tri(prob_mat)] = odds / (1 + odds)
  prob_mat = prob_mat + t(prob_mat)
  return(prob_mat) 
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelLSM"), getEdgeProbMat.NetworkModelLSM)

