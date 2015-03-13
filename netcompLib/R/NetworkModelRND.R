# Defines the NetworkModelRND class
setClass("NetworkModelRND", representation(counts = "numeric", prob = "numeric", ids = "list"), contains = "NetworkModel")

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
  rnd_Ngroups = model_param$random_ngroups
  
  # Compute index numbers of the lower-diagonal portion of the matrix
  idm = matrix(1:(Nnodes^2), nrow = Nnodes)
  idm[lower.tri(x = idm, diag = TRUE)] = 0
  idm = as.vector(idm)
  good_ids = idm[idm != 0]
  if (length(good_ids) < rnd_Ngroups) { stop("Too many groups for too few edges") }
  
  # Sample ids to assign into groups
  rand_order = sample(good_ids, size = length(good_ids), replace = FALSE)
  sizes = rep(floor(length(good_ids) / rnd_Ngroups), times = rnd_Ngroups)
  m = length(good_ids) %% rnd_Ngroups
  if (m > 0) {
    sizes[1:m] = sizes[1:m] + 1
  }
  
  counted = 0
  id_list = list()
  for(k in 1:rnd_Ngroups) {
    id_list[[k]] = sort(rand_order[seq(from = counted+1, length.out = sizes[k])])
    counted = counted + sizes[k]
  }
  
  netm = new("NetworkModelRND", Nnodes = Nnodes, counts = sizes, ids = id_list, 
             prob = runif(n = rnd_Ngroups, min = model_param$pmin, max = model_param$pmax))
  return(netm)
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
  pmat = matrix(0, nrow = getNnodes(NetM), ncol = getNnodes(NetM))
  for(j in seq_along(NetM@counts)) {
    pmat[NetM@ids[[j]]] = NetM@prob[j]
  }
  return(pmat + t(pmat))
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelRND"), getEdgeProbMat.NetworkModelRND)


#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelRND = function(NetM) {
  nets = new("NetworkStructRND", Nnodes = getNnodes(NetM), counts = NetM@counts, ids = NetM@ids)
  return(nets)
}
setMethod("extractStruct", signature = (NetM = "NetworkModelRND"), extractStruct.NetworkModelRND)

