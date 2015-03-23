
setGeneric("getEdgeProbMat", function(NetM) standardGeneric("getEdgeProbMat"))



#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM Network model object
#' 
#' @return Matrix -- the edge probability matrix
#' 
#' @export
#' 
getEdgeProbMat = function(NetM) {
  return(NULL) 
}


#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM Network Model object
#' 
#' @return NULL - no edge probability matrix for generic network model object
#' 
#' @export
#' 
getEdgeProbMat.NetworkModel = function(NetM) {
  NULL 
}


#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM Network Model object
#' 
#' @return NULL - no edge probability matrix for generic network model object
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelPair = function(NetM) {
  return(list(getEdgeProbMat(NetM@m1),getEdgeProbMat(NetM@m2)))
}


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model Object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelSBM = function(NetM) {
  res = matrix(0, nrow = NetM@Nnodes, ncol = NetM@Nnodes) 
  for(j in 1:(NetM@Nnodes-1)) { for(k in (j+1):NetM@Nnodes) { 
    res[j,k] = NetM@probmat[NetM@assign[j], NetM@assign[k]]
    res[k,j] = res[j,k]
  }}
  
  return(res)
}


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


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelHRG = function(NetM) {
  nn = getNnodes(NetM)
  ## TODO: Rewrite closest_ancestor, so this section of code isn't this ugly. 
  
  tm = list(prob = NetM@prob, children = NetM@children, parents = NetM@parents, nodes = nn)
  clo.anc = closest_ancestor(tm)$anc.table
  
  out <- matrix(0, nn, nn)
  series <- lower_diag(nn)
  out[series] <- tm$prob[clo.anc[series]-nn]
  out = out + t(out)
  return(out)
}




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






# setMethod ---------------------------------------------------------------

setMethod("getEdgeProbMat", signature = (NetM = "NetworkModel"), getEdgeProbMat.NetworkModel)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelHRG"), getEdgeProbMat.NetworkModelHRG)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelLSM"), getEdgeProbMat.NetworkModelLSM)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelPair"), getEdgeProbMat.NetworkModelPair)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelSBM"), getEdgeProbMat.NetworkModelSBM)
