##@S Generic Function that computes the edge probability matrix of a given network model

setGeneric("getEdgeProbMat", function(NetM, mode) standardGeneric("getEdgeProbMat"))

#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Network Model
#' @param mode [char] :: Either "prob" or "group": if "prob", this function returns an edge probability matrix. Else, it returns a matrix where each cell gets a group id number
#' 
#' @return [matrix] :: Edge probability matrix (or group index matrix)
#' 
#' @export
#' 
getEdgeProbMat = function(NetM, mode = "prob") {
  return(NULL) 
}


getEdgeProbMat.NetworkModel = function(NetM, mode) {
  NULL 
}


getEdgeProbMat.NetworkModelPair = function(NetM, mode) {
  return(list(getEdgeProbMat(NetM@m1, mode),getEdgeProbMat(NetM@m2, mode)))
}


getEdgeProbMat.NetworkModelSBM = function(NetM, mode) {
  res = matrix(0, nrow = NetM@Nnodes, ncol = NetM@Nnodes) 
  for(j in 1:(NetM@Nnodes-1)) { for(k in (j+1):NetM@Nnodes) { 
    if (mode == "prob") {
      res[j,k] = NetM@probmat[NetM@assign[j], NetM@assign[k]]
    } else if (mode == "group") {
      res[j,k] = max(NetM@assign[j] + NetM@assign[k] * NetM@Nnodes, NetM@assign[k] + NetM@assign[j] * NetM@Nnodes)
    }
    res[k,j] = res[j,k]
  }}
  
  return(res)
}


getEdgeProbMat.NetworkModelLSM = function(NetM, mode) {
  if (mode != "prob") { stop("Invalid mode (for LSM)") }
  Nnodes = getNnodes(NetM)
  odds = exp(NetM@alpha - dist(NetM@locs))
  prob_mat = matrix(0, nrow = Nnodes, ncol = Nnodes)
  prob_mat[lower.tri(prob_mat)] = odds / (1 + odds)
  prob_mat = prob_mat + t(prob_mat)
  return(prob_mat) 
}


getEdgeProbMat.NetworkModelHRG = function(NetM, mode) {
  nn = getNnodes(NetM)
  ## TODO: Rewrite closest_ancestor, so this section of code isn't this ugly. 
  
  tm = list(prob = NetM@prob, children = NetM@children, parents = NetM@parents, nodes = nn)
  clo.anc = closest_ancestor(tm)$anc_table
  
  out <- matrix(0, nn, nn)
  series <- lower_diag(nn)
  if (mode == "prob") {
    out[series] <- tm$prob[clo.anc[series]-nn]
  } else if (mode == "group") {
    out[series] <- clo.anc[series]
  }  
  out = out + t(out)
  return(out)
}


getEdgeProbMat.NetworkModelRND = function(NetM, mode) {
  pmat = matrix(0, nrow = getNnodes(NetM), ncol = getNnodes(NetM))
  for(j in seq_along(NetM@counts)) {
    if (mode == "prob") {
      pmat[NetM@ids[[j]]] = NetM@prob[j]
    } else if (mode == "group") {
      pmat[NetM@ids[[j]]] = j
    }
  }
  return(pmat + t(pmat))
}


# setMethod ---------------------------------------------------------------
setMethod("getEdgeProbMat", signature(NetM = "NetworkModel"), getEdgeProbMat.NetworkModel)
setMethod("getEdgeProbMat", signature(NetM = "NetworkModelHRG"), getEdgeProbMat.NetworkModelHRG)
setMethod("getEdgeProbMat", signature(NetM = "NetworkModelLSM"), getEdgeProbMat.NetworkModelLSM)
setMethod("getEdgeProbMat", signature(NetM = "NetworkModelPair"), getEdgeProbMat.NetworkModelPair)
setMethod("getEdgeProbMat", signature(NetM = "NetworkModelSBM"), getEdgeProbMat.NetworkModelSBM)
setMethod("getEdgeProbMat", signature(NetM = "NetworkModelRND"), getEdgeProbMat.NetworkModelRND)

