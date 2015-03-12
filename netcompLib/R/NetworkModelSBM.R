# Defines the SBMNetworkModel class
setClass("NetworkModelSBM", representation(assign = "numeric", probmat = "matrix"), contains = "NetworkModel")

#' Constructor for SBM network model
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
NetworkModelSBM = function(Nnodes = 10, model_param = set_model_param()) {
  
  ## Helper function for adjusting block model probabilities  
  adjust_blockprobs = function(mod, avgden = 0.4, plimit = c(0.05, 0.95)) {
    classct = table(mod$assign)
    params = length(classct) * (length(classct) - 1 ) / 2 + length(classct)
    coordmat = matrix(NA, nrow = params, ncol = 4)
    colnames(coordmat) = c("row", "col", "ID", "count")
    
    ct = 1
    for(j in 1:length(classct)) { for(k in j:length(classct)) {
      coordmat[ct,] = c(j,k,ct, ifelse(j == k, classct[j] * (classct[j]-1)/2, classct[j] * classct[k]))
      ct = ct + 1
    }}
    
    samp_probvec = function(counts, avgden, plimit) {
      notvalid = TRUE
      while(notvalid) {
        test = runif(n = length(counts), min = plimit[1], max = plimit[2])
        test[1] = (sum(counts) * avgden - sum(counts[-1] * test[-1]))/counts[1]
        if (test[1] > plimit[1] & test[1] < plimit[2]) { notvalid = FALSE }
      }
      return(test)
    }
    probvec = samp_probvec(counts = coordmat[,4], avgden = avgden, plimit = plimit)
    
    pmat = matrix(NA, ncol = length(classct), nrow = length(classct))
    ct = 1
    for(j in 1:length(classct)) { for(k in j:length(classct)) {
      pmat[j,k] = probvec[ct]
      if (j != k) { pmat[k,j] = probvec[ct] }
      ct = ct + 1
    }}
    return(pmat)
  }
  
  # if block assignments are not pre-specified: 
  K = model_param$block_nclass
  if (is.null(model_param$block_assign)) {
    group_assign = sample(1:K, size = Nnodes, replace = TRUE)
  } else { # Use prespecified block assignments
    group_assign = model_param$block_assign
  }
  
  # if block probability matrix is not pre-specified: 
  if (is.null(model_param$block_probs)) {
    # If we want to control average density: 
    if (!is.null(model_param$block_avgdensity)) {
      prob_matrix = adjust_blockprobs(mod = res, avgden = model_param$block_avgdensity, plimit = c(model_param$pmin, model_param$pmax))
    } else {
      prob_matrix = matrix(0, nrow = K, ncol = K)
      for(j in 1:K) { for(i in 1:K) {
        if (i <= j) {
          prob_matrix[i,j] = runif(1, model_param$pmin, model_param$pmax)
          if (i != j) {prob_matrix[j,i] = prob_matrix[i,j]}
        }
      }}
    }
  } else { # Use prespecified block probabilities
    prob_matrix = model_param$block_probs
  }
  
  netm = new("NetworkModelSBM", Nnodes = Nnodes, assign = group_assign, probmat = prob_matrix)
  return(netm)
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelSBM objects, this returns "block"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'block'
#' 
#' @export
#' 
getNetType.NetworkModelSBM = function(NetM) { "block" }
setMethod("getNetType", signature(NetM = "NetworkModelSBM"), getNetType.NetworkModelSBM)


#' Computes the edge probability matrix
#' 
#' @param NetM temp
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
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelSBM"), getEdgeProbMat.NetworkModelSBM)

