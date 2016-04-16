##@S Function to sample a network from a specific model

setGeneric("sampleNetwork", function(NetM, Nobs = 1, Nsim = 1) standardGeneric("sampleNetwork"))


#' Sample a network from a network model
#' 
#' Nobs is the number of network observations to simulate per simulation. This is built in here since callling sampleNetwork many times is slow (since it requires re-computing the edge probability matrix, which is the same each time given the same network model). There are several cases: 
#' If Nsim = 1 -> Result is an array (with third dimension equal to Nobs)
#' If Nsim > 1 -> Result is a list of arrays
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model to sample network from
#' @param Nobs [int; DEFAULT = 1] :: Number of network observations to simulate
#' @param Nsim [int; DEFAULT = 1] :: Number of simulations
#' 
#' @return [array OR list-array OR list-(array OR list-array)] :: List or array of adjacency matrices. If calling on a \code{\link{NetworkModelPair}}, the output is a list of output from calling this function on each of the two specific network models. 
#' 
#' @export
#' 
sampleNetwork = function(NetM, Nobs = 1, Nsim = 1) { 
  return(NULL) 
}


sampleNetwork.NetworkModel = function(NetM, Nobs = 1, Nsim = 1) { 
  
  Nnodes = getNnodes(NetM)
  epmat = getEdgeProbMat(NetM)
  
  gen_one_network = function() { 
    ## NOTE: uses variables out of scope (epmat, Nnodes)
    res = matrix(0, nrow = Nnodes, ncol = Nnodes) 
    for(k in 1:Nnodes) { for(j in 1:Nnodes) {
      if (k < j) {
        v = rbinom(n = 1, size = 1, prob = epmat[k,j])
        res[j,k] = v; res[k,j] = v
      }
    }}
    return(res)
  }
  
  if (Nsim == 1) {
    
    res = array(0, dim = c(Nnodes, Nnodes, Nobs))
    for(d in 1:Nobs) {
      res[,,d] = gen_one_network()
    }
    return(res)
    
  } else {
    
    reslist = list()
    for(f in 1:Nsim) {
      res = array(0, dim = c(Nnodes, Nnodes, Nobs))
      for(d in 1:Nobs) {
        res[,,d] = gen_one_network()
      }
      reslist[[f]] = res
    }
    return(reslist)
    
  }
}


sampleNetwork.NetworkModelPair = function(NetM, Nobs = 1, Nsim = 1) { 
  return(list(sampleNetwork(NetM@m1, Nobs = Nobs, Nsim = Nsim), sampleNetwork(NetM@m2, Nobs = Nobs, Nsim = Nsim)))
}



# setMethod ---------------------------------------------------------------

setMethod("sampleNetwork", signature(NetM = "NetworkModel"), sampleNetwork.NetworkModel)
setMethod("sampleNetwork", signature(NetM = "NetworkModelPair"), sampleNetwork.NetworkModelPair)

