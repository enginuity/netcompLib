
setGeneric("sampleNetwork", function(NetM, Nobs = 1, ...) standardGeneric("sampleNetwork"))

#' Generic -- Samples an observed networks from a given network model(s)
#' 
#' Nobs is the number of network observations to simulate per simulation. This is built in here since callling sampleNetwork many times is bad (since it requires re-computing the edge probability matrix, which is the same each time given the same network model). There are several cases: 
#' If Nsim = 1, Nobs = 1 -> Result is a matrix
#' If Nsim = 1, Nobs > 1 -> Result is an array (with third dimension equal to Nobs)
#' If Nsim > 1 -> Result is a list of matrices or arrays (depends on Nobs)
#' 
#' @param NetM NetworkModel object
#' @param Nobs Number of network observations to simulate
#' @param Nsim Number of simulations
#' 
#' @return List or array of adjacency matrices
#' 
#' @export
#' 
sampleNetwork = function(NetM, Nobs = 1, Nsim = 1) { 
  return(NULL) 
}




#' Samples any number of observed networks from a given network model
#' 
#' Note -- this implementation requires the edge probability matrix. It's not really efficient (since it re-computes the edge probability matrix each time)
#' 
#' @param NetM Network Model object
#' @param Nobs Number of observed networks to sample
#' @param Nsim Number of simulations (number of times to do Nobs samples)
#' 
#' @return List or array of simulated networks
#' 
#' @export
#' 
sampleNetwork.NetworkModel = function(NetM, Nobs = 1, Nsim = 1) { 
  #Nsim -- if you want multiple sets of networks from the same model, generate it all now; will be given in a list if Nsim > 1. 
  
  # TODO: Use this style of code to improve on network simulation
  #   # Written by Andrew Thomas
  #   
  #   nn <- tm$nodes
  #   
  #   clo.anc <- closest_ancestor(tm)$anc.table
  #   out <- array(0, c(nn, nn, Nobs))
  #   series <- lower_diag(nn)
  #   for (kk in 1:Nobs) {
  #     out[series+nn^2*(kk-1)] <- rbinom(nn*(nn-1)/2, 1, tm$prob[clo.anc[series]-nn])
  #     out[,,kk] <- out[,,kk]+t(out[,,kk])
  #   }
  
  # uses variables out of scope (epmat, Nnodes)
  gen_one_network = function() {
    res = matrix(0, nrow = Nnodes, ncol = Nnodes) 
    for(k in 1:Nnodes) { for(j in 1:Nnodes) {
      if (k < j) {
        v = rbinom(n = 1, size = 1, prob = epmat[k,j])
        res[j,k] = v
        res[k,j] = v
      }
    }}
    return(res)
  }
  
  Nnodes = getNnodes(NetM)
  
  epmat = getEdgeProbMat(NetM)
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


#' Samples any number of observed networks from a given network model
#' 
#' Note -- this implementation requires the edge probability matrix. It's not really efficient (since it re-computes the edge probability matrix each time)
#' 
#' @param NetM Network Model object
#' @param Nobs Number of observed networks to sample
#' @param Nsim Number of simulations (number of times to do Nobs samples)
#' 
#' @return List of two lists/arrays of network models. 
#' 
#' @export
#' 
sampleNetwork.NetworkModelPair = function(NetM, Nobs = 1, Nsim = 1) { 
  return(list(sampleNetwork(NetM@m1, Nobs = Nobs, Nsim = Nsim), sampleNetwork(NetM@m2, Nobs = Nobs, Nsim = Nsim)))
}





# setMethod ---------------------------------------------------------------

setMethod("sampleNetwork", signature = (NetM = "NetworkModel"), sampleNetwork.NetworkModel)
setMethod("sampleNetwork", signature = (NetM = "NetworkModelPair"), sampleNetwork.NetworkModelPair)

