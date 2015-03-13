# This file defines the NetworkModel class.

## TODO: Add validity checks?



# Defines the NetworkModel Class and its subclasses -----------------------
setClass("NetworkModel", representation(Nnodes = "numeric"))

#' Instantiates an object of class NetworkModel
#' 
#' @param Nnodes Number of nodes in network
#' @param type Type of network model (accepts 'block', 'tree', 'random', 'latent') ('none' is also accepted, but the resulting object isn't really interesting; its mainly for testing the class)
#' @param model_param Model parameters specified by set_model_param()
#' 
#' @return Object of class NetworkModel
#' 
#' @export
#' 
NetworkModel = function(Nnodes = 10, type = "none", model_param = set_model_param()) {
  # creates a default networkmodel object -- this is the default constructor, but probably should never be used. the specific ones for a specific model should be used. 
  # maybe want to make this eventually call the network generation methods
  
  ## currently this does nothing but return a lame object, that satisfies the generic methods (although the generic methods would not do much?)
  if (type == "none") { return(new("NetworkModel", Nnodes = Nnodes)) }
  if (type == "block") { return(NetworkModelSBM(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "tree") { return(NetworkModelHRG(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "latent") { return(NetworkModelLSM(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "random") { return(NetworkModelRND(Nnodes = Nnodes, model_param = model_param)) }
  stop("Invalid 'type' specified")
}


# Define Generic Methods --------------------------------------------------


### New method to extract the number of nodes from a network
setGeneric("getNnodes", function(NetM) standardGeneric("getNnodes"))

#' Generic Function -- Extracts the number of nodes in the network model
#' 
#' @param NetM NetworkModel object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes = function(NetM) { NULL }

#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkModel
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkModel = function(NetM) {
  # NetM should be object of type NetworkModel
  return(NetM@Nnodes)
}
setMethod("getNnodes", signature("NetworkModel"), getNnodes.NetworkModel)




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
sampleNetwork = function(NetM, Nobs = 1, Nsim = 1) { return(NULL) }




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
setMethod("sampleNetwork", signature = (NetM = "NetworkModel"), sampleNetwork.NetworkModel)




# Define Default print/summary methods ------------------------------------





## Methods for network-model specific information
## These are the functions that need to be defined for each specific network model. 


# Define generic functions that are model-specific ------------------------


setGeneric("getNetType", function(NetM) standardGeneric("getNetType"))

#' Generic Function -- Extracts the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType = function(NetM) { 
  NULL 
}




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






# Define generic null-model values so this doesn't crash ------------------

#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkModel = function(NetM) {
  "none" 
}
setMethod("getNetType", signature(NetM = "NetworkModel"), getNetType.NetworkModel)

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
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModel"), getEdgeProbMat.NetworkModel)


