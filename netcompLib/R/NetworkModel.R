# This file defines the NetworkModel class.

## TODO: Add validity checks?



# Defines the NetworkModel Class and its subclasses -----------------------
setClass("NetworkModel", representation(Nnodes = "numeric"))
setClass("NetworkModelLSM", contains = "NetworkModel")
setClass("NetworkModelHRG", contains = "NetworkModel")
setClass("NetworkModelRND", contains = "NetworkModel")

# 
# setMethod("getNetType", signature(NetM = "NetworkModelLSM"), function(NetM) { "latent" })
# getNetType(new("NetworkModelLSM"))


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (NetworkModel)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (NetworkModel)
#' Instantiates an object of class NetworkModel
#' 
#' @param model_params temp
#' @param Nnodes temp
#' @param type temp
#' @param model_param temp
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

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sampleNetwork)
#' Generic -- Samples an observed networks from a given network model(s)
#' 
#' @param NetM temp
#' @param Nobs temp
#' 
#' @return temp
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
#' 
#' @return temp
#' 
#' @export
#' 
sampleNetwork.NetworkModel = function(NetM, Nobs = 1, Nsim = 1) { 
  #Nsim -- if you want multiple sets of networks from the same model, generate it all now; will be given in a list if Nsim > 1. 
  
  # uses variables out of scope (epmat, Nnodes)
  gen_one_network = function() {
    res = matrix(0, nrow = Nnodes, Ncol = Nnodes) 
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
getNetType = function(NetM) { NULL }




setGeneric("getEdgeProbMat", function(NetM) standardGeneric("getEdgeProbMat"))

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getEdgeProbMat)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getEdgeProbMat = function(NetM) { return(NULL) }






# Define generic null-model values so this doesn't crash ------------------

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkModel)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkModel = function(NetM) { "none" }
setMethod("getNetType", signature(NetM = "NetworkModel"), getNetType.NetworkModel)

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getEdgeProbMat.NetworkModel)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getEdgeProbMat.NetworkModel = function(NetM) { NULL }
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModel"), getEdgeProbMat.NetworkModel)


