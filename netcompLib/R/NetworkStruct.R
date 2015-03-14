# This file defines the NetworkStruct class.

## TODO: Add validity checks?



# Defines the NetworkStruct Class and its subclasses -----------------------
setClass("NetworkStruct", representation(Nnodes = "numeric"))

#' Instantiates an object of class NetworkStruct
#' 
#' @param Nnodes Number of nodes in network
#' @param type Type of network model (accepts 'block', 'tree', 'random', 'latent') ('none' is also accepted, but the resulting object isn't really interesting; its mainly for testing the class)
#' @param model_param Model parameters specified by set_model_param()
#' 
#' @return Object of class NetworkStruct
#' 
#' @export
#' 
NetworkStruct = function(Nnodes = 10, type = "none", model_param = set_model_param()) {
  # creates a default NetworkStruct object -- this is the default constructor, but probably should never be used. the specific ones for a specific model should be used. 
  # maybe want to make this eventually call the network generation methods
  
  ## currently this does nothing but return a lame object, that satisfies the generic methods (although the generic methods would not do much?)
  if (type == "none") { return(new("NetworkStruct", Nnodes = Nnodes)) }
  if (type == "block") { return(NetworkStructSBM(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "tree") { return(NetworkStructHRG(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "latent") { stop("Not valid with latent space models") }
  if (type == "random") { return(NetworkStructRND(Nnodes = Nnodes, model_param = model_param)) }
  stop("Invalid 'type' specified")
}


# Define Generic Methods --------------------------------------------------

#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkStruct
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkStruct = function(NetM) {
  # NetM should be object of type NetworkStruct
  return(NetM@Nnodes)
}
setMethod("getNnodes", signature("NetworkStruct"), getNnodes.NetworkStruct)




# Define Default print/summary methods ------------------------------------





## Methods for network-model specific information
## These are the functions that need to be defined for each specific network model. 

#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkStruct = function(NetM) {
  "none" 
}
setMethod("getNetType", signature(NetM = "NetworkStruct"), getNetType.NetworkStruct)





# Define new generics -----------------------------------------------------

setGeneric("computePval", function(NetS, adja1, adja2, Nobs, pl) standardGeneric("computePval"))
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval)
#' <What does this function do>
#' 
#' @param NetS temp
#' @param adja1 temp
#' @param adja2 temp
#' @param Nobs temp
#' @param pl temp
#' 
#' @return temp
#' 
#' @export
#' 
computePval = function(NetS, adja1, adja2, Nobs, pl) {
  stop("Placeholder for documentation purposes")
}

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStruct)
#' <What does this function do>
#' 
#' @param NetS temp
#' @param adja1 temp
#' @param adja2 temp
#' @param Nobs temp
#' @param pl temp
#' 
#' @return temp
#' 
#' @export
#' 
computePval.NetworkStruct = function(NetS, adja1, adja2, Nobs, pl) {
  stop("Not implemented for this case")
}
setMethod("computePval", signature(NetS = "NetworkStruct"), computePval.NetworkStruct)
