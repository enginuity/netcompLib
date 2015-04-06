##@S Function to create a NetworkModel object from a NetworkStruct object

setGeneric("fitModel", function(NetS, adja) standardGeneric("fitModel"))

## TODO: Figure out why roxygen documentation doesn't seem to work ?!?!?! it behaves differently for whatever reason. there's something strange with this file ???


fitModel = function(NetS, adja) {
  # if adja is null then fill in random uniform probabilities. 
  stop("No implementation in template")
  return(NULL)
}


fitModel.NetworkStruct = function(NetS, adja) { 
  stop("No implmentation for template NetworkStruct")
}


fitModel.NetworkStructList = function(NetS, adja) {
  return(lapply(NetS@models, function(x) { fitModel(x, adja) } ))
}


fitModel.NetworkStructSBM = function(NetS, adja) {
  res = NetworkModel(Nnodes = getNnodes(NetS), type = "block", model_param = set_model_param(block_assign = NetS@groups))
  if (is.null(adja)) {
  } else {
    ## TODO: Fill in code for estimating probabilities. 
    stop("Not implemented for non-null adjacency array")
  }
  return(res)
}


fitModel.NetworkStructRND = function(NetS, adja) {
  if (is.null(adja)) {
    stop("Not implemented")
    ## TODO: Fill in code for conversion
  } else {
    ## TODO: Fill in code for estimating probabilities. 
    stop("Not implemented for non-null adjacency array")
  }
}


fitModel.NetworkStructHRG = function(NetS, adja) {
  if (is.null(adja)) {
    stop("Not implemented")
    ## TODO: Fill in code for conversion
  } else {
    ## TODO: Fill in code for estimating probabilities. 
    stop("Not implemented for non-null adjacency array")
  }
}



# setMethod ---------------------------------------------------------------

setMethod("fitModel", signature(NetS = "NetworkStruct"), fitModel.NetworkStruct)
setMethod("fitModel", signature(NetS = "NetworkStructList"), fitModel.NetworkStructList)
setMethod("fitModel", signature(NetS = "NetworkStructSBM"), fitModel.NetworkStructSBM)
setMethod("fitModel", signature(NetS = "NetworkStructRND"), fitModel.NetworkStructRND)
setMethod("fitModel", signature(NetS = "NetworkStructHRG"), fitModel.NetworkStructHRG)

