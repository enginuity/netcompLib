##@S Function to create a NetworkModel object from a NetworkStruct object

setGeneric("extractModel", function(NetS, model_param) standardGeneric("extractModel"))


#' Generates a random network model from a fixed structure
#' 
#' @param NetS temp
#' @param model_param temp
#' 
#' @return temp
#' 
#' @export
#' 
extractModel = function(NetS, model_param = set_model_param()) {
  stop("No implementation for template")
}


extractModel.NetworkStruct = function(NetS, model_param = set_model_param()) { 
  stop("No implmentation for template NetworkStruct")
}


extractModel.NetworkStructList = function(NetS, model_param = set_model_param()) {
  return(lapply(NetS@models, extractModel))
}


extractModel.NetworkStructSBM = function(NetS, model_param = set_model_param()) {
  model_param$Nnodes = getNnodes(NetS)
  model_param$type = "block"
  model_param$block_assign = NetS@assign
  
  res = NetworkModel(model_param)
  
  return(res)
}


extractModel.NetworkStructRND = function(NetS, model_param = set_model_param()) {
  model_param$Nnodes = getNnodes(NetS)
  model_param$type = "random"
  model_param$random_ngroups = length(NetS@ids)
  
  res = NetworkModel(model_param)
  res@ids = NetS@counts
  res@counts = NetS@counts
  
  return(res)
}


extractModel.NetworkStructHRG = function(NetS, model_param = set_model_param()) {
  model_param$Nnodes = getNnodes(NetS)
  model_param$type = "tree"
  
  res = NetworkModel(model_param)
  res@parents = NetS@tree_list$parents
  res@children = NetS@tree_list$children
  
  return(res)
}


# setMethod ---------------------------------------------------------------
setMethod("extractModel", signature(NetS = "NetworkStruct"), extractModel.NetworkStruct)
setMethod("extractModel", signature(NetS = "NetworkStructList"), extractModel.NetworkStructList)
setMethod("extractModel", signature(NetS = "NetworkStructSBM"), extractModel.NetworkStructSBM)
setMethod("extractModel", signature(NetS = "NetworkStructRND"), extractModel.NetworkStructRND)
setMethod("extractModel", signature(NetS = "NetworkStructHRG"), extractModel.NetworkStructHRG)
