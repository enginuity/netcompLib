setClass("NetworkStructList", representation(models = "list"), contains = "NetworkStruct")


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (NetworkStructList)
#' <What does this function do>
#' 
#' @param Nnodes temp
#' @param Nmodels temp
#' @param type temp
#' @param model_param temp
#' 
#' @return temp
#' 
#' @export
#' 
NetworkStructList = function(Nnodes = 10, Nmodels = 10, type = "none", model_param = set_model_param()) {
  # creates a default NetworkStruct object -- this is the default constructor, but probably should never be used. the specific ones for a specific model should be used. 
  # maybe want to make this eventually call the network generation methods
  
  res = replicate(n = Nmodels, expr = NetworkStruct(Nnodes = Nnodes, type = type, model_param = model_param))
  netsl = new("NetworkStructList", Nnodes = Nnodes, models = res)
  return(netsl)
}


# generalize gettype
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructList)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructList = function(NetM) {
  return(sapply(NetM@models, getNetType))
}
setMethod("getNetType", signature(NetM = "NetworkStructList"), getNetType.NetworkStructList)

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructList)
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
computePval.NetworkStructList = function(NetS, adja1, adja2, Nobs, pl) {
  # note that NetS here is a netstructlist. 
  res = lapply(NetS@models, function(x) { computePval(x, adja1, adja2, Nobs, pl) } )
  return(res)
}
setMethod("computePval", signature(NetS = "NetworkStructList"), computePval.NetworkStructList)


