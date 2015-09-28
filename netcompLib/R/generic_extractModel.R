##@S Function to create a NetworkModel object from a NetworkStruct object

setGeneric("extractModel", function(NetS, model_param) standardGeneric("extractModel"))


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
#   
#   
#   ## mode -- default
#   ## if mode = densitydiff, then, assume global density difference parameter
#   ## TODO: Re-implement this?!?! the modes may not be the best way to deal with this... 
#   
#   res = NetworkModel(set_model_param(Nnodes = getNnodes(NetS), type = "block", block_assign = NetS@groups))
#   ## TODO: Implement this as: [issue 1]
#   ## res = extractModel(NetS)
#   
#   Nobs = dim(adja)[3]
#   if (length(dim(adja)) == 2) { Nobs = 1}
#   
#   NG = max(res@assign)
#   gr = res@assign
#   res@probmat = matrix(NA, NG, NG)
#   if (mode == "default") {
#     if (Nobs == 1) { adjm = adja[,,1] }
#     if (Nobs > 1) { adjm = apply(adja, c(1,2), sum) }
#     for(i in 1:NG) { for(j in 1:NG) {
#       ii = which(gr == i); ij = which(gr == j)
#       ci = sum(gr == i); cj = sum(gr == j)
#       if (i != j) { res@probmat[i,j] = sum(adjm[ii,ij]) / (Nobs * ci * cj) }
#       if (i == j) { res@probmat[i,j] = sum(adjm[ii,ij]) / (Nobs * ci * (ci - 1)) }
#     }}
#     return(res) 
#   }
#   
}


extractModel.NetworkStructRND = function(NetS, model_param = set_model_param()) {
  stop("Not implemented yet")
  ## TODO: Fill in code for conversion
}


extractModel.NetworkStructHRG = function(NetS, model_param = set_model_param()) {
  stop("Not implemented")
  ## TODO: Fill in code for conversion
}


# setMethod ---------------------------------------------------------------
setMethod("extractModel", signature(NetS = "NetworkStruct"), extractModel.NetworkStruct)
setMethod("extractModel", signature(NetS = "NetworkStructList"), extractModel.NetworkStructList)
setMethod("extractModel", signature(NetS = "NetworkStructSBM"), extractModel.NetworkStructSBM)
setMethod("extractModel", signature(NetS = "NetworkStructRND"), extractModel.NetworkStructRND)
setMethod("extractModel", signature(NetS = "NetworkStructHRG"), extractModel.NetworkStructHRG)
