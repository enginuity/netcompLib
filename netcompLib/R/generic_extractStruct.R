##@S Function to create a NetworkStruct object from a NetworkModel object

setGeneric("extractStruct", function(NetM) standardGeneric("extractStruct"))

#' Extracts the Edge Partition from a Network Model
#' 
#' @param NetM NetworkModel object
#' 
#' @return NetworkStruct object that corresponds with the input model
#' 
#' @export
#' 
extractStruct = function(NetM) {
  return(NULL) 
}


extractStruct.NetworkModel = function(NetM) {
  stop("This usage case is not defined. ")
}


extractStruct.NetworkModelPair = function(NetM) {
  netsl = new("NetworkStructList", Nnodes = getNnodes(NetM), 
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:32 2015--
              models = list(extractStruct(NetM@m1), extractStruct(NetM@m2)))
  return(netsl)
}


extractStruct.NetworkModelSBM = function(NetM) {
  # group assignments
  ga = NetM@assign
  NClass = length(unique(ga))
  
  counts = rep(0, times = NClass + NClass * (NClass - 1) / 2)
  correction = rep(0, times = NClass + NClass * (NClass - 1) / 2)
  
  expanded = list()
  cur = 1
  for(k in 1:NClass) { for(m in 1:NClass) {
    if (k >= m) {
      if (k == m) { 
        counts[cur] = sum(ga == k) * (sum(ga == k) - 1) / 2 
        correction[cur] = 2
      } else {
        counts[cur] = sum(ga == k) * sum(ga == m)
        correction[cur] = 1
      }
      expanded[[cur]] = list(which(ga == k), which(ga == m))
      cur = cur + 1
    }
  }}
  
  nets = new("NetworkStructSBM", Nnodes = getNnodes(NetM), groups = ga, counts = counts, expand = expanded, correct = correction)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:33 2015--
  return(nets)
}


extractStruct.NetworkModelHRG = function(NetM) {
  tr = list(prob = NetM@prob, children = NetM@children, parents = NetM@parents, nodes = getNnodes(NetM))
  expc = expanded_children_from_tree(tr)
  
  nets = new("NetworkStructHRG", Nnodes = getNnodes(NetM), tree_list = tr, expand = expc, 
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:33 2015--
             counts = sapply(expc, function(x) {length(x[[1]]) * length(x[[2]]) }))
  return(nets)
}


extractStruct.NetworkModelRND = function(NetM) {
  nets = new("NetworkStructRND", Nnodes = getNnodes(NetM), counts = NetM@counts, ids = NetM@ids)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:32 2015--
  return(nets)
}


# setMethod ---------------------------------------------------------------
setMethod("extractStruct", signature(NetM = "NetworkModel"), extractStruct.NetworkModel)
setMethod("extractStruct", signature(NetM = "NetworkModelHRG"), extractStruct.NetworkModelHRG)
setMethod("extractStruct", signature(NetM = "NetworkModelPair"), extractStruct.NetworkModelPair)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:28 2015--
setMethod("extractStruct", signature(NetM = "NetworkModelSBM"), extractStruct.NetworkModelSBM)
setMethod("extractStruct", signature(NetM = "NetworkModelRND"), extractStruct.NetworkModelRND)
