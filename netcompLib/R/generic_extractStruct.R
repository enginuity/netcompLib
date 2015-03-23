
setGeneric("extractStruct", function(NetM) standardGeneric("extractStruct"))






#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct = function(NetM) {
  return(NULL) 
}



#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModel = function(NetM) {
  stop("This usage case is not defined. ")
}


#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelPair = function(NetM) {
  netsl = new("NetworkStructList", Nnodes = getNnodes(NetM), 
              models = list(extractStruct(NetM@m1), extractStruct(NetM@m2)))
  return(netsl)
}




#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
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
  return(nets)
}





#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelHRG = function(NetM) {
  tr = list(prob = NetM@prob, children = NetM@children, parents = NetM@parents, nodes = getNnodes(NetM))
  expc = expanded_children_from_tree(tr)
  
  nets = new("NetworkStructHRG", Nnodes = getNnodes(NetM), tree_list = tr, expand = expc, 
             counts = sapply(expc, function(x) {length(x[[1]]) * length(x[[2]]) }))
  return(nets)
}



#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelRND = function(NetM) {
  nets = new("NetworkStructRND", Nnodes = getNnodes(NetM), counts = NetM@counts, ids = NetM@ids)
  return(nets)
}
setMethod("extractStruct", signature = (NetM = "NetworkModelRND"), extractStruct.NetworkModelRND)



# setMethod ---------------------------------------------------------------
setMethod("extractStruct", signature = (NetM = "NetworkModel"), extractStruct.NetworkModel)
setMethod("extractStruct", signature = (NetM = "NetworkModelHRG"), extractStruct.NetworkModelHRG)
setMethod("extractStruct", signature = (NetM = "NetworkModelPair"), extractStruct.NetworkModelPair)
setMethod("extractStruct", signature = (NetM = "NetworkModelSBM"), extractStruct.NetworkModelSBM)
