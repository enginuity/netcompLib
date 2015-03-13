# Defines the NetworkModelHRG class
setClass("NetworkModelHRG", representation(parents = "numeric", children = "list", prob = "numeric"), contains = "NetworkModel")

#' Constructor for HRG network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelHRG object
#' 
#' @export
#' 
NetworkModelHRG = function(Nnodes = 10, model_param = set_model_param()) {
  ## TODO: - fill this in eventually 
  
  # helper function that generates a fixed structure tree (as close to binary tree as possible)
  starter_tree = function(Nnodes = 10) {
    # Written by Andrew Thomas
    
    #format: id, left child, right child, value.
    id1 <- Nnodes+1:(Nnodes-1)  #internal node ids.
    lcrc <- array(NA, c(2,Nnodes-1)); lcrc[1:(Nnodes-2)] <- 2:(Nnodes-1)+Nnodes; lcrc[(Nnodes-1):(2*(Nnodes-1))] <- 1:Nnodes
    children <- list(); for (kk in 1:dim(lcrc)[2]) children[[kk]] <- lcrc[,kk]
    
    val <- runif(Nnodes-1)
    
    parents <- rep(0, 2*Nnodes - 1)
    for (kk in 1:length(children)) parents[children[[kk]]] <- kk+Nnodes
    
    output <- list(prob=val,
                   children=children,
                   parents=parents,
                   Nnodes=Nnodes)
    
    return(output)
  }
  
  ## Helper functions
  node_to_add = function(inp_clist, cur_node, NN) {
    lr = 1
    if (runif(1) > 0.5) { lr = 2 }
    
    if (inp_clist[[cur_node - NN]][lr] == -1) {
      return(c(cur_node, lr))
    } else {
      return(node_to_add(inp_clist, inp_clist[[cur_node - NN]][lr], NN))
    }
  }
  
  find_first_empty = function(inp_clist) {
    nores = TRUE
    j = 0
    while(nores) {
      j = j + 1
      z = which(inp_clist[[j]] == -1)
      nores = (length(z) == 0)
    }
    return(c(j,z[1]))
  }
  
  
  ## Note that random_plimit only applies in the case of "random" random_type. (and for "left" random_type also). 
  
  ## TODO: [TEMP] Using these two variables shouldn't be necessary. 
  random_type = model_param$tree_type
  random_plimit = c(model_param$pmin, model_param$pmax)
  
  res_tree = list()
  
  if (random_type == "original") {
    res_tree = starter_tree(Nnodes)  
  } else if (random_type == "random") {
    clist = list()
    for(j in 1:(Nnodes-1)) {
      clist[[j]] = c(-1,-1)
    }
    
    for(j in (Nnodes+2):(2*Nnodes-1)) {
      z = node_to_add(clist, cur_node = (Nnodes+1), Nnodes)
      clist[[z[1] - Nnodes]][z[2]] = j
    }
    
    ## fill remaining children :
    for(j in sample(1:Nnodes, size = Nnodes)) {
      z = find_first_empty(clist)
      clist[[z[1]]][z[2]] = j
    }
    
    pars = rep(0, times = (2*Nnodes - 1))
    for(j in (Nnodes+1):(2*Nnodes-1)) {
      pars[clist[[j-Nnodes]]] = j
    }
    
  } else if (random_type == "left") {
    pars = rep(0, times = 2*Nnodes - 1)
    pars[(Nnodes+2):(2*Nnodes-1)] = (Nnodes+1):(2*Nnodes-2)
    pars[1:Nnodes] = sample(c( (Nnodes+1):(2*Nnodes-1), (2*Nnodes-1)), size = Nnodes)
    
  } else {
    stop(paste("ERROR: no such type", random_type))
  }
  
  res_tree$parents = pars
  res_tree$children = tree_from_parents(pars)
  res_tree$prob = runif(Nnodes-1, min = random_plimit[1], max = random_plimit[2])
  
  netm = new("NetworkModelHRG", Nnodes = Nnodes, parents = res_tree$parents, 
             children = res_tree$children, prob = res_tree$prob)
  return(netm)
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelHRG objects, this returns "tree"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'tree'
#' 
#' @export
#' 
getNetType.NetworkModelHRG = function(NetM) { "tree" }
setMethod("getNetType", signature(NetM = "NetworkModelHRG"), getNetType.NetworkModelHRG)


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelHRG = function(NetM) {
  nn = getNnodes(NetM)
  ## TODO: Rewrite closest_ancestor, so this section of code isn't this ugly. 
  
  tm = list(prob = NetM@prob, children = NetM@children, parents = NetM@parents, nodes = nn)
  clo.anc = closest_ancestor(tm)$anc.table
  
  out <- matrix(0, nn, nn)
  series <- lower_diag(nn)
  out[series] <- tm$prob[clo.anc[series]-nn]
  out = out + t(out)
  return(out)
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelHRG"), getEdgeProbMat.NetworkModelHRG)

