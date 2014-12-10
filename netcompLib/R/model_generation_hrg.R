## Code to generate HRG models

## TODO: [Obselete]


#' Generate a tree: fixed structure, node probabilities at random (uniform)
#' 
#' @param Nnodes Number of nodes in network
#' 
#' @return Tree object (list)
#' 
#' @export
#' 
starter_tree = function(Nnodes=10) {
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
                 nodes=Nnodes)
  
  return(output)
}

#' Generate a random tree structure / model
#' 
#' @param Nnodes Number of nodes in the network
#' @param random_type Type of randomization of structure
#' @param random_plimit Limits for uniform distribution for the edge probabilities
#' 
#' @return HRG tree object
#' 
#' @export
#' 
gen_tree = function(Nnodes, random_type = c("original", "random", "left"), random_plimit= c(0,1)) {
  
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
  
  
  ## Note that random_plimit only applies in the case of "random" random_type. (and for "left" random_type also. 
  
  res_tree = list()
  res_tree$nodes = Nnodes
  res_tree$prob = runif(Nnodes-1, min = random_plimit[1], max = random_plimit[2])
  
  if (random_type == "original") {
    return(starter_tree(Nnodes))
    
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
  
  return(res_tree)
}



#' Generates observations from a CMN network model
#' 
#' @param tree CMN network object
#' @param replicates Number of edge matrices to generate
#' 
#' @return Array of generated edge matrices
#' 
#' @export
#' 
cmn_network = function(tree, replicates=1) {
  # Written by Andrew Thomas
  
  
  #model.out <- trial.block$tree[[20000]; replicates=1
  #tree=net.sol.comb$choice.tree.ret
  
  nn <- tree$nodes
  
  clo.anc <- closest_ancestor(tree)$anc.table
  out <- array(0, c(nn, nn, replicates))
  series <- lower_diag(nn)
  for (kk in 1:replicates) {
    out[series+nn^2*(kk-1)] <- rbinom(nn*(nn-1)/2, 1, tree$prob[clo.anc[series]-nn])
    out[,,kk] <- out[,,kk]+t(out[,,kk])
  }
  return(out)
  #plot (electrograph(out[,,1]))
}


