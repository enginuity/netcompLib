##@S Starts off as Andrew's codefile. Not sure if changes have been made to this file... presumably yes. 



#' Quick output of tree
#' 
#' @param tree Tree to display
#' 
#' @return Display of the tree (to be 'print'ed)
#' 
#' @export
#' 
display_tree = function(tree) {
  # Written by Andrew Thomas
  
  max.length <- max(sapply(tree$children, length))
  object <- sapply(tree$children, function (ll) c(ll, rep(0, max.length-length(ll))),simplify="array")
  cbind(tree$nodes+1:length(tree$children), t(object), tree$prob)
}


#' Generate edge matrix from a simple block model
#' 
#' @param nodes Number of nodes
#' @param blocks Number of blocks
#' @param block.matrix Block probability matrix (randomly generated if not passed in)
#' @param replicates Number of edge matrices to generate
#' 
#' @return Array of generated edge matrices
#' 
#' @export
#' 
simple_block_model = function(nodes=100, blocks=3, block.matrix = {out <- array(0.05, rep(blocks,2)); diag(out) <- 0.3; out}, replicates=10) {
  # Written by Andrew Thomas
  
  membership <- sample(1:blocks, nodes, replace=TRUE)
  piece <- membership[rep(1:nodes, nodes)]+blocks*(membership[sort(rep(1:nodes, nodes))]-1)
  hold <- array (rbinom(nodes*nodes*replicates, 1, block.matrix[piece]), c(nodes, nodes, replicates))
  for (kk in 1:replicates) diag(hold[,,kk]) <- 0
  return(hold)
}




#' Extract edge probabilities (in matrix format)
#' 
#' @param tree CMN tree object
#' 
#' @return Edge probability matrix
#' 
#' @export
#' 
edge_probs = function(tree) {
  # Written by Andrew Thomas
  
  #tree=trial$choice.tree.ret
  
  nn <- tree$nodes  #n-1 internal nodes.
  
  clo.anc <- closest_ancestor(tree)$anc.table
  out <- array(0, c(nn, nn))
  series <- lower_diag(nn)
  
  out[series] <- tree$prob[clo.anc[series]-nn]
  out <- out+t(out)
  return(out)
}


#' Computes distance between two trees
#'   
#' Formula is the sum of cell-wise absolute differences for edge probability matrices
#' 
#' @param tree1 CMN tree 1
#' @param tree2 CMN tree 2
#' 
#' @return Number: distance
#' 
#' @export
#' 
tree_distance = function(tree1, tree2) {
  # Written by Andrew Thomas
  
  sum(abs(edge_probs(tree1)-edge_probs(tree2)))
}


#' Computes the closest ancestor for each pair of nodes
#' 
#' @param parents.vector Vector of parent nodes
#' 
#' @return Matrix giving closest ancestor for each pair of nodes
#' 
#' @export
#' 
closest_ancestor_from_parents = function(parents.vector) {
  # Written by Andrew Thomas
  
  # super slow
  #parents.vector=c(6,5,7,7,0,5,6)
  nn <- (length(parents.vector)+1)/2
  history <- function(kk) {
    out <- kk
    while (parents.vector[kk] > 0) {kk <- parents.vector[kk]; out <- c(kk,out)}
    return(out)
  }
  history.all <- lapply(1:length(parents.vector), history)
  length.set <- sapply(history.all, length)
  
  out.table <- array(0, c(nn,nn))
  trial <- sapply(1:(nn-1),function(ii)
                  sapply ((ii+1):nn, function(jj) {
    kk <- 1;
    while (kk < length.set[ii] & kk < length.set[jj] & history.all[[ii]][kk]==history.all[[jj]][kk]) kk <- kk+1
    pick <- history.all[[ii]][kk-1]
    out.table[ii,jj] <- out.table[jj,ii] <- pick
  }))
  diag(out.table) <- 1:nn

  return(out.table)   
}


#' Compute parent vector when it's a generated binary tree (from starter_tree)
#' 
#' @param tree CMN tree object
#' 
#' @return Parent vector
#' 
#' @export
#' 
parents_from_binary_tree = function(tree) {
  # Written by Andrew Thomas
  ## TODO: Is this function even necessary?
  
  total.nn <- dim(tree)[1]*2+1
  out <- rep(0, total.nn)
  out[tree[,2]] <- out[tree[,3]] <- tree[,1]
  return(out)
}







## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (plot_cmn)
#' <<BasicInfo>> 
#' 
#' @param tree.obj temp
#' @param ... temp
#' 
#' @return temp
#' 
#' @export
#' 
plot_cmn = function(tree.obj, ...) {
  # Written by Andrew Thomas
  
  tree.obj <- matrix(tree.obj, ncol=4)
  edgelist <- rbind(tree.obj[,1:2], tree.obj[,c(1,3)]) 
  plot(electrograph(edgelist),...)
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (trial_steps)
#' <<BasicInfo>> 
#' 
#' @param  temp
#' 
#' @return temp
#' 
#' @export
#' 
trial_steps = function() {
  # Written by Andrew Thomas
  
  #source("allcode.R")
  #library(ElectroGraph)
  blocks <- 3
  netdraw <- simple_block_model(block.matrix={out <- array(0.05, rep(blocks,2)); diag(out) <- 0.7; out})
  x11(); first.try <- plot(electrograph(netdraw[,,1]), node.colors=5)
  x11(); new.try <- plot(electrograph(netdraw[,,1]), manual.coords=first.try$coord, node.colors=2+1*((1:100)>33)+1*((1:100)>66))


  trial <- cmn_mcmc(netdraw, iterations=100) #,tree=tree.start)
#|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
 
}

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (single_procedure_test)
#' <<BasicInfo>> 
#' 
#' @param tree.start temp
#' @param id temp
#' @param iter temp
#' 
#' @return temp
#' 
#' @export
#' 
single_procedure_test = function(tree.start, id=1, iter=30000) {
  # Written by Andrew Thomas
  
  #tree.start=trial.block.tog$tree[,100]; id=1; iter=2000

  
  gen.net <- cmn_network (tree.start)
  trial <- cmn_mcmc(gen.net, iterations=iter) #,tree=tree.start)
#|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
  save(gen.net, trial, file=paste(id,"-cmn.RData",sep=""))
  
  pick <- max(which(trial$loglik==max(trial$loglik)))
  out <- list(gen.net=gen.net, tree=trial$tree[,pick])
  return(out)
}



