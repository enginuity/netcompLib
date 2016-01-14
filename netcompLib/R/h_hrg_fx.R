
#' Compute depths of each node
#' 
#' @param parents [vector-int] :: Index location of parent node. Value of 0 indicates that this is the root node.
#' 
#' @return [vector-int] :: Depth of each node (root node has depth 1)
#' 
#' @export
#' 
depth_from_parents = function(parents) {
  n = (length(parents)+1)/2
  
  depths = rep(NA, times = n)
  cur_depth = 0
  cur_nodes = 0
  
  while(any(is.na(depths))) {
    cur_depth = cur_depth + 1
    cur_nodes = which(parents %in% cur_nodes)
    depths[cur_nodes] = cur_depth
  }
  
  return(depths)
}



#' Obtains the lower diagonal area of a square matrix
#' 
#' @param nn [int] :: Number of rows in square matrix
#' 
#' @return [vector-int] :: Indices for lower-diagonal area
#' 
#' @export
#' 
lower_diag = function(nn) {
  inter = 0:(nn^2-1)
  return(which(inter %% nn > floor(inter/nn)))
}


#' Generate tree object from parent vector
#' 
#' @param parents [vector-int] :: Index location of parent node. Value of 0 indicates that this is the root node.
#' 
#' @return [list] :: A tree in list format (left child index, right child index)
#' 
#' @export
#' 
tree_from_parents = function(parents) {
  ## Asume internal nodes have indeices greater than the leaf nodes' indices
  ## Eg. parents=c(6,6,7,7,0,5,5)
  
  n = (length(parents) + 1)/2
  if (any(parents %in% seq_len(n))) { stop("Improper tree form.") }
  children = list()
  for (k in seq_len(n-1)) { children[[k]] = numeric(0) }
  for (k in seq_along(parents)) { 
    ## Iterate through all nodes -- if the node's parents are appropriate, add children
    if (parents[k] > n) { children[[parents[k] - n]] = c(children[[parents[k] - n]], k)}
  }
  return(children)
}


#' Extract expanded children for a tree
#' 
#' @param tree [\code{\link{NetworkModelHRG}}] :: Input network model
#' 
#' @return [list] :: A tree in list format (left child index, right child index)
#' 
#' @export
#' 
expanded_children_from_tree = function(tree) {
  res = lapply(tree@children, as.list)
  n = (length(tree@parents)+1)/2
  
  ## Replace children with their unlisted leaf children, starting with deepest nodes. 
  depths = depth_from_parents(tree@parents)
  for (i in order(depths, decreasing = TRUE)) {
    if (i > n) {
      parentIndex = tree@parents[i] - n
      if (parentIndex > 0) {
        childSide = (i == tree@children[[parentIndex]][2]) + 1 ## Left child or Right child? Left child = 1, Right = 2
        removeIndex = which(i == res[[parentIndex]][[childSide]])
        res[[parentIndex]][[childSide]] = c(res[[parentIndex]][[childSide]][-removeIndex], c(res[[i-n]], recursive = TRUE))
      }
    }
  }
  return(res)
}


#' Compute ancestor_table from a tree
#' 
#' @param tree [\code{\link{NetworkModelHRG}}] :: Input tree
#' 
#' @return [list] :: expanded_children & ancestor_table
#' 
#' @export
#' 
closest_ancestor = function(tree) {
  N = tree@Nnodes
  expCdn = expanded_children_from_tree(tree)
  anc_table = matrix(NA, N, N)
  diag(anc_table) = seq_len(N)
  for(k in seq_along(expCdn)) {
    anc_table[expCdn[[k]][[1]], expCdn[[k]][[2]]] = k + N
    anc_table[expCdn[[k]][[2]], expCdn[[k]][[1]]] = k + N
  }
  return(list(expanded_children=expCdn, anc_table=anc_table))  
}

