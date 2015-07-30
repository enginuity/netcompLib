
#' Compute depths of each node
#' 
#' @param parents [vector-int] :: Valued from 0 to length; these indicate which index the parent node is at
#' 
#' @return [vector-int] :: Depth of each node
#' 
#' @export
#' 
depth_from_parents = function(parents) {
  nn = (length(parents)+1)/2
  
  depths = rep(NA, times = nn)
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
#' @param parents [vector-int] :: Indices of parents
#' 
#' @return [list] :: A tree in list format (left child index, right child index)
#' 
#' @export
#' 
tree_from_parents = function(parents) {
  #Assume internal nodes are greater than leaves.
  # parents=c(6,6,7,7,0,5,5)
  total_nn = length(parents)
  nn = (total_nn+1)/2
  if (any(parents %in% 1:nn)) stop("Improper tree form.")
  children = list()
  for (kk in 1:(nn-1)) children[[kk]] = numeric(0)
  for (kk in 1:total_nn) {
    if (parents[kk]>nn) children[[parents[kk]-nn]] = c(children[[parents[kk]-nn]], kk)
  }
  return(children)
}


#' Extract expanded children for a tree
#' 
#' @param tree [list] :: CMN tree object
#' 
#' @return [list] :: A tree in list format (left child index, right child index)
#' 
#' @export
#' 
expanded_children_from_tree = function(tree) {
  ## TODO: Modify this for different tree format?
  
  expanded.children <- lapply(tree$children, as.list)
  #put them in order.
  depths <- depth_from_parents(tree$parents)[-(1:tree$nodes)]
  
  #from the bottom up, replace entries greater than nn with their (unlisted) children.
  for (pick in rev(order(depths))) {
    
    for (kk in 1:length(expanded.children[[pick]])) { # lists.
      keeper <- expanded.children[[pick]][[kk]][expanded.children[[pick]][[kk]] <= tree$nodes]
      repl <- expanded.children[[pick]][[kk]][which(expanded.children[[pick]][[kk]] > tree$nodes)]
      if (length(repl)>0) for (jj in 1:length(repl))
        keeper <- c(keeper, unlist(expanded.children[[repl[jj]-tree$nodes]]))
      expanded.children[[pick]][[kk]] <- keeper
    }
  }

  return(expanded.children)  
}


#' Compute ancestor table
#' 
#' @param expanded.children List of expanded children
#' @param nn Number nodes in tree
#' @param nodes.to.update Only update a subset of nodes 
#' 
#' @return Ancestor table (matrix)
#' 
#' @export
#' 
anc_table_from_expanded_children = function(expanded.children, nn, nodes.to.update=NULL) {
  # Written by Andrew Thomas
  
  
  anc.table <- array(NA, rep(nn, 2))
  diag(anc.table) <- 1:nn
  int.nodes <- length(expanded.children)
  
  if (is.null(nodes.to.update)) nodes.to.update <- 1:int.nodes
  
  #build the p matrix.
  for (kk in 1:int.nodes) {
    for (ii in 1:(length(expanded.children[[kk]])-1))
      for (jj in (ii+1):length(expanded.children[[kk]])) {
        anc.table[expanded.children[[kk]][[ii]],
                  expanded.children[[kk]][[jj]]] <-
          anc.table[expanded.children[[kk]][[jj]],
                    expanded.children[[kk]][[ii]]] <- kk+nn
      }
  }
  
  return(anc.table)
}


#' Compute ancestor_table from a tree
#' 
#' @param tree [list] :: CMN tree model object
#' 
#' @return [list] :: expanded_children & ancestor_table
#' 
#' @export
#' 
closest_ancestor = function(tree) {
  nn = tree$nodes
  int_nodes = length(tree$prob)
  tot_nodes  nn + int_nodes
  
  expanded_children = expanded_children_from_tree(tree)
  anc_table = anc_table_from_expanded_children(expanded_children, nn)
  
  return(list(expanded_children=expanded_children, anc_table=anc_table))
}
