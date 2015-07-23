
#' Compute depths of each node
#' 
#' @param parents Parents vector 
#' 
#' @return Vector of depths for each node
#' 
#' @export
#' 
depth_from_parents = function(parents) {
  # Written by Andrew Thomas
  
  nn <- (length(parents)+1)/2
  history <- function(kk) {
    out <- kk
    while (parents[kk] > 0) {kk <- parents[kk]; out <- c(kk,out)}
    return(out)
  }
  history.all <- lapply(1:length(parents), history)
  sapply(history.all, length)
}



#' Obtains the lower diagonal area of a square matrix
#' 
#' @param nn Number of rows in square matrix
#' 
#' @return Indices for lower-diagonal area
#' 
#' @export
#' 
lower_diag = function(nn) {
  # Written by Andrew Thomas
  
  inter <- 0:(nn^2-1)
  which(inter %% nn > floor(inter/nn))
}


#' Generate tree object from parent vector
#' 
#' @param parents Parent vector
#' 
#' @return CMN tree object
#' 
#' @export
#' 
tree_from_parents = function(parents) {
  #Assume internal nodes are greater than leaves.
  
  # Written by Andrew Thomas
  
  #parents=c(6,6,7,7,0,5,5)
  total.nn <- length(parents)
  nn <- (total.nn+1)/2
  if (any(parents %in% 1:nn)) stop("Improper tree form.")
  children <- list(); for (kk in 1:(nn-1)) children[[kk]] <- numeric(0)
  for (kk in 1:total.nn) if (parents[kk]>nn) children[[parents[kk]-nn]] <- c(children[[parents[kk]-nn]], kk)
  return(children)
}


#' Extract expanded children for a tree
#' 
#' @param tree CMN tree object
#' 
#' @return List of terminal leaves for each child
#' 
#' @export
#' 
expanded_children_from_tree = function(tree) {
  # Written by Andrew Thomas
  
  
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
#' @param tree CMN tree model object
#' 
#' @return List of expanded_children & ancestor_table
#' 
#' @export
#' 
closest_ancestor = function(tree) {
  # Written by Andrew Thomas
  
  #tree=tree.prop
  #iterate it.
  
  nn <- tree$nodes
  int.nodes <- length(tree$prob)
  tot.nodes <- nn+int.nodes
  
  
  expanded.children <- expanded_children_from_tree (tree)
  
  anc.table <- anc_table_from_expanded_children(expanded.children, nn)
  
  out <- list(expanded.children=expanded.children,
              anc.table=anc.table)
  return(out)
  
}
