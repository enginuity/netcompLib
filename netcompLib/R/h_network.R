## Network helper functions


#' Hides a random set of edges
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param frac [double] :: Proportion of dyads to hide, rounded down to nearest number of dyads. 
#' @param template [matrix] :: Matrix with NA's -- If this is provided, then the NA's in this matrix are copied into 'adjm'
#' @param invert_template [logical] :: If TRUE, the input template is inverted. 
#' 
#' @return [matrix-int] :: Adjacency matrix with hidden edges
#' 
#' @export
#' 
hide_edges = function(adjm, frac = 0.1, template = NULL, invert_template = FALSE) {
  
  ## TODO: Make this work for arrays
  ## copied from netcompSBM
  
  if (is.null(template)) {
    tm = matrix(1:(nrow(adjm)^2), nrow = nrow(adjm))
    vals = tm[upper.tri(tm)]
    vals = sample(vals, size = floor(length(vals) * frac))
    adjm[vals] = NA
    adjm = symmetrize_mat(adjm)
  } else {
    if (invert_template) {
      adjm[!is.na(template)] = NA
    } else {
      adjm[is.na(template)] = NA
    }
    diag(adjm) = 0
  }
  return(adjm)
}


#' Symmetrizes matrix by using upper-triangular portion
#' 
#' @param mat [matrix] :: Input matrix
#' 
#' @return [matrix] :: Filled matrix (with lower triangular portion replaced by the upper triangular portion). 
#' 
#' @export
#' 
symmetrize_mat = function(mat) {
  ## copied from netcompSBM
  diag_mat = diag(mat)
  mat[lower.tri(mat, diag = TRUE)] = 0
  mat = mat + t(mat)
  diag(mat) = diag_mat
  return(mat)
}


graph_laplacian = function(adjm) {
  ## Given an adjacency matrix, this computes the graph laplacian matrix. 
  ## Note, this also works with missing data, by estimating the degree of all nodes by scaling up relative to the amount of missing dyads. 
  
  ## TODO: Extend this to work for adjacency arrays too? 
  
  ## Zero out the diagonal
  diag(adjm) = 0
  N = length(diag(adjm))
  
  ## Estimate/compute the degree: 
  degs = rowSums(adjm, na.rm = TRUE)
  na_prop = rowSums(is.na(adjm)) / (N-1)
  
  ## Get negative of adjacency matrix, and set diagonal appropriately
  adjm = -adjm
  diag(adjm) = degs/(1-na_prop)
  
  return(adjm)
}

