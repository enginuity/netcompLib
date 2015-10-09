## Network helper functions


#' Hides a random set of edges
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param frac [double] :: Fraction of edges to hide
#' 
#' @return [matrix-int] :: Adjacency matrix with hidden edges
#' 
#' @export
#' 
hide_edges = function(adjm, frac = 0.1) {
  ## TODO: Make this work for arrays
  ## copied from netcompSBM
  tm = matrix(1:(nrow(adjm)^2), nrow = nrow(adjm))
  vals = tm[upper.tri(tm)]
  vals = sample(vals, size = floor(length(vals) * frac))
  adjm[vals] = NA
  return(symmetrize_mat(adjm))
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

