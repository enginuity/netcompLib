## Network helper functions


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (hide_edges)
#' Hides a random set of edges
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param frac [double] :: Fraction of edges to hide
#' @param template temp
#' @param invert_template temp
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

