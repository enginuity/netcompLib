##@S Functions to compute various measures of distances between two probability edge matrices (or sometimes trees)


#' Compute scaled absolute distance between two trees
#' Distance matrix: abs|p1-p2|/se(p1-p2)
#' 
#' @param tree1 Tree structure 1
#' @param tree2 Tree structure 2
#' 
#' @return Matrix of distances
#' 
#' @export
#' 
m_tree_dist_absscale = function(tree1, tree2) {
  p1 = edge_probs(tree1)
  p2 = edge_probs(tree2)
  
  m_mat = edge_prob_n(tree1)
  n_mat = edge_prob_n(tree2)
  res_mat = abs(p1 - p2) / sqrt(p1*(1-p1) / m_mat + p2 * (1 - p2)/n_mat)
  diag(res_mat) = 0
  return(res_mat)
}


#' Computes distance between two trees using squared scaled distances
#' Distance matrix: (p1-p2)^2/se(p1-p2)
#' 
#' @param tree1 Tree structure 1
#' @param tree2 Tree structure 2
#' 
#' @return Distance matrix
#' 
#' @export
#' 
m_tree_dist_sqscale = function(tree1, tree2) {
  p1 = edge_probs(tree1)
  p2 = edge_probs(tree2)
  
  m_mat = edge_prob_n(tree1)
  n_mat = edge_prob_n(tree2)
  res_mat = (p1 - p2)^2 / sqrt(p1*(1-p1) / m_mat + p2 * (1 - p2)/n_mat)
  diag(res_mat) = 0
  return(res_mat)
}


#' Computes absolute distance between edge probability matrices
#' Distance matrix: |p1-p2|
#' 
#' @param tree1 Tree structure 1
#' @param tree2 Tree structure 2
#' 
#' @return Absolute edge probability distance matrix
#' 
#' @export
#' 
m_tree_dist_abs = function (tree1, tree2) {
  p1 = edge_probs(tree1)
  p2 = edge_probs(tree2)
  
  return(abs(p1-p2))
}


#' Computes distance between two trees
#' 
#' @param tree1 Tree structure 1
#' @param tree2 Tree structure 2
#' 
#' @return Value of distance measure
#' 
#' @export
#' 
tree_distance = function(tree1, tree2) {
  counts = mle_est(tree1, cmn_network(tree1))$count
  
  d = abs(tree1$prob - tree2$prob)
  p = (tree1$prob + tree2$prob) / 2
  se = sqrt( p * (1-p) * (2 / counts))
  return(sum( (d/se)[se > 0]))
}


#' Compute the KL distance between two trees
#' 
#' @param tree1 Tree structure 1
#' @param tree2 Tree structure 2
#' 
#' @return KL distance of the two trees
#' 
#' @export
#' 
compute_KLdist = function(tree1, tree2) {
  p1 = edge_probs(tree1)
  p2 = edge_probs(tree2)
  diag(p1) <- 0.5
  diag(p2) <- 0.5
  
  q1 = 1 - p1
  q2 = 1 - p2
  
  res = log(q1) - log(q2) + p1 * (log( (p1 / q1) / (p2 / q2 ) ) )
  diag(res) <- 0
  return(sum(res))
}



compute_KLdist_epmat = function(p1 = NULL, p2 = NULL, plist = NULL) {
  ## computes the KL distance if p1, p2 are edge prob matrices (or if plist is a list of two edge probability matrices)
  if (is.null(p1)) {
    p1 = plist[[1]]
    p2 = plist[[2]]
  }
  diag(p1) <- 0.5
  diag(p2) <- 0.5
  
  q1 = 1 - p1
  q2 = 1 - p2
  
  res = log(q1) - log(q2) + p1 * (log( (p1 / q1) / (p2 / q2 ) ) )
  diag(res) <- 0
  return(sum(res))
}


compute_sym_KLd = function(p1, p2) {
  ## computes the symmetric kl distance
  
  k1 = compute_KLdist_epmat(p1, p2)
  k2 = compute_KLdist_epmat(p2, p1)
  return(0.5 * (k1 + k2))
}



