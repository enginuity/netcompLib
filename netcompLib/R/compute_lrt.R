##@S Functions related to computation of the likelihood ratio test


# Functions -- computing MLE estimate for fixed structure -----------------

#' Computes the cell-wise MLE estimate
#' 
#' @param tree Tree object (structure)
#' @param obsadjs Input observed network
#' @param exp_child Expanded_children information (for runtime improvements)
#' 
#' @return Cell-wise MLE estimates
#' 
#' @export
#' 
mle_est = function(tree, obsadjs, exp_child = NULL) {
  NN = dim(obsadjs)[3]
  obsadjm = apply(obsadjs, c(1,2), sum)
  
  if (is.null(exp_child)) {
    exp_child = expanded_children_from_tree(tree)
  }
  
  probs = rep(NA, times = length(exp_child))
  counts = rep(NA, times = length(exp_child))
  for(j in 1:length(exp_child)) {
    lcs = exp_child[[j]][[1]]
    rcs = exp_child[[j]][[2]]
    counts[j] = length(lcs) * length(rcs) * NN
    probs[j] = sum(obsadjm[lcs,rcs]) / counts[j]
  }
  return(list(probs = probs, counts = counts))
}


#' Computes MLE probabilities (with convex weighting)
#' 
#' @param tree Tree structure
#' @param obsadjs Observed network
#' @param wt Convex weights
#' 
#' @return Tree-object with probabilities filled in
#' 
#' @export
#' 
cvx_est = function(tree, obsadjs, wt) {
  temp = mle_est(tree,obsadjs)
  N = (length(tree$parents) + 1) / 2
  ests = rep(NA, times = N-1)
  ests[1] = temp$probs[1]
  
  to_try = which(tree$parents == (N+1))
  to_try = to_try[to_try > N]
  while(length(to_try) > 0) {
    j = to_try[1]
    to_try = to_try[-1]
    par = tree$parents[j]
    ests[j-N] = (wt*ests[par-N] + (temp$counts[j-N])*temp$probs[j-N])/ (wt + temp$counts[j-N])
    
    to_add = which(tree$parents == j)
    to_try = c(to_try, to_add[to_add > N])
  }
  tree$prob = ests
  return(tree)
}


#' Compute the likelihood of an observed network (only one observation allowed), given a model
#' 
#' @param tree Tree model to compute likelihood for
#' @param obsadjm Observed adjacency matrix
#' @param epmat Edge probability matrix (tree is ignored IF this is non-NULL)
#' 
#' @return Log-Likelihood
#' 
#' @export
#' 
compute_likelihood = function(tree, obsadjm, epmat = NULL) {
  if (is.null(epmat)) epmat = edge_probs(tree)
  
  epmat[epmat == 0] = NA
  epmat[epmat == 1] = NA
  
  return(sum(obsadjm * log(epmat) + (1 - obsadjm) * log(1 - epmat), na.rm = TRUE)/2)
}


#' Compute likelihoods for a sequence of adjacency matrices
#' 
#' @param tree Tree model to compute likelihood for
#' @param obsadjs Observed adjacency matrices
#' @param epmat Edge probability matrix (tree is ignored IF this is non-NULL)
#' 
#' @return Vector of log-likelihoods
#' 
#' @export
#' 
compute_likelihood_v = function(tree, obsadjs, epmat = NULL) {
  obs = dim(obsadjs)[3]
  logliks = rep(0, times = obs)
  
  for(j in 1:obs) {
    logliks[j] = compute_likelihood(tree, obsadjs[,,j], epmat = epmat)
  }
  return(logliks)
}


#' Computes the number of data points used to estimate each edge group probability
#' 
#' @param tree Fixed tree
#' 
#' @return Vector of edge group sizes
#' 
#' @export
#' 
edge_prob_n = function(tree) {
  n = tree$nodes
  mat = matrix(NA, nrow = n, ncol = n)
  
  ex_c = expanded_children_from_tree(tree)
  for(j in 1:(n-1)) {
    c1 = ex_c[[j]][[1]]
    c2 = ex_c[[j]][[2]]
    NN = length(c1) * length(c2)
    mat[c1,c2] = NN
    mat[c2,c1] = NN
  }
  return(mat)
}


#' Given a fixed tree structure and an observed network, compute the MLE log likelihood 
#' 
#' @param tree fixed tree structure
#' @param obsadjs observed network
#' @param cellwise T/F: Return cell-wise MLE results?
#' @param exp_child expanded children (to make runtime faster)
#' 
#' @return cell-wise log likelihood or tree loglik
#' 
#' @export
#' 
mle_loglik = function(tree, obsadjs, cellwise = FALSE, exp_child = NULL) {
  z = mle_est(tree, obsadjs, exp_child)
  b = which((z$probs > 0) & (z$probs < 1))
  p = z$probs[b]
  c = z$counts[b]
  x = round(c * p, 0)
  if (!cellwise) {
    return(sum( x * log(p) + (c - x) * log(1 - p)))
  } else {
    p = z$probs
    c = z$counts
    x = round(c*p, 0)
    res = x * log(p) + (c - x) * log(1 - p)
    res[    is.nan(res)] = 0
    return(res)
  }
}

