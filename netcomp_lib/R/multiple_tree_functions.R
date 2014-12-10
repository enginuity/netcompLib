

#' Computes probability matrix from averaging trees
#' 
#' @param params Weights for averaging trees
#' @param list_trees List of trees to average
#' 
#' @return Matrix of average probabilities
#' 
#' @export
#' 
tree_av = function(params, list_trees) {
  ## returns probability matrix, for averaging the list of trees, using params.
  params = params / sum(params)
  res_mat = 0 * edge_probs(list_trees[[1]])
  for(j in 1:length(params)) {
    res_mat = res_mat + params[j] * edge_probs(list_trees[[j]])
  }
  
  return(res_mat)
}


#' Compute combined likelihood under average tree model
#' 
#' @param params Parameters for weighted average
#' @param list_trees List of fixed trees
#' @param adj Observed network
#' 
#' @return Likelihood
#' 
#' @export
#' 
comp_lik_comb = function(params, list_trees, adj) {
  if (is.null(params)) {
    params = c(lapply(list_trees, FUN = mle_loglik, graf = adj), recursive = TRUE)
    #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
    params = params + -1*max(params)
    params = exp(params)/sum(exp(params))
  }
  pm = tree_av(params, list_trees)
  keep = upper.tri(pm, diag = FALSE)
  keep = keep * (pm > 0) * (pm < 1)
  
  res = 0
  
  for(k in 1:dim(adj)[3]) {
    res = res +
      sum(log(pm[keep*adj[,,k] == 1])) +
      sum(log(1-pm[keep*(1-adj[,,k]) == 1]))
  }
  return(res)
}


#' Compute the likelihood ratio: Null means model parameters are the same
#' 
#' @param adjs Observed adjacency matrix
#' @param tr Tree object (tree_structure)
#' 
#' @return Likelihood ratio
#' 
#' @export
#' 
compute_likratio = function(adjs, tr) {  
  null = mle_loglik(tr, adjs)
  #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
  alt = mle_loglik(tr,adjs[,,1, drop = FALSE]) + mle_loglik(tr, adjs[,,2, drop = FALSE])
  #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
  test_statistic = 2 * (alt - null)
  return(test_statistic)
}


#' Compute the likelihood ratio of the averaged model
#' 
#' @param adjs Observed adjacency matrices
#' @param trs Tree objects
#' @param av_method Averaging method: Accepts "avg" or "weighted" now
#' 
#' @return Likelihood ratio
#' 
#' @export
#' 
compute_likratio_mult = function(adjs, trs, av_method = "avg") {
  ntrees = length(trs)
  null_trees = list()
  alt1_trees = list()
  alt2_trees = list()
  for(j in 1:ntrees) {
    null_trees[[j]] = cvx_est(trs[[j]], adjs, wt = 0)
    alt1_trees[[j]] = cvx_est(trs[[j]], adjs[,,1, drop = FALSE], wt = 0)
    alt2_trees[[j]] = cvx_est(trs[[j]], adjs[,,2, drop = FALSE], wt = 0)
  }
  
  if (av_method == "avg") { 
    nullp = tree_av(rep(1, times = ntrees), list_trees = null_trees)
    altp1 = tree_av(rep(1, times = ntrees), list_trees = alt1_trees)
    altp2 = tree_av(rep(1, times = ntrees), list_trees = alt2_trees)
  } else if (av_method == "weighted") {
    weights = sapply(null_trees, function(x) compute_likelihood(base_tree = x, adj_mat = apply(adjs, c(1,2), sum)))
    weights = exp(weights - max(weights))
    
    nullp = tree_av(weights, list_trees = null_trees)
  
    altp1 = tree_av(weights, list_trees = alt1_trees)
    altp2 = tree_av(weights, list_trees = alt2_trees)
  } else { 
    return ("Error in av_method: nonexisting method")
  }
  
  nondiag_1s = (1 - diag(dim(adjs)[1]))
  
  obs_sum = apply(adjs, c(1,2), sum)
  null_LL = (2*nondiag_1s - obs_sum) * log(1 - nullp)
  diag(nullp) <- 1
  null_LL = null_LL + obs_sum * log(nullp)
  
  alt1_LL = (nondiag_1s - adjs[,,1]) * log(1 - altp1)
  diag(altp1) <- 1
  alt1_LL = alt1_LL + adjs[,,1] * log(altp1)
  
  alt2_LL = (nondiag_1s - adjs[,,2]) * log(1 - altp2)
  diag(altp2) <- 1
  alt2_LL = alt2_LL + adjs[,,2] * log(altp2)
  
  return(2 * sum((alt1_LL + alt2_LL - null_LL)))
  
}

