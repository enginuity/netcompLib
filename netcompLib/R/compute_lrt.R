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



# Functions relating to cell-wise correlation adjustment ------------------

#' Compute true correlation adjustment (df adjustment)
#' 
#' @param btree base tree
#' @param ftree fixed tree
#' 
#' @return correlation adjustment (cell-wise)
#' 
#' @export
#' 
true_cor = function(btree, ftree) {
  ep = edge_probs(btree)
  N = btree$nodes
  est_struct = closest_ancestor(ftree)$anc.table - N
  diag(est_struct) <- 0
  
  cors = rep(0, times = N-1)
  for(j in 1:(N-1)) {
    ps = ep[est_struct == j]
    mp = mean(ps)
    cors[j] = (mean(ps^2) - mp^2) / (mp - mp^2)
  }
  return(cors)
}


#' Estimates cell-wise correlations (for df adjustment)
#' 
#' @param ftree Fitting fixed tree
#' @param obsadjm1 Adj matrix 1
#' @param obsadjm2 Adj matrix 2
#' @param exp_child Expanded_child of tree
#' 
#' @return Cell-wise df corrections from correlations
#' 
#' @export
#' 
compute_corrs = function(ftree, obsadjm1, obsadjm2, exp_child = NULL) {
  ## computes cell-wise correlations for tree, using obsadjm1/obsadjm2. these can be arrays...
  
  if(is.null(exp_child)) {
    ec = expanded_children_from_tree(ftree)
  } else {
    ec = exp_child
  }
  N = ftree$nodes
  res = rep(NA, times = 2*N - 1)
  for(j in (N+1):(2*N-1)) {
    
    res[j] = cor(
      as.vector(obsadjm1[ec[[j-N]][[1]],ec[[j-N]][[2]],]),
      as.vector(obsadjm2[ec[[j-N]][[1]],ec[[j-N]][[2]],]))
  }
  return(res)
}



# Code to perform hypothesis test -----------------------------------------


#' Estimate the needed df adjustment given an adjacency matrix and a tree structure
#' 
#' @param ftree Fixed tree structure
#' @param obsadjs Observed adjacency array
#' @param exp_child Expanded_child
#' 
#' @return DF adjustments
#' 
#' @export
#' 
compute_df_adjustment = function(ftree, obsadjs, exp_child = NULL) {
  ## Takes the two observed edge matrix, and computes the df correction (cell-wise)
  
  ## ftree = bt
  ## ftree = fts[[1]]
  ## obsadjs = adjmat
  
  est = mle_est(ftree, obsadjs, exp_child)
  
  ## Cell-counts
  counts = est$counts / 2
  
  ## Estimate cell-correlations
  cell_corr = compute_corrs(ftree, obsadjs[,,1, drop = FALSE], obsadjs[,,2, drop = FALSE], exp_child)[-1:-ftree$nodes]
  
  cell_est_se = 1/sqrt(counts)
  ## low estimate of correlation is the cell_correlation - 2 * SE
  low_est_corr = cell_corr - 2 * cell_est_se
  
  ## Threshold this (negative correlations are not allowed)
  cor_adjust = sapply(low_est_corr, function(x) {ifelse(is.na(x), 0, max(0, x))})
  
  ## Apply small-sample correction
  ss_adjust = sapply(counts, function(x) {
    ifelse((x > 265), 1.01, small_samp_dfcorr[x]) } )
  
  ## Combine corrections
  df_adj = ss_adjust - cor_adjust
  
  return(df_adj)
}


#' Compute the p-value given a fixed tree structure and observed adjacency matrices
#' 
#' @param ftree Fixed testing tree structure
#' @param obsadjs Observed adjacency matrix
#' @param thres_ignore Ignore all edge groups with these many values (or fewer)
#' @param exp_child expanded children of fixed tree structure
#' 
#' @return p-value
#' 
#' @export
#' 
compute_pval = function(ftree, obsadjs, thres_ignore = 2, exp_child = NULL) {
  ## This function runs the likelihood ratio test, applying the df adjustment, and ignoring all cells with 'thres_ignore' data points or fewer.
  
  ## Compute loglikelihoods
  null_ll = mle_loglik(ftree, obsadjs, cellwise = TRUE, exp_child = exp_child)
  #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
  alt1_ll = mle_loglik(ftree, obsadjs[,,1, drop = FALSE], cellwise = TRUE, exp_child = exp_child)
  #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
  alt2_ll = mle_loglik(ftree, obsadjs[,,2, drop = FALSE], cellwise = TRUE, exp_child = exp_child)
  #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
  
  ## Compute df adjustment
  df_adj = compute_df_adjustment(ftree, obsadjs, exp_child = exp_child)
  
  ## Find cell-counts, threshold
  counts = mle_est(ftree, obsadjs[,,1, drop = FALSE], exp_child = exp_child)$counts
  to_keep = which(counts > thres_ignore)
  
  ## Compute test statistic, and appropriate conservative df. 
  cellwise_TS = -2 * (null_ll - alt1_ll - alt2_ll)
  csq = sum(cellwise_TS[to_keep])
  df = sum(df_adj[to_keep])
  
  ## Compute p-value and return
  pval = pchisq(csq, df, lower.tail = FALSE)
  return(pval)
}


# Faster versions of hypothesis testing -----------------------------------


## Helper functions


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (compute_samp_corrs)
#' <What does this function do>
#' 
#' @param pc temp
#' @param p1 temp
#' @param p2 temp
#' 
#' @return temp
#' 
#' @export
#' 
compute_samp_corrs = function(pc, p1, p2) {
  ## Update of 'compute_correlation_fromprobs', for lists instead of simple vectors.
  ## If input are not lists, convert into a list. 
  if (any(c(!is.list(pc), !is.list(p1), !is.list(p2)))) {
    pc = list(pc); p1 = list(p1); p2 = list(p2)
  }
  res = lapply(seq_along(pc), function(i) {
    num = pc[[i]] - (p1[[i]] * p2[[i]])
    den = sqrt(p1[[i]] * (1 - p1[[i]]) * p2[[i]] * (1 - p2[[i]]))
    return(num/den)
  })
  return(res)
}



#' Compute the cell-correlation given cell probabilities
#' 
#' @param pc Probability of both 1s
#' @param p1 Probability of 1s in adjmatrix 1
#' @param p2 Probability of 1s in adjmatrix 2
#' 
#' @return Correlation estimate
#' 
#' @export
#' 
compute_correlation_fromprobs = function(pc, p1, p2) {
  # pc, p1, p2 are all vectors; applies this vectorized
  numer = pc - (p1 * p2)
  denom = sqrt(p1 * (1 - p1) * p2 * (1 - p2))
  return(numer/denom)
}


#' Compute the loglikelihood from probabilities and counts
#' 
#' @param x observed count
#' @param n Size of edge groups (cells)
#' @param p Estimated cell probabilities (is just x / n)
#' 
#' @return Log likelihood
#' 
#' @export
#' 
compute_loglik_fromPC = function(x, n, p) {
  ## x is the observed count; n is total count, p is x / n
  ## this function returns vectorized output
  res = x * log(p) + (n - x) * log(1 - p)
  res[is.nan(res)] = 0
  return(res)
}

