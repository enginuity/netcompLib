## Helper functions

#' Compute the loglikelihood from probabilities and counts
#' 
#' @param x [matrix/array-int] :: observed counts in each cell
#' @param n [int] :: Number observations
#' @param p [vector-double] :: Corresponding estimated cell probabilities (is just x / n)
#' 
#' @return [vector-double] :: Vectorized version of log-likelihood (per edge group)
#' 
#' @export
#' 
compute_cellwise_loglik = function(x, n, p) {
  res = x * log(p) + (n - x) * log(1 - p)
  diag(res) = NA
  
  ## If probability is 0, but also observe no dyad, log-lik = 0. Likewise if probability is 1, but observe all dyads
  bad0s = which(p == 0)
  bad1s = which(p == 1)
  res[c(bad0s, bad1s)] = 0
  
  ## Else, likelihood of data is 0 (so log-lik is -Inf)
  res[bad0s][x[bad0s] != 0] = -Inf
  res[bad1s][n != x[bad1s]] = -Inf
  
  return(res)
}


compute_small_samp_dfadj = function(n, p) {
  
  compute_exp1 = function(n, p) {
    x = 1:n # Ignore x = 0 since that gives a log-likelihood of 0 (likelihood of 1)
    exp_split = dbinom(x, n, p) * x * log(x/n)
    return(sum(exp_split))
  }
  
  compute_expsym = function(n,p) { return(compute_exp1(n,p) + compute_exp1(n, 1-p)) }  
  
  compute_expLLR = function(n, p) {
    Enull_ll = compute_expsym(2*n, p) 
    Ealt_ll = 2*compute_expsym(n, p) 
    Ellr = -2 * (Enull_ll - Ealt_ll)
    return(Ellr)
  }
  
  return(compute_expLLR(n,p))
}


#' Compute the loglikelihood from probabilities and counts
#' 
#' @param x [vector-int] :: observed counts in each edge group
#' @param n [vector-int] :: Corresponding size of edge groups (cells)
#' @param p [vector-double] :: Corresponding estimated cell probabilities (is just x / n)
#' 
#' @return [vector-double] :: Vectorized version of log-likelihood (per edge group)
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



#' Compute the df adjustment
#' 
#' Degree of freedom adjustment for the likelihood ratio test due to mis-specified model, and also due to small-sample issues
#' 
#' @param n [vector-int] :: Edge group sizes
#' @param cell_corr [vector-double] :: Estimated within-edge-group correlations
#' @param cc_adj [double] :: Number of SE's away for conservativeness in cell_correlation estimate
#' 
#' @return [vector-double] :: Vector of df adjustments (per-edgegroup)
#' 
#' @export
#' 
compute_df_adjustment2 = function(n, cell_corr, cc_adj = 2) {
  ## TODO: [Rename] Fix the name of this function. 
  
  ## Assumes null is true:                                       
  
  ## Cell-correlation adjustment #####################
  ## Attach standard error to cell-correlations
  cell_est_se = 1/sqrt(n)
  low_est_corr = cell_corr - cc_adj * cell_est_se
  
  ## Threshold these at least 0. 
  cor_adjust = sapply(low_est_corr, function(x) {ifelse(is.na(x), 0, max(0, x))})
  
  ## Small-sample correction ########################
  ## Apply small-sample correction
  ss_adjust = sapply(n, function(x) {
    ifelse((x > 265), 1.01, small_samp_dfcorr[x]) } )
  
  ## Combine corrections
  df_adj = ss_adjust - cor_adjust
  return(df_adj)
}


#' Checks vector for valid probabilities
#' 
#' @param x [vector-double] :: Vector to check for probabilities
#' 
#' @return [vector-logical] :: ith entry is TRUE if x[i] is a valid probability
#' 
#' @export
#' 
is_prob = function(x) { 
  return(x >= 0 & x <= 1)
}
