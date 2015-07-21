
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
  ## TODO: [Obselete?] -- only used in fast_compute_pval2 and its friends
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



#' Different version of computing the df adjustment
#' 
#' @param n cell sizes
#' @param cell_corr temp
#' @param cc_adj Number of SE's away for conservativeness in cell_correlation estimate
#' 
#' @return Df adjustment
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

