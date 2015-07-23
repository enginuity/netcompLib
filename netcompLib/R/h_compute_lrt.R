## Helper functions


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

