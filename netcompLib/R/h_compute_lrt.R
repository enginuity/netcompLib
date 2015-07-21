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

