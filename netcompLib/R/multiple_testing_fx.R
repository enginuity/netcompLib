##@S Functions that compute multiple testing p-values (or a statistic)

#' Compute higher-critiism test statistic
#' 
#' @param x vector of p-values
#' @param alpha_min Minimal alpha to maximize over
#' @param alpha_max Maximal alpha to maximize over
#' 
#' @return Test statistic
#' 
#' @export
#' 
mult_highcrit = function(x, alpha_min = 10^(-8), alpha_max = 0.2) {
  x[x < alpha_min] = alpha_min
  m = ceiling(alpha_max * length(x))
  n = length(x)
  xsub = sort(x, decreasing = FALSE)[1:m]
  stat = max(sqrt(n) * abs(seq_along(xsub) / n - xsub)/sqrt(xsub * (1 - xsub)))
  
  return(stat)
}


#' Computes pearson-like test statistic (or p-value)
#' 
#' Computes p-value by assuming normality (using larger variance estimate)
#' Larger variance estiamte is bounded by 1/n + cor, but uses 1/n * adjust as variance. 
#' 
#' @param x Vector of pvalues
#' @param adjust If NULL, function returns test statistic. Otherwise, this takes on a multiplier to variance. 
#' 
#' @return temp
#' 
#' @export
#' 
mult_pearson = function(x, adjust = NULL) { 
  stat = sum(-1 * log(x)) / length(x)
  if (is.null(adjust)) { return(stat) }
  return(1 - pnorm(q = ( (stat - 1) / sqrt(1/n * adjust))))
}


#' Computes Bonferroni correction
#' 
#' @param x Vector of p-values
#' 
#' @return Returns minimal p-value
#' 
#' @export
#' 
mult_bonferroni = function(x) {
  return(min(x))
}



