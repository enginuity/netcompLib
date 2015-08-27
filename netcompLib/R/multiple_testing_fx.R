##@S Functions that compute multiple testing p-values (or a statistic)

#' Compute higher-critiism test statistic
#' 
#' @param x [vector-double] :: P-values
#' @param alpha_min [double] :: Minimal alpha to maximize over
#' @param alpha_max [double] :: Maximal alpha to maximize over
#' 
#' @return [double] :: Test statistic
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
#' @param x [vector-double] :: P-values
#' @param adjust [logical] :: If NULL, function returns test statistic. Otherwise, this takes on a multiplier to variance. 
#' 
#' @return [double] :: Adjusted test statistic
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
#' @param x [vector-double] :: P-values
#' 
#' @return [double] :: Minimal p-value
#' 
#' @export
#' 
mult_bonferroni = function(x) {
  return(min(x))
}



