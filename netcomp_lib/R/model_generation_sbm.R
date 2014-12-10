##@S Code to generate 

## TODO: [Obselete]

#' Generate block model parameters
#' 
#' @param NN Number of nodes
#' @param K Number of classes
#' @param pmin Probability minimum
#' @param pmax Probability maximum
#' 
#' @return List with block model information
#' 
#' @export
#' 
gen_block = function(NN, K, pmin = 0, pmax = 1) {
  group_assign = sample(1:K, size = NN, replace = TRUE)
  prob_matrix = matrix(0, nrow = K, ncol = K)
  for(j in 1:K) { for(i in 1:K) {
    if (i <= j) {
      prob_matrix[i,j] = runif(1, pmin, pmax)
      if (i != j) {prob_matrix[j,i] = prob_matrix[i,j]}
    }
  }}
  
  return(list(assign = group_assign, probmat = prob_matrix))
}


#' Generate data from a block model
#' 
#' @param bm Block model information
#' @param n Number of observations
#' 
#' @return Array of network observations
#' 
#' @export
#' 
block_network = function(bm, n = 1) {
  NN= length(bm$assign)
  res = array(0, dim = c(NN, NN, n))
  for(d in 1:n) {
    for(k in 1:NN) {
      for(j in 1:NN) {
        if (k < j) {
          v = rbinom(n = 1, size = 1, prob = bm$probmat[bm$assign[k], bm$assign[j]])
          res[j,k,d] = v
          res[k,j,d] = v
        }
      }}
  }
  return(res)
}


