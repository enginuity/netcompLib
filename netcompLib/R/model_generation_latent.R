
## TODO: [Obselete]

#' Generate model parameters for a latent space model. 
#' 
#' @param NN Number of nodes
#' @param K Number of centers for latent space model
#' @param gen_norm T/F: Use normal distribution to generate latent locations (otherwise, uniform)
#' @param D Dimensionality of latent space
#' @param sd_center SD of location distribution
#' 
#' @return Latent space model object
#' 
#' @export
#' 
gen_latent = function(NN, K,gen_norm = TRUE, D = 2, sd_center = 5) {
  group_assign = sample(1:K, size = NN, replace = TRUE)
  centers = list()
  
  for(j in 1:K) {
    if (gen_norm) {
      centers[[j]] = rnorm(n = D, mean = 0, sd = sd_center)
    } else {
      centers[[j]] = runif(n = D, min = -2 * sd_center, max = 2 * sd_center)
    }
  }
  
  locs = matrix(nrow= NN, ncol = D)
  for(j in 1:NN) {
    locs[j,] = centers[[group_assign[j]]] + rnorm(n = D, mean = 0, sd = 1)
  }
  return(list(locs = locs, alpha = runif(1, 0, 10)))
}


#' Simulate networks from a latent space model object. 
#' 
#' @param lm Latent space model object
#' @param n Number of networks drawn
#' 
#' @return Array of networks
#' 
#' @export
#' 
latent_network = function(lm, n = 1) {
  NN = nrow(lm$locs)
  odds = exp(lm$alpha - dist(lm$locs))
  prob_mat = matrix(nrow = NN, ncol = NN)
  prob_mat[lower.tri(prob_mat)] = odds / (1 + odds)
  
  res = array(0, dim = c(NN, NN, n))
  for(d in 1:n) {
    for(k in 1:NN) {
      for(j in 1:NN) {
        if (k < j) {
          v = rbinom(n = 1, size = 1, prob = prob_mat[j,k])
          res[j,k,d] = v
          res[k,j,d] = v
        }
      }}
  }
  return(res) 
}



