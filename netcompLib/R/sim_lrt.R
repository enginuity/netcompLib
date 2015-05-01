##@S Functions for performing the LRT simulation


#' Obtain simulated distribution of the loglikelihood ratio
#' 
#' @param btree Generating tree structure
#' @param ftree Fitting tree structure
#' @param Nsim Number of simulations
#' 
#' @return Vector of log likelihood ratios
#' 
#' @export
#' 
sim_mle_loglik_ratio = function(btree, ftree, Nsim = 1) {
  graf2 = cmn_network(btree,2*Nsim)
  tree = ftree
  sameparam = mle_loglik(tree,graf2)
  #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
  diffparam = mle_loglik(tree,graf2[,,1:Nsim,drop = FALSE]) + mle_loglik(tree, graf2[,,(Nsim+1):(2*Nsim),drop = FALSE])
  #|----##Need to Use updated version of mle_loglik in lrt_functions.R; more parameters to pass in --Mon Sep 15 02:36:03 2014--
  return(-2 * (sameparam - diffparam))
}

