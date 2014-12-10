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


#' Wrapper for compute_pval_multtrees, runs this for a random adjmat. Uses global variables. 
#' 
#' Requires existence of: fl_LIST, params, gen_model
#' where params is a list with: mode, cc_adj, thres_ignore, alphas, n_models, gen_model
#' -- gen_model is a list with: $mode, $m1, $m2
#' 
#' @param noparam DELETE THIS PARAMTER once i fix the auto param generation code. 
#' 
#' @return p-value matrix/array
#' 
#' @export
#' 
sim_one_set = function(noparam) {

  if (gen_model$mode == "tree") {
    adjmat = abind(cmn_network(gen_model$m1, 1), cmn_network(gen_model$m2,1))
  } else if (gen_model$mode == "block" | gen_model$mode == "blockmodel") {
    adjmat = abind(block_network(gen_model$m1, 1), block_network(gen_model$m2, 1))
  } else if (gen_model$mode == "latent") {
    adjmat = abind(latent_network(gen_model$m1, 1), latent_network(gen_model$m2, 1))
  }
  
  res = compute_pval_multtrees(
    adj1 = adjmat[,,1,drop = FALSE], adj2 = adjmat[,,2,drop = FALSE],
    mode = params$mode,
    pl = list(thres_ignore = params$thres_ignore,
              cc_adj = params$cc_adj, alphas = params$alphas),
    fl_list = fl_LIST, n_models = params$n_models)
  
  return(res)
}

