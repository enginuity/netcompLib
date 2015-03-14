## Collection of functions to set parameters... 


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (set_model_param)
#' Create a list of model parameters
#' 
#' This function sets up the model parameters to be passed into model generation functions (eg. sets the max number of blocks in a SBM). There are default parameters that are used if this function is called with no arguments. 
#' 
#' @param pmin Minimal possible edge probability
#' @param pmax Maximal possible edge probability
#' @param block_nclass Number of blocks in block model
#' @param block_avgdensity Set average density in block model (ignored if NULL)
#' @param block_assign temp
#' @param block_probs temp
#' @param random_ngroups Number of groups for completely random edge partition
#' @param tree_type Randomization on structure of tree: can be "left", "random", or "original"
#' @param latent_dim Dimension of latent space in latent space models
#' @param latent_nclass Number of clusters in latent space model
#' @param latent_sdcenter SD on centers of latent space model
#' @param latent_isgennorm If TRUE: Uses normal distribution for latent locations. Otherwise, uses uniform distribution. 
#' 
#' @return Return list of parameters
#' 
#' @export
#' 
set_model_param = function(pmin = 0, pmax = 1, block_nclass = 3, block_avgdensity = NULL, block_assign = NULL, block_probs = NULL, random_ngroups = 10, tree_type = "random", latent_dim = 3, latent_nclass = 3, latent_sdcenter = 5, latent_isgennorm = TRUE) {
  return(list(pmin = pmin, pmax = pmax, block_nclass = block_nclass, block_avgdensity = block_avgdensity, block_assign = block_assign, block_probs = block_probs, random_ngroups = random_ngroups, tree_type = tree_type, latent_dim = latent_dim, latent_nclass = latent_nclass, latent_sdcenter = latent_sdcenter, latent_isgennorm = latent_isgennorm))
}


#' Creates a list of simulation parameters
#' 
#' @param cc_adj Amount of SE's away from the correlation estimate used (how conservative? 0 means no adjustment (and requires large sample for guarantees; larger values give a conservative p-value))
#' @param thres_ignore Ignore edge groups with fewer than this many edges
#' @param alphas Significance level of test
#' @param n_models Number of edge partitions to use for testing
#' 
#' @return A list of parameters (often asked for as pl in functions)
#' 
#' @export
#' 
set_sim_param = function(cc_adj = c(0,2), thres_ignore = c(5, 10), alphas = 0.05, n_models = c(1,25,50,100)) {
  return(list(cc_adj = cc_adj, thres_ignore = thres_ignore, alphas = alphas, n_models = n_models))
}

  
