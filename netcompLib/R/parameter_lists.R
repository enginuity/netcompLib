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

set_sim_param = function(cc_adj = c(0,2), thres_ignore = c(5, 10), alphas = 0.05, n_models = c(1,25,50,100)) {
  
  ## pl is a list of parameters, and the constitutient values could be vectors if more than one is to be tested. This returns a matrix of pvalues
  ## -- thres_ignore: Ignore all cells concern a smaller number of edges than this value. Set to 0 to use every single cell. 
  ## -- cc_adj: The amount of SE's away from correlation estimate used.
  ## -- alphas: The cutoffs for p-values to be used
  ## -- -- this is set to 2 by default; setting to 0 -> simply use correlation estimate (thresholded to 0)
  return(list(cc_adj = cc_adj, thres_ignore = thres_ignore, alphas = alphas, n_models = n_models))
}

  
