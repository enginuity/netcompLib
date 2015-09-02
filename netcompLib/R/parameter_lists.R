## Collection of functions to set parameters... 

## TODO: [Examine Default Parameters]

#' Set network model parameters
#' 
#' This function sets up the (list of) model parameters to be passed into model generation functions (eg. sets the max number of blocks in a SBM). There are default parameters that are used if this function is called with no arguments. 
#' 
#' @param Nnodes [int] :: Number of nodes in the network model
#' @param type [char; ALLOWED = c("none", "block", "tree", "latent", "random")] :: type of network model
#' \itemize{
#'  \item none -- Probably not useful; can be used to construct relatively empty model class for whatever purpose... 
#'  \item block -- Stochastic Block Model
#'  \item tree -- Hierarchical Random Graph
#'  \item latent -- Latent Space Model. Note that this type of model CANNOT generate random edge partitions since edge groups will be size 1 each (almost surely)
#'  \item random -- A model where each edge belongs to one of a number of classes, and the edge probability depends solely on this assignment. Should really only be used to generate random edge partitions. 
#' }
#' @param pmin [double] :: Minimal possible edge probability
#' @param pmax [double] :: Maximal possible edge probability
#' @param block_nclass [int] :: Number of blocks in block model
#' @param block_avgdensity [double] :: Set average density in block model (ignored if NULL)
#' @param block_assign [vector-int] :: block assignments
#' @param block_probs [matrix-double] :: block model probabilities
#' @param random_ngroups [int] :: Number of groups for completely random edge partition
#' @param tree_type [char; ALLOWED = c("left", "random", "original")] :: Randomization on structure of tree
#' \itemize{
#'  \item left -- This is a left-branching tree (all splits are on the left branch), but the nodes are in a random order
#'  \item random -- This is a tree where a random terminal node is chosen for a split. Tends to randomize into more balanced trees
#'  \item original -- This is the most balanced possible tree, with node order shuffled. 
#' }
#' @param latent_dim [int] :: Dimension of latent space in latent space models
#' @param latent_nclass [int] :: Number of clusters in latent space model
#' @param latent_sdcenter [double] :: SD on centers of latent space model
#' @param latent_isgennorm [logical] :: If TRUE: Uses normal distribution for latent locations. Otherwise, uses uniform distribution. 
#' 
#' @return [list] :: A list of parameters
#' 
#' @export
#' 
set_model_param = function(Nnodes = 30, type = 'block', pmin = 0.03, pmax = 0.97, block_nclass = 3, block_avgdensity = NULL, block_assign = NULL, block_probs = NULL, random_ngroups = 10, tree_type = "random", latent_dim = 3, latent_nclass = 3, latent_sdcenter = 5, latent_isgennorm = TRUE) {
  
  return(list(Nnodes = Nnodes, type = type, pmin = pmin, pmax = pmax, block_nclass = block_nclass, block_avgdensity = block_avgdensity, block_assign = block_assign, block_probs = block_probs, random_ngroups = random_ngroups, tree_type = tree_type, latent_dim = latent_dim, latent_nclass = latent_nclass, latent_sdcenter = latent_sdcenter, latent_isgennorm = latent_isgennorm))
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (set_sim_param)
#' Set up simulation parameters
#' 
#' This function creates a list of simulation parameters that are desired. 
#' 
#' @param cc_adj [vector-double] :: Amount of SE's away from the correlation estimate used (how conservative? 0 means no adjustment (and requires large sample for guarantees; larger values give a conservative p-value))
#' @param thres_ignore [vector-int] :: Ignore edge groups with fewer than this many edges
#' @param alphas [vector-double] :: Size(s) of the hypothesis test (probability of reject given true null)
#' @param n_models [vector-int] :: Number(s) of edge partitions to use for testing
#' @param pval_fx_names [] :: temp
#' @param pval_sim_null [] :: temp
#' @param recycle_fitstructs [logical] :: Should fitting structures be re-used for different number of structures used?
#' 
#' @return [list] :: A list of parameters, structured as follows
#' \itemize{
#' \item cc_adj -- [vector-double] :: Amount of SE's away from the correlation estimate used (how conservative? 0 means no adjustment (and requires large sample for guarantees; larger values give a conservative p-value))
#' \item thres_ignore -- [vector-int] :: Ignore edge groups with fewer than this many edges
#' \item alphas -- [vector-double] :: Size(s) of the hypothesis test (probability of reject given true null)
#' \item n_models -- [vector-int] :: Number(s) of edge partitions to use for testing
#' \item pval_adj -- [list] :: Information about p-value adjustment functions
#'   \itemize{
#'    \item fx -- [named_list-functions] :: This is a named list of functions
#'    \item simnull -- [vector-logical] :: This is a corresponding vector -- TRUE if the null distribution needs to be simulated.
#'   }
#' \item struct_needed -- [int] :: Total number of random structures needed
#' \item struct_indices -- [list-vector-int] :: Indices for the structure list for each iteration (each value in n_models)
#' }
#' 
#' @export
#' 
set_sim_param = function(cc_adj = c(0,2), thres_ignore = c(5, 10), alphas = 0.05, n_models = c(1,25,50,100), pval_fx_names = c("mult_bonferroni", "mult_highcrit", "mult_pearson"), pval_sim_null = c(FALSE, TRUE, TRUE), recycle_fitstructs = TRUE) {
  
  ## TODO: [Rename] n_models into n_structs or something
  
  ## Extract appropriate information from function inputs
  fxlist = list()
  for (j in seq_along(pval_fx_names)) {
    fxlist[[pval_fx_names[j]]] = eval(parse(text = pval_fx_names[j]))
  }
  
  ## Extract fitting structure indices
  if (recycle_fitstructs) { 
    struct_needed = max(n_models)
    struct_indices = lapply(n_models, function(x) { 1:x })
  } else {
    struct_needed = sum(n_models)
    temp = c(0, cumsum(n_models))
    struct_indices = lapply(seq_along(temp[-1]), function(x) { (temp[x]+1):temp[x+1] })
  }
  
  return(list(cc_adj = cc_adj, thres_ignore = thres_ignore, alphas = alphas, n_models = n_models,
              pval_adj = list(fx = fxlist, simnull = pval_sim_null), 
              struct_needed = struct_needed, struct_indices = struct_indices))
}

