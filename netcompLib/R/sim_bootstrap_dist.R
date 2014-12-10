##@S Functions for simulation of some network distances using bootstrap



#' Generate a boostrapped distance
#' 
#' @param rngseed Seed for random number generator
#' @param fix_tree_est Fixed tree structure
#' 
#' @return Simulated distance measure
#' 
#' @export
#' 
gen_boot_dists = function(rngseed, fix_tree_est) {
  set.seed(rngseed)
  adjs = cmn_network(fix_tree_est, 2)
  est_dist = tree_distance(cvx_est(fix_tree_est, adjs[,,1,drop=FALSE], 5),
                           cvx_est(fix_tree_est, adjs[,,2,drop= FALSE], 5))
  return(est_dist)
}


#' Generate bootstrap distribution p-value
#' 
#' @param seed Random number generation seed
#' @param fix_tree_est Estimating tree
#' 
#' @return List of boostrapped distances & true distance
#' 
#' @export
#' 
gen_boot_dist_pval = function(seed, fix_tree_est) {
  set.seed(seed)
  adjs = cmn_network(fix_tree_est,2)
  boot_dist = tree_distance(cvx_est(fix_tree_est, adjs[,,1,drop=FALSE], 5),
                            cvx_est(fix_tree_est, adjs[,,2,drop= FALSE], 5))
  
  new_comb_est = cvx_est(fix_tree_est, adjs, 5)
  bootp = sapply(sample(1:500000, size = 100), function(x) {gen_boot_dists(x, new_comb_est)})
  return(list(bootdist = bootp, true_distance = boot_dist))
}


