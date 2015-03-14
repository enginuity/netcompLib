# Defines the NetworkStructHRG class
setClass("NetworkStructHRG", representation(tree_list = "list", expand = "list", counts = "numeric"), contains = "NetworkStruct")

#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkStructHRG = function(Nnodes = 10, model_param = set_model_param()) {
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelHRG(Nnodes = Nnodes, model_param = model_param)
  return(extractStruct(NetM))
}







# Generic Function Definitions --------------------------------------------

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructHRG)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructHRG = function(NetM) {
  return("tree")
}
setMethod("getNetType", signature(NetM = "NetworkStructHRG"), getNetType.NetworkStructHRG)


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructHRG)
#' <What does this function do>
#' 
#' @param NetS temp
#' @param adja1 temp
#' @param adja2 temp
#' @param Nobs temp
#' @param pl temp
#' 
#' @return temp
#' 
#' @export
#' 
computePval.NetworkStructHRG = function(NetS, adja1, adja2, Nobs, pl) {
  if (FALSE) {
    NetM = NetworkModel(Nnodes = 30, type = "block")
    adja1 = sampleNetwork(NetM)
    adja2 = sampleNetwork(NetM)
    Nobs = 1
    pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_models = c(1,20))
    NetS = NetworkStructHRG(Nnodes = 30, model_param = set_model_param(block_nclass = 3))
    load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
  }
  
  ## Compute total count for each edge (and edgesumc is the sum of products, to be used in computing cell-wise correlations)
  if (Nobs > 1) {
    edgesum1 = apply(adja1, c(1,2), sum); edgesum2 = apply(adja2, c(1,2), sum); 
    edgesumc = apply(adja1*adja2, c(1,2), sum)
  } else if (Nobs == 1) {
    if (length(dim(adja1)) == 2) { edgesum1 = adja1 } else { edgesum1 = adja1[,,1] }
    if (length(dim(adja2)) == 2) { edgesum2 = adja2 } else { edgesum2 = adja2[,,1] }
    edgesumc = edgesum1 * edgesum2
  }
  
  ## Do counting
  obs1_count = sapply(NetS@expand, function(x) { sum(edgesum1[x[[1]], x[[2]]]) })
  obs2_count = sapply(NetS@expand, function(x) { sum(edgesum2[x[[1]], x[[2]]]) })
  obsc_count = obs1_count + obs2_count
  obsp_count = sapply(NetS@expand, function(x) { sum(edgesumc[x[[1]], x[[2]]]) })
  cell_sizes = NetS@counts * Nobs
  
  
  ## Compute MLE Estimates
  mle_p1 = obs1_count / cell_sizes
  mle_p2 = obs2_count / cell_sizes
  mle_pc = obsc_count / cell_sizes / 2
  
  ## Compute cell-correlation sample estimates
  mle_pxy = obsp_count / cell_sizes
  num = mle_pxy - (mle_p1 * mle_p2)
  den = sqrt(mle_p1 * (1 - mle_p1) * mle_p2 * (1 - mle_p2))
  cell_corrs = num/den
  
  ## Compute cellwise log-likelihoods
  LL_null = compute_loglik_fromPC(x = obsc_count, n = 2 * cell_sizes, p = mle_pc)
  LL_alt1 = compute_loglik_fromPC(x = obs1_count, n = cell_sizes, p = mle_p1)
  LL_alt2 = compute_loglik_fromPC(x = obs2_count, n = cell_sizes, p = mle_p2)
  
  cellwise_TS = -2 * (LL_null - LL_alt1 - LL_alt2)
  
  pval_matrix = matrix(-1, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
  df_adjs = sapply(pl$cc_adj, 
                   function(x) {compute_df_adjustment2(n = cell_sizes, cell_corr = cell_corrs, cc_adj = x)})
  for(i in seq_along(pl$thres_ignore)) {
    to_keep = which(cell_sizes >= pl$thres_ignore[i])
    csq = sum(cellwise_TS[to_keep])
    dfs = apply(df_adjs[to_keep,,drop = FALSE], 2, sum)
    pval_matrix[,i] = pchisq(csq, dfs, lower.tail = FALSE)
  }
  
  return(pval_matrix)  
}
setMethod("computePval", signature(NetS = "NetworkStructHRG"), computePval.NetworkStructHRG)



