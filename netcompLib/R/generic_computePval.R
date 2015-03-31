##@S Generic function that performs the hypothesis test (computes the p-value for the likelihood ratio test)

## TODO: [Fully Documented] (remove this marking eventually)
## TODO: [Implement] The 'mode' parameter only works for SBMs right now. 

setGeneric("computePval", function(NetS, adja1, adja2, Nobs, pl, mode) standardGeneric("computePval"))

#' Computes p-value for likelihood ratio test
#' 
#' This function takes the input network structure (or a list of network structures), and computes the likelihood ratio test for the two input adjacency matrices (or arrays). 
#' 
#' @param NetS Input NetworkStruct object
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of network observations per class (default = 1)
#' @param pl A list of parameters, set by set_sim_param
#' @param mode How to output results? 'default' gives the standard p-values; 'nodewise' gives the chi-square contributions per node
#' 
#' @return A matrix (or a list of matrices) of p-values (depending on the testing parameters)
#' 
#' @export
#' 
computePval = function(NetS, adja1, adja2, Nobs, pl, mode) {
  stop("Placeholder for documentation purposes")
}


# Define new generics -----------------------------------------------------

#' Computes p-value for likelihood ratio test
#' 
#' This function takes the input network structure (or a list of network structures), and computes the likelihood ratio test for the two input adjacency matrices (or arrays). 
#' 
#' @param NetS Input NetworkStruct object
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of network observations per class (default = 1)
#' @param pl A list of parameters, set by set_sim_param
#' @param mode How to output results? 'default' gives the standard p-values; 'nodewise' gives the chi-square contributions per node
#' 
#' @return A matrix (or a list of matrices) of p-values (depending on the testing parameters)
#' 
#' @export
#' 
computePval.NetworkStruct = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  stop("Not implemented for this case")
}



#' Computes p-value for likelihood ratio test
#' 
#' This function takes the input network structure (or a list of network structures), and computes the likelihood ratio test for the two input adjacency matrices (or arrays). 
#' 
#' @param NetS Input NetworkStruct object
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of network observations per class (default = 1)
#' @param pl A list of parameters, set by set_sim_param
#' @param mode How to output results? 'default' gives the standard p-values; 'nodewise' gives the chi-square contributions per node
#' 
#' @return A list of matrices of p-values (depending on the testing parameters)
#' 
#' @export
#' 
computePval.NetworkStructList = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  res = lapply(NetS@models, function(x) { computePval(x, adja1, adja2, Nobs, pl, mode = mode) } )
  return(res)
}


#' Computes p-value for likelihood ratio test
#' 
#' This function takes the input network structure (or a list of network structures), and computes the likelihood ratio test for the two input adjacency matrices (or arrays). 
#' 
#' @param NetS Input NetworkStruct object
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of network observations per class (default = 1)
#' @param pl A list of parameters, set by set_sim_param
#' @param mode How to output results? 'default' gives the standard p-values; 'nodewise' gives the chi-square contributions per node
#' 
#' @return A matrix (or a list of matrices) of p-values (depending on the testing parameters)
#' 
#' @export
#' 
computePval.NetworkStructRND = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  if (FALSE) {
    NetM = NetworkModel(Nnodes = 30, type = "block")
    adja1 = sampleNetwork(NetM)
    adja2 = sampleNetwork(NetM)
    Nobs = 1
    pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_models = c(1,20))
    NetS = NetworkStructRND(Nnodes = 30, model_param = set_model_param(block_nclass = 3))
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
  vedgesum1 = as.vector(edgesum1); vedgesum2 = as.vector(edgesum2); vedgesumc = as.vector(edgesumc)
  obs1_count = sapply(NetS@ids, function(x) { sum(vedgesum1[x])} )
  obs2_count = sapply(NetS@ids, function(x) { sum(vedgesum2[x])} )
  obsc_count = obs1_count + obs2_count
  obsp_count = sapply(NetS@ids, function(x) { sum(vedgesumc[x])} )
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


#' Computes p-value for likelihood ratio test
#' 
#' This function takes the input network structure (or a list of network structures), and computes the likelihood ratio test for the two input adjacency matrices (or arrays). 
#' 
#' @param NetS Input NetworkStruct object
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of network observations per class (default = 1)
#' @param pl A list of parameters, set by set_sim_param
#' @param mode How to output results? 'default' gives the standard p-values; 'nodewise' gives the chi-square contributions per node
#' 
#' @return A matrix (or a list of matrices) of p-values (depending on the testing parameters)
#' 
#' @export
#' 
computePval.NetworkStructHRG = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
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


#' Computes p-value for likelihood ratio test
#' 
#' This function takes the input network structure (or a list of network structures), and computes the likelihood ratio test for the two input adjacency matrices (or arrays). 
#' 
#' @param NetS Input NetworkStruct object
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of network observations per class (default = 1)
#' @param pl A list of parameters, set by set_sim_param
#' @param mode How to output results? 'default' gives the standard p-values; 'nodewise' gives the chi-square contributions per node
#' 
#' @return A matrix (or a list of matrices) of p-values (depending on the testing parameters)
#' 
#' @export
#' 
computePval.NetworkStructSBM = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  if (FALSE) {
    NetM = NetworkModel(Nnodes = 30, type = "block")
    adja1 = sampleNetwork(NetM)
    adja2 = sampleNetwork(NetM)
    Nobs = 1
    pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_models = c(1,20))
    NetS = NetworkStructSBM(Nnodes = 30, model_param = set_model_param(block_nclass = 3))
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
  COR = NetS@correct # correction factor for block models
  obs1_count = sapply(NetS@expand, function(x) { sum(edgesum1[x[[1]], x[[2]]], na.rm = TRUE) }) / COR
  obs2_count = sapply(NetS@expand, function(x) { sum(edgesum2[x[[1]], x[[2]]], na.rm = TRUE) }) / COR
  obsc_count = obs1_count + obs2_count
  obsp_count = sapply(NetS@expand, function(x) { sum(edgesumc[x[[1]], x[[2]]], na.rm = TRUE) }) / COR
  cell_sizes = sapply(NetS@expand, function(x) { sum(!is.na(edgesumc[x[[1]], x[[2]]])) }) / COR
  
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
  
  ## Compute edgewise LL's if desired
  if (mode == "nodewise") {
    ## TODO: Make this faster if necessary. 
    ncs = rep(0, times = nrow(edgesum1))
    mlenull = 0 * edgesum1
    mlealt1 = 0 * edgesum1
    mlealt2 = 0 * edgesum1
    for(pp in seq_along(NetS@expand)) {
      mlenull[NetS@expand[[pp]][[1]], NetS@expand[[pp]][[2]]] = (mle_pc[pp])
      mlealt1[NetS@expand[[pp]][[1]], NetS@expand[[pp]][[2]]] = (mle_p1[pp])
      mlealt2[NetS@expand[[pp]][[1]], NetS@expand[[pp]][[2]]] = (mle_p2[pp])
      mlenull[NetS@expand[[pp]][[2]], NetS@expand[[pp]][[1]]] = (mle_pc[pp])
      mlealt1[NetS@expand[[pp]][[2]], NetS@expand[[pp]][[1]]] = (mle_p1[pp])
      mlealt2[NetS@expand[[pp]][[2]], NetS@expand[[pp]][[1]]] = (mle_p2[pp])
    }
    llnull = edgesum1 * log(mlenull) + (Nobs - edgesum1) * log(1 - mlenull) + edgesum2 * log(mlenull) + (Nobs - edgesum2) * log(1 - mlenull)
    llalt1 = edgesum1 * log(mlealt1) + (Nobs - edgesum1) * log(1 - mlealt1)
    llalt2 = edgesum2 * log(mlealt2) + (Nobs - edgesum2) * log(1 - mlealt2)
    diag(llnull) = 0; diag(llalt1) = 0; diag(llalt2) = 0
    
    res = -2 * (llnull - llalt1 - llalt2)
    ncs = apply(res, 1, sum)
    
    return(list(pvals = pval_matrix, nodecontrib = ncs))
  } else if (mode == "default") {
    return(pval_matrix)  
  }
}





# setMethod ---------------------------------------------------------------
setMethod("computePval", signature(NetS = "NetworkStruct"), computePval.NetworkStruct)
setMethod("computePval", signature(NetS = "NetworkStructList"), computePval.NetworkStructList)
setMethod("computePval", signature(NetS = "NetworkStructSBM"), computePval.NetworkStructSBM)
setMethod("computePval", signature(NetS = "NetworkStructRND"), computePval.NetworkStructRND)
setMethod("computePval", signature(NetS = "NetworkStructHRG"), computePval.NetworkStructHRG)