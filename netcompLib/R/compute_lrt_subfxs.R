## Revision of LRT computations. Hopefully is faster! Is list-ified for new parametrization of models. 


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (fast_compute_pval_v2)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (fast_compute_pval_v2)
#' Runs test for a set of fixed trees. Returns matrix of p-values based on parameter settings. 
#' 
#' @param adja1 Adjacency array 1
#' @param adja2 Adjacency array 2
#' @param Nobs temp
#' @param pl Parameter list
#' @param fit_models Fixed models list for testing
#' @param return_chisq temp
#' 
#' @return Matrix of p-values (based on parameter list)
#' 
#' @export
#' 
fast_compute_pval_v2 = function(adja1, adja2, Nobs, pl, fit_models, return_chisq = FALSE) {
  ## Every time this function needs adjustment, adjust one of the following tester functions; can't really write a test for this large function. 
  
  ## TODO: [Obselete] This function obseletes fast_compute_pval. 
  
  ########## slow_compute_general_cellwise_mles ##########
  Nmodels = fit_models$Nmodels
  
  ## Compute total count for each edge (and edgesumc is the sum of products, to be used in computing cell-wise correlations)
  if (Nobs > 1) {
    edgesum1 = apply(adja1, c(1,2), sum); edgesum2 = apply(adja2, c(1,2), sum); 
    edgesumc = apply(adja1*adja2, c(1,2), sum)
  } else if (Nobs == 1) {
    if (length(dim(adja1)) == 2) { edgesum1 = adja1 } else { edgesum1 = adja1[,,1] }
    if (length(dim(adja2)) == 2) { edgesum2 = adja2 } else { edgesum2 = adja2[,,1] }
    edgesumc = edgesum1 * edgesum2
  }
  
  obs1_count_list = list(); obs2_count_list = list(); obsc_count_list = list(); obsp_count_list = list()
  cell_sizes_list = list()
  
  ## Compute cell-wise observed counts, depending on the fixed model type. 
  for(j in 1:Nmodels) {
    if (fit_models$mode == "tree") {
      obs1_count_list[[j]] = sapply(fit_models$model_list[[j]]$ft_expand, function(x) { sum(edgesum1[x[[1]], x[[2]]])})
      obs2_count_list[[j]] = sapply(fit_models$model_list[[j]]$ft_expand, function(x) { sum(edgesum2[x[[1]], x[[2]]])})
      obsc_count_list[[j]] = obs1_count_list[[j]] + obs2_count_list[[j]]
      obsp_count_list[[j]] = sapply(fit_models$model_list[[j]]$ft_expand, function(x) { sum(edgesumc[x[[1]], x[[2]]])})
      cell_sizes_list[[j]] = fit_models$model_list[[j]]$ft_counts * Nobs
    } else if (fit_models$mode %in% c( "block", "blockmodel")) {
      COR = fit_models$model_list[[j]]$bm_correct
      obs1_count_list[[j]] = sapply(fit_models$model_list[[j]]$bm_expand, function(x) { sum(edgesum1[x[[1]], x[[2]]])}) / COR
      obs2_count_list[[j]] = sapply(fit_models$model_list[[j]]$bm_expand, function(x) { sum(edgesum2[x[[1]], x[[2]]])}) / COR
      obsc_count_list[[j]] = obs1_count_list[[j]] + obs2_count_list[[j]]
      obsp_count_list[[j]] = sapply(fit_models$model_list[[j]]$bm_expand, function(x) { sum(edgesumc[x[[1]], x[[2]]])}) / COR
      cell_sizes_list[[j]] = fit_models$model_list[[j]]$bm_counts * Nobs
    } else if (fit_models$mode == "random") {
      vedgesum1 = as.vector(edgesum1); vedgesum2 = as.vector(edgesum2); vedgesumc = as.vector(edgesumc)
      obs1_count_list[[j]] = sapply(fit_models$model_list[[j]]$rs_ids, function(x) { sum(vedgesum1[x])})
      obs2_count_list[[j]] = sapply(fit_models$model_list[[j]]$rs_ids, function(x) { sum(vedgesum2[x])})
      obsc_count_list[[j]] = obs1_count_list[[j]] + obs2_count_list[[j]]
      obsp_count_list[[j]] = sapply(fit_models$model_list[[j]]$rs_ids, function(x) { sum(vedgesumc[x])})
      cell_sizes_list[[j]] = fit_models$model_list[[j]]$rs_counts * Nobs
    }
  }
  
  ## Compute MLE estimates (Count / Potential Edges)
  mle_p1 = lapply(1:Nmodels, function(j) {obs1_count_list[[j]] / cell_sizes_list[[j]]})
  mle_p2 = lapply(1:Nmodels, function(j) {obs2_count_list[[j]] / cell_sizes_list[[j]]})
  mle_pc = lapply(1:Nmodels, function(j) {obsc_count_list[[j]] / 2 / cell_sizes_list[[j]]})
  
  ## Compute cell-correlation sample estimates
  mle_pxy = lapply(1:Nmodels, function(j) {obsp_count_list[[j]] / cell_sizes_list[[j]]})
  cell_corrs = compute_samp_corrs(mle_pxy, mle_p1, mle_p2)
  
  ########## slow_compute_cellwise_logliks ##########
  ## Compute cellwise log-likelihoods
  LL_null = lapply(1:Nmodels, function(i) { compute_loglik_fromPC(x = obsc_count_list[[i]], n = 2 * cell_sizes_list[[i]], p = mle_pc[[i]]) })
  LL_alt1 = lapply(1:Nmodels, function(i) { compute_loglik_fromPC(x = obs1_count_list[[i]], n = cell_sizes_list[[i]], p = mle_p1[[i]]) })
  LL_alt2 = lapply(1:Nmodels, function(i) { compute_loglik_fromPC(x = obs2_count_list[[i]], n = cell_sizes_list[[i]], p = mle_p2[[i]]) })
  
  ## Compute cellwise test statistic from cell-wise log likelihoods
  cellwise_TS = lapply(1:Nmodels, function(i) { -2 * (LL_null[[i]] - LL_alt1[[i]] - LL_alt2[[i]]) })
  
  ########## slow_compute_adjust_pvals ##########
  pval_list = list()
  for(j in seq_along(cell_sizes_list)) {
    pval_list[[j]] = matrix(-1, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
    
    df_adj_matrix = matrix(-1, nrow = length(pl$cc_adj), ncol = length(cell_sizes_list[[j]]))
    for(i in seq_along(pl$cc_adj)) { df_adj_matrix[i,] = compute_df_adjustment2(n = cell_sizes_list[[j]], cell_corr = cell_corrs[[j]], cc_adj = pl$cc_adj[i]) } 
    
    for(i in seq_along(pl$thres_ignore)) {
      to_keep = which(cell_sizes_list[[j]] >= pl$thres_ignore[i])
      
      csq = sum(cellwise_TS[[j]][to_keep])
      dfs = apply(df_adj_matrix[,to_keep, drop = FALSE], 1, sum)
      
      ## Compute p-value and return
      if (return_chisq) {
        pval_list[[j]][,i] = csq
      } else {
        pval_list[[j]][,i] = pchisq(csq, dfs, lower.tail = FALSE)  
      }
    }
  }
  return(pval_list)
}



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (slow_compute_pval)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (slow_compute_pval)
#' Runs test for a set of fixed trees. Returns matrix of p-values based on parameter settings. 
#' 
#' This function is linked with fast_compute_pval_v2 -- this function just breaks the faster version into sub-functions for ease of testing. 
#' 
#' @param adja1 Adjacency array 1
#' @param adja2 Adjacency array 2
#' @param Nobs temp
#' @param pl Parameter list
#' @param fit_models Fixed models list for testing
#' @param return_chisq temp
#' 
#' @return Matrix of p-values (based on parameter list)
#' 
#' @export
#' 
slow_compute_pval = function(adja1, adja2, Nobs, pl, fit_models, return_chisq = FALSE) {
  ## This function uses the new structures... 
  ## This is just the 'slow' version of the code, as it calls test-able sub functions. Any changes to the subfunctions *SHOULD* correspond to a change in fast_compute_pval(version2, for now?)
  
  ## Compute cell-wise MLEs
  mle_lists = slow_compute_general_cellwise_mles(adja1 = adja1, adja2 = adja2, fit_models = fit_models, Nobs = Nobs)
  
  mle_p1 = mle_lists$mle_p1; mle_p2 = mle_lists$mle_p2; mle_pc = mle_lists$mle_pc
  mle_pxy = mle_lists$mle_pxy; cell_corrs = mle_lists$cell_corrs
  obs1_count_list = mle_lists$obs1_count_list; obs2_count_list = mle_lists$obs2_count_list
  obsc_count_list = mle_lists$obsc_count_list; cell_sizes_list = mle_lists$cell_sizes_list
  
  ## Compute log-likelihoods
  ll_lists = slow_compute_cellwise_logliks(obs1_count_list = obs1_count_list, obs2_count_list = obs2_count_list, obsc_count_list = obsc_count_list, cell_sizes_list = cell_sizes_list, mle_p1 = mle_p1, mle_p2 = mle_p2, mle_pc = mle_pc)
  
  LL_null = ll_lists$LL_null; LL_alt1 = ll_lists$LL_alt1; LL_alt2 = ll_lists$LL_alt2; cellwise_TS = ll_lists$cellwise_TS
  
  ## Compute p-values after adjustment
  
  pval_list = slow_compute_adjust_pvals(pl = pl, cell_corrs = cell_corrs, cell_sizes_list = cell_sizes_list, cellwise_TS = cellwise_TS, return_chisq = return_chisq)
  
  return(pval_list)
}


  
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (slow_compute_general_cellwise_mles)
#' Compute cell-wise MLEs and return all of the components needed
#' 
#' @param adja1 Adjacency array 1
#' @param adja2 Adjacency array 2
#' @param Nobs temp
#' @param fit_models Fixed models list for testing
#' 
#' @return List of many results
#' 
#' @export
#' 
slow_compute_general_cellwise_mles = function(adja1, adja2, Nobs, fit_models) {
  
  Nmodels = fit_models$Nmodels
  
  ## Compute total count for each edge (and edgesumc is the sum of products, to be used in computing cell-wise correlations)
  if (Nobs > 1) {
    edgesum1 = apply(adja1, c(1,2), sum); edgesum2 = apply(adja2, c(1,2), sum); 
    edgesumc = apply(adja1*adja2, c(1,2), sum)
  } else if (Nobs == 1) {
    if (length(dim(adja1)) == 2) { edgesum1 = adja1 } else { edgesum1 = adja1[,,1] }
    if (length(dim(adja2)) == 2) { edgesum2 = adja2 } else { edgesum2 = adja2[,,1] }
    edgesumc = edgesum1 * edgesum2
  }
  
  
  obs1_count_list = list(); obs2_count_list = list(); obsc_count_list = list(); obsp_count_list = list()
  cell_sizes_list = list()
  
  ## Compute cell-wise observed counts, depending on the fixed model type. 
  for(j in 1:Nmodels) {
    if (fit_models$mode == "tree") {
      obs1_count_list[[j]] = sapply(fit_models$model_list[[j]]$ft_expand, function(x) { sum(edgesum1[x[[1]], x[[2]]])})
      obs2_count_list[[j]] = sapply(fit_models$model_list[[j]]$ft_expand, function(x) { sum(edgesum2[x[[1]], x[[2]]])})
      obsc_count_list[[j]] = obs1_count_list[[j]] + obs2_count_list[[j]]
      obsp_count_list[[j]] = sapply(fit_models$model_list[[j]]$ft_expand, function(x) { sum(edgesumc[x[[1]], x[[2]]])})
      cell_sizes_list[[j]] = fit_models$model_list[[j]]$ft_counts * Nobs
    } else if (fit_models$mode %in% c( "block", "blockmodel")) {
      COR = fit_models$model_list[[j]]$bm_correct
      obs1_count_list[[j]] = sapply(fit_models$model_list[[j]]$bm_expand, function(x) { sum(edgesum1[x[[1]], x[[2]]])}) / COR
      obs2_count_list[[j]] = sapply(fit_models$model_list[[j]]$bm_expand, function(x) { sum(edgesum2[x[[1]], x[[2]]])}) / COR
      obsc_count_list[[j]] = obs1_count_list[[j]] + obs2_count_list[[j]]
      obsp_count_list[[j]] = sapply(fit_models$model_list[[j]]$bm_expand, function(x) { sum(edgesumc[x[[1]], x[[2]]])}) / COR
      cell_sizes_list[[j]] = fit_models$model_list[[j]]$bm_counts * Nobs
    } else if (fit_models$mode == "random") {
      vedgesum1 = as.vector(edgesum1); vedgesum2 = as.vector(edgesum2); vedgesumc = as.vector(edgesumc)
      obs1_count_list[[j]] = sapply(fit_models$model_list[[j]]$rs_ids, function(x) { sum(vedgesum1[x])})
      obs2_count_list[[j]] = sapply(fit_models$model_list[[j]]$rs_ids, function(x) { sum(vedgesum2[x])})
      obsc_count_list[[j]] = obs1_count_list[[j]] + obs2_count_list[[j]]
      obsp_count_list[[j]] = sapply(fit_models$model_list[[j]]$rs_ids, function(x) { sum(vedgesumc[x])})
      cell_sizes_list[[j]] = fit_models$model_list[[j]]$rs_counts * Nobs
    }
  }
  
  ## Compute MLE estimates (Count / Potential Edges)
  mle_p1 = lapply(1:Nmodels, function(j) {obs1_count_list[[j]] / cell_sizes_list[[j]]})
  mle_p2 = lapply(1:Nmodels, function(j) {obs2_count_list[[j]] / cell_sizes_list[[j]]})
  mle_pc = lapply(1:Nmodels, function(j) {obsc_count_list[[j]] / 2 / cell_sizes_list[[j]]})
  
  ## Compute cell-correlation sample estimates
  mle_pxy = lapply(1:Nmodels, function(j) {obsp_count_list[[j]] / cell_sizes_list[[j]]})
  cell_corrs = compute_samp_corrs(mle_pxy, mle_p1, mle_p2)
  
  return(list(obs1_count_list = obs1_count_list, obs2_count_list = obs2_count_list, obsc_count_list = obsc_count_list, cell_sizes_list = cell_sizes_list, mle_p1 = mle_p1, mle_p2 = mle_p2, mle_pc = mle_pc, mle_pxy = mle_pxy, cell_corrs = cell_corrs))
}



#' Compute cellwise log-likelihoods
#' 
#' @param obs1_count_list List of cell-counts for adjacency matrix 1
#' @param obs2_count_list List of cell-counts for adjacency matrix 2
#' @param obsc_count_list List of combined cell-counts
#' @param cell_sizes_list List of cell-sizes
#' @param mle_p1 List of MLEs for adjacency matrix 1
#' @param mle_p2 List of MLEs for adjacency matrix 2
#' @param mle_pc List of combined MLEs (null hypothesis MLE)
#' 
#' @return List of log-likelihoods
#' 
#' @export
#' 
slow_compute_cellwise_logliks = function(obs1_count_list, obs2_count_list, obsc_count_list, cell_sizes_list, mle_p1, mle_p2, mle_pc) {
  Nmodels = length(mle_p1)
  
  ## Compute cellwise log-likelihoods
  LL_null = lapply(1:Nmodels, function(i) { compute_loglik_fromPC(x = obsc_count_list[[i]], n = 2 * cell_sizes_list[[i]], p = mle_pc[[i]]) })
  LL_alt1 = lapply(1:Nmodels, function(i) { compute_loglik_fromPC(x = obs1_count_list[[i]], n = cell_sizes_list[[i]], p = mle_p1[[i]]) })
  LL_alt2 = lapply(1:Nmodels, function(i) { compute_loglik_fromPC(x = obs2_count_list[[i]], n = cell_sizes_list[[i]], p = mle_p2[[i]]) })
  
  ## Compute cellwise test statistic from cell-wise log likelihoods
  cellwise_TS = lapply(1:Nmodels, function(i) { -2 * (LL_null[[i]] - LL_alt1[[i]] - LL_alt2[[i]]) })
  
  return(list(LL_null = LL_null, LL_alt1 = LL_alt1, LL_alt2 = LL_alt2, cellwise_TS = cellwise_TS))
}



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (slow_compute_adjust_pvals)
#' Compute adjusted p-values
#' 
#' @param pl Parameter list
#' @param cell_corrs List of cell-wise correlation estimates
#' @param cell_sizes_list List of cell-sizes
#' @param cellwise_TS List of cell-wise test statistics (-2 loglikelihoodratio)
#' @param return_chisq temp
#' 
#' @return Matrix of p-values
#' 
#' @export
#' 
slow_compute_adjust_pvals = function(pl, cell_corrs, cell_sizes_list, cellwise_TS, return_chisq = FALSE) {
  pval_list = list()
  for(j in seq_along(cell_sizes_list)) {
    pval_list[[j]] = matrix(-1, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
    
    df_adj_matrix = matrix(-1, nrow = length(pl$cc_adj), ncol = length(cell_sizes_list[[j]]))
    for(i in seq_along(pl$cc_adj)) { df_adj_matrix[i,] = compute_df_adjustment2(n = cell_sizes_list[[j]], cell_corr = cell_corrs[[j]], cc_adj = pl$cc_adj[i]) } 
    
    for(i in seq_along(pl$thres_ignore)) {
      to_keep = which(cell_sizes_list[[j]] >= pl$thres_ignore[i])
      
      csq = sum(cellwise_TS[[j]][to_keep])
      dfs = apply(df_adj_matrix[,to_keep, drop = FALSE], 1, sum)
      
      ## Compute p-value and return
      if (return_chisq) {
        pval_list[[j]][,i] = csq
      } else {
        pval_list[[j]][,i] = pchisq(csq, dfs, lower.tail = FALSE)  
      }
    }
  }
  return(pval_list)
}

