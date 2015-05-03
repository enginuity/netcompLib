
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sim_hyptest)
#' Do simulations for the hypothesis testing, using the new code. 
#' 
#' @param gen_NetMPair NetworkModelPair object
#' @param fit_NetSList NetworkStructList object
#' @param adjm_list temp
#' @param Nobs Number network observations (default = 1)
#' @param Nsim Number of simulations to do
#' @param param_list Parameter list for testing procedure
#' @param pval_adj_fx List of p-value adjustment functions
#' @param verbose T/F: output?
#' 
#' @return List of results
#' 
#' @export
#' 
sim_hyptest = function(gen_NetMPair, fit_NetSList, adjm_list = NULL, 
                       Nobs = 1, Nsim = 100, param_list = set_sim_param(), 
                       pval_adj_fx = list(mult_bonferroni, mult_pearson, mult_highcrit), verbose = TRUE) {
  ## TODO: add to documentation
  ## adjm_list -- if null, ignore. if non-null, we can ignore gen_NetMPair (since the adjacency matrix list is already passed in)
  
  ## Test case
  if(FALSE) {
    library(netcompLib); library(abind)
    load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
    
    gen_NetMPair = NetworkModelPair(NetworkModelSBM(Nnodes = 30), NetworkModelSBM(Nnodes = 30), is_null = FALSE)
    fit_NetSList = NetworkStructList(Nnodes = 30, Nmodels = 100, type = "block")
    Nobs = 1; Nsim = 100; verbose = TRUE;
    param_list = set_sim_param(cc_adj = c(0,1,2), thres_ignore = c(5, 10), n_models = c(1, 25, 50, 100))
  }
  
  # Check that enough models were generated
  Nfit = max(param_list$n_models)
  if (Nfit > length(fit_NetSList@models)) { stop("Not enough fitted models generated") }
  partial_indices = lapply(param_list$n_models, function(x) { 1:x })
  
  ## Setup storage
  result_list = list(); for(f in seq_along(pval_adj_fx)) { result_list[[f]] = list() } 
  pval_reslist = list()  
  
  ## Generate two lists of adjacency arrays
  if (is.null(adjm_list)) { adjm_list = sampleNetwork(gen_NetMPair, Nsim = Nsim) }
  
  ## Do simulations
  for(j in 1:Nsim) {  
    if (verbose) { cat("."); if (j %% floor(Nsim/10) == 0) { print(paste("Simulation number:", j)) } }
    pval_results = computePval(fit_NetSList, adja1 = adjm_list[[1]][[j]], adja2 = adjm_list[[2]][[j]], pl = param_list, Nobs = 1)
    pval_reslist[[j]] = abind(pval_results, along = 3)
    
    for(f in seq_along(pval_adj_fx)) {
      result_list[[f]][[j]] = array(0, dim = c(length(param_list$cc_adj), length(param_list$thres_ignore), length(param_list$n_models), length(param_list$alphas)), dimnames = list(cc_adj = param_list$cc_adj, thres_ignore = param_list$thres_ignore, n_models = param_list$n_models, alphas = param_list$alphas))
      for (k in seq_along(param_list$n_models)) {
        for (a in seq_along(param_list$alphas)) {
          result_list[[f]][[j]][,,k,a] = apply(pval_reslist[[j]][,,partial_indices[[k]], drop = FALSE], c(1,2),pval_adj_fx[[f]])
        }
      }
    }
  }
  return(result_list)
}



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sim_critvals)
#' <What does this function do>
#' 
#' @param NetMPair temp
#' @param Nsim temp
#' @param Nobs temp
#' @param fit_models_type temp
#' @param fit_models_params temp
#' @param param_list temp
#' @param pval_adj_fx temp
#' 
#' @return temp
#' 
#' @export
#' 
sim_critvals = function(NetMPair, Nsim = 500, Nobs = 1, fit_models_type, fit_models_params = set_model_param(), param_list = set_sim_param(), pval_adj_fx = list(mult_pearson, mult_highcrit)) { 
  ## netMPair = generating model pair, can be null hypothesis
  # fit_models_type = 'block', 'tree', 'random' -- should add this parameter to set_model_param?
  
  sim_vals = sim_hyptest(
    gen_NetMPair = NetMPair, fit_NetSList = NetworkStructList(
      Nnodes = getNnodes(NetMPair), Nmodels = max(param_list$n_models), 
      type = fit_models_type, model_param = fit_models_params), 
    Nobs = Nobs, Nsim = Nsim, param_list = param_list, 
    pval_adj_fx = pval_adj_fx, verbose = verbose)
  
  res_list = list()
  for(j in seq_along(pval_adj_fx)) {
    res_list[[j]] = sapply(seq_along(param_list$n_models), function(y) { quantile(x = sapply(sim_vals[[2]], function(x) {x[1,1,y,1]}), probs = 0.95) })
  }
  return(res_list)
}

