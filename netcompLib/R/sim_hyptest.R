
#' Do simulations for the hypothesis testing, using the new code. 
#' 
#' @param gen_NetMPair NetworkModelPair object
#' @param fit_NetSList NetworkStructList object
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
sim_hyptest = function(gen_NetMPair, fit_NetSList, Nobs = 1, Nsim = 100, param_list, pval_adj_fx = list(mult_bonferroni, mult_pearson, mult_highcrit), verbose = TRUE) {
  if(FALSE) {
    library(netcompLib)
    library(abind)
    load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
    
    gen_NetMPair = NetworkModelPair(NetworkModelSBM(Nnodes = 30), NetworkModelSBM(Nnodes = 30), is_null = FALSE)
    fit_NetSList = NetworkStructList(Nnodes = 30, Nmodels = 100, type = "block")
    Nobs = 1; Nsim = 100; verbose = TRUE;
    param_list = set_sim_param(cc_adj = c(0,1,2), thres_ignore = c(5, 10), n_models = c(1, 25, 50, 100))
    pval_adj_fx = list(mult_bonferroni, mult_pearson, mult_highcrit)
  }
  # Does simulation for a given pair of NetM, NetS
  Nfit = max(param_list$n_models)
  dir.create("netcomp_logfile", showWarnings = FALSE)
  logfile = paste("netcomp_logfile/SIMNULL_", gsub("[^0-9]", "", Sys.time()), "_log.txt", sep = "")  
  if (verbose) {cat("Simulation log \n", file = logfile, append = TRUE)}
  
  if (Nfit > length(fit_NetSList@models)) { stop("Not enough fitted models generated") }
  partial_indices = lapply(param_list$n_models, function(x) { 1:x })
  
  ## Setup storage
  result_list = list(); for(f in seq_along(pval_adj_fx)) { result_list[[f]] = list() } 
  pval_reslist = list()  
  
  ## Generate adjacency matrices
  adjm_list = sampleNetwork(gen_NetMPair, Nsim = Nsim) # This is two lists of adjacency arrays. 
  # TODO: allow inputting an adjm_list instead of having to generate a random one
  
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

