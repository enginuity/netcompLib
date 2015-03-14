
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sim_hyptest)
#' <What does this function do>
#' 
#' @return temp
#' 
#' @export
#' 
sim_hyptest = function(gen_NetMPair, fit_NetSList, Nobs = 1, Nsim = 100, param_list, pval_adj_fx) {
  # Does simulation for a given pair of NetM, NetS
  Nfit = max(param_list$n_models)
  dir.create("netcomp_logfile", showWarnings = FALSE)
  logfile = paste("netcomp_logfile/SIMNULL_", gsub("[^0-9]", "", Sys.time()), "_log.txt", sep = "")  
  if (verbose) {cat("Simulation log \n", file = logfile, append = TRUE)}
  
  if (Nfit > length(fit_NetSList@models)) { stop("Not enough fitted models generated") }
  partial_indices = lapply(param_list$n_models, function(x) { 1:x })
  
  ## Setup storage
  adjm_list = list() ## the implementation of this changes. This is now two separate lists of adjacency matrices. 
  result_list = list(); for(f in seq_along(pval_adj_fx)) { result_list[[f]] = list() } 
  pval_reslist = list()  
  
  
  

  ## Do simulations
  for(j in 1:Nsim) {
    
    if (j %% floor(Nsim/10) == 0) { if (verbose) {print(paste("Simulation number:", j))} }
    
    
    ## fix starting here -- 
    
    if (is.null(adjmats)) {
      adjm_list[[j]] = sample_network_pair(gen_model = gen_models, Nobs = Nobs)
    } else {
      adjm_list[[j]] = adjmats[[j]]
    }
    
    pval_results = fast_compute_pval_v2(adja1 = adjm_list[[j]][[1]], adja2 = adjm_list[[j]][[2]], Nobs = Nobs, 
                                        pl = param_list, fit_models = fit_models, return_chisq = return_chisq)
    
    pval_results = abind(pval_results, along = 3)
    pval_reslist[[j]] = pval_results
    
    for(f in seq_along(pval_adj_fx)) {
      result_list[[f]][[j]] = array(0, dim = c(length(param_list$cc_adj), length(param_list$thres_ignore), length(param_list$n_models), length(param_list$alphas)))
      for (k in seq_along(param_list$n_models)) {
        for (a in seq_along(param_list$alphas)) {
          result_list[[f]][[j]][,,k,a] = apply(pval_results[,,partial_indices[[k]], drop = FALSE], c(1,2),pval_adj_fx[[f]])
        }
      }
    }
  }
  
  ## Output full or partial depending on 'complete'
  if (return_chisq) { return(pval_reslist) }
  if (complete) { 
    full_output = list(gen_mode = gen_mode, gen_models = gen_models, 
                       fit_mode = fit_mode, fit_models = fit_models, 
                       adjm_list = adjm_list, param_list = param_list,
                       result_list = result_list, pval_reslist = pval_reslist)
    return(full_output) 
  } else { 
    return(result_list) 
  }
}
