

#' Do simulations for the hypothesis testing, using the new code. 
#' 
#' @param gen_NetMPair NetworkModelPair object
#' @param fit_NetSList NetworkStructList object (default = NULL)
#' @param fitm_params set_model_params() -- use this to specify the type of fitting model
#' @param adjm_list Can input a list of adjacency matrices. If done, then the generating model pair is ignored. 
#' @param Nobs Number network observations (default = 1)
#' @param Nsim Number of simulations to do
#' @param param_list Parameter list for testing procedure
#' @param pval_adj_fx List of p-value adjustment functions
#' @param verbose level: 0 -- no output. > 3 -- full output
#' 
#' @return List of results
#' 
#' @export
#' 
sim_hyptest = function(gen_NetMPair, fit_NetSList = NULL, fitm_params = set_model_param(), adjm_list = NULL, 
                       Nobs = 1, Nsim = 100, param_list = set_sim_param(), 
                       pval_adj_fx = list(mult_bonferroni, mult_pearson, mult_highcrit), verbose = 0) {
  
  
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
  if (is.null(fit_NetSList)) { fit_NetSList = NetworkStructList(Nmodels = Nfit, model_param = fitm_params) }
  if (Nfit > length(fit_NetSList@models)) { stop("Not enough fitted models generated") }
  partial_indices = lapply(param_list$n_models, function(x) { 1:x })
  
  ## Setup storage
  result_list = list(); for(f in seq_along(pval_adj_fx)) { result_list[[f]] = list() } 
  pval_reslist = list()  
  
  ## Generate two lists of adjacency arrays
  if (is.null(adjm_list)) { adjm_list = sampleNetwork(gen_NetMPair, Nsim = Nsim) }
  
  ## Do simulations
  for(j in 1:Nsim) {  
    if (verbose > 3) { cat("."); if (j %% floor(Nsim/10) == 0) { print(paste("Simulation number:", j)) } }
    pval_results = computePval(fit_NetSList, adja1 = adjm_list[[1]][[j]], adja2 = adjm_list[[2]][[j]], pl = param_list, Nobs = 1)
    pval_reslist[[j]] = abind(pval_results, along = 3)
    
    for(f in seq_along(pval_adj_fx)) {
      result_list[[f]][[j]] = setup_array(param_list)
      for (k in seq_along(param_list$n_models)) {
        for (a in seq_along(param_list$alphas)) {
          result_list[[f]][[j]][,,a,k] = apply(pval_reslist[[j]][,,partial_indices[[k]], drop = FALSE], c(1,2),pval_adj_fx[[f]])
        }
      }
    }
  }
  return(result_list)
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sim_critvals)
#' Simulates critical values for each type of multiple testing adjustment
#' 
#' @param NetMPair pair of network models
#' @param Nsim number of simulations
#' @param Nobs number of network observations
#' @param fit_NetSList temp
#' @param fit_models_params parameter settings
#' @param param_list list of parameters for testing
#' @param pval_adj_fx functions for p-value adjustment
#' @param verbose verbose level: if > 2, output stuff. 
#' 
#' @return list of arrays of critical values
#' 
#' @export
#' 
sim_critvals = function(NetMPair, Nsim = 500, Nobs = 1, fit_NetSList = NULL, fit_models_params = set_model_param(), param_list = set_sim_param(), pval_adj_fx = list(mult_pearson, mult_highcrit), verbose = 0) { 
  if (is.null(fit_NetSList)) {
    fit_NetSList = NetworkStructList(model_param = fit_models_params, Nmodels = max(param_list$n_models))
  }
  
  if (verbose > 2) { print("** Calling sim_hyptest **")}
  sim_vals = sim_hyptest(
    gen_NetMPair = NetMPair, fit_NetSList = fit_NetSList, 
    Nobs = Nobs, Nsim = Nsim, param_list = param_list, 
    pval_adj_fx = pval_adj_fx, verbose = verbose)
  
  if (verbose > 2) { print("** Aggregating results **")}
  res_list = list()
  cases = expand.grid(param_list)
  for(j in seq_along(pval_adj_fx)) {
    res_list[[j]] = setup_array(param_list)
    for(k in seq_len(nrow(cases))) {
      res_list[[j]][cases$cc_adj[k] == param_list$cc_adj, 
                    cases$thres_ignore[k] == param_list$thres_ignore,
                    cases$alphas[k] == param_list$alphas, 
                    cases$n_models[k] == param_list$n_models] = quantile(
                      sim_subsetresults(sim_vals[j], param_list, 
                                        cases$cc_adj[k], cases$thres_ignore[k], cases$alphas[k], cases$n_models[k]), 
                      probs = 1 - cases$alphas[k])
    }
  }
  return(res_list)
}


#' Create a named array matching a list
#' 
#' @param pl parameter list for which array is desired
#' 
#' @return NA-filled named array
#' 
#' @export
#' 
setup_array = function(pl) {
  ## creates an array that names things properly
  if (!is.list(pl)) { stop("Input is not a list") }
  
  lengths = sapply(pl, length)
  namelist = list()
  for (j in seq_along(pl)) {
    namelist[[names(pl)[j]]] = pl[[j]]
  }
  
  res = array(NA, dim = lengths, dimnames = namelist)
  return(res)
}


#' Extracts results matching a subset of parameters
#' 
#' @param rl result list -- should be the output of sim_hyptest (or something similar)
#' @param pl parameter list -- as given by set_sim_param()
#' @param cc_adj Amount of SE's away from the correlation estimate used (how conservative? 0 means no adjustment (and requires large sample for guarantees; larger values give a conservative p-value))
#' @param thres_ignore Ignore edge groups with fewer than this many edges
#' @param alphas Size of test
#' @param n_models Number of edge partitions to use for testing
#' 
#' @return a matrix of results of a certain type
#' 
#' @export
#' 
sim_subsetresults = function(rl, pl, cc_adj, thres_ignore, alphas, n_models) { 
  ## is similar to extract_result_list (but not sure to remove extract_result_list yet)
  ## extracts results of a certain type
  
  i1 = which(pl$cc_adj == cc_adj)
  i2 = which(pl$thres_ignore == thres_ignore)
  i3 = which(pl$alphas == alphas)
  i4 = which(pl$n_models == n_models)
  
  res = matrix(-1, nrow = length(rl[[1]]), ncol = length(rl))
  for(j in seq_along(rl)) {
    res[,j] = sapply(rl[[j]], function(x) { x[i1, i2, i3, i4]} )
  }
  return(res)
}


sim_power_rpart = function(GL, NL, FL, Nsim = 500, Nsim_crit = 100, Nobs = 1, verbose = 0, pl, pval_adj_fx = set_pval_fx()) {
  ## then, will also need a function that helps generate the lists of models? maybe? 
  ## for simulations -- given a list of generating models, a list of null-hypothesis models (to generate crit values from), and a corresspnding list of fitting models: 
  
  ## for each pair, simulate the critical values
  ## then -- do simulation
  
  ## input list of generating models is GL, null-hyp -> NL, fitting models -> FL
  if (FALSE) {
    GL = list(NetworkModelPair(m1 = NetworkModelSBM(Nnodes = 30), is_null = TRUE), 
              NetworkModelPair(m1 = NetworkModelSBM(Nnodes = 30), is_null = TRUE))
    NL = GL
    FL = list(set_model_param(Nnodes = 30), set_model_param(Nnodes = 30))
    Nsim = 200
    Nsim_crit = 200
    Nobs = 1
    verbose = 5
    pl = set_sim_param()
    pval_adj_fx = set_pval_fx()
  }
  
  
  #pl = param list for simulations
  if (length(GL) != length(NL) || length(NL) != length(FL)) { stop("The model lists are not compatible") }  
  power_list = list()
  cases = expand.grid(pl)
  
  ## Do simulations: 
  for (S in seq_along(GL)) {
    
    ## If FL contains networkstructlists, then use those for everything. otherwise, generate a new one each time...
    if (is(FL[[S]], "NetworkStructList")) { 
      if (verbose > 0) { cat("===== Simulation Number: ", S, " --- (Simulating Critical Values) =====\n", sep = "") }
      critvals = sim_critvals(NetMPair = NL[[S]], Nsim = Nsim_crit, Nobs = Nobs, fit_NetSList = FL[[S]], param_list = pl, pval_adj_fx = pval_adj_fx, verbose = verbose)
      
      if (verbose > 0) { cat("===== Simulation Number: ", S, " --- (Simulating Power) =====\n", sep = "") }
      sim_res = sim_hyptest(gen_NetMPair = GL[[S]], fit_NetSList = FL[[S]], Nobs = Nobs, Nsim = Nsim, param_list = pl, pval_adj_fx = pval_adj_fx, verbose = verbose)
    } else {
      if (verbose > 0) { cat("===== Simulation Number: ", S, " --- (Simulating Critical Values) =====\n", sep = "") }
      critvals = sim_critvals(NetMPair = NL[[S]], Nsim = Nsim_crit, Nobs = Nobs, fit_models_params = FL[[S]], param_list = pl, pval_adj_fx = pval_adj_fx, verbose = verbose)
      
      if (verbose > 0) { cat("===== Simulation Number: ", S, " --- (Simulating Power) =====\n", sep = "") }
      sim_res = sim_hyptest(gen_NetMPair = GL[[S]], fitm_params = FL[[S]], Nobs = Nobs, Nsim = Nsim, param_list = pl, pval_adj_fx = pval_adj_fx, verbose = verbose)
    }
    
    ## Use critical values to reject if necessary
    if (verbose > 0) { cat("===== Simulation Number: ", S, " --- (Computing Power) =====\n", sep = "") }
    power_list[[S]] = list()
    
    
    for(j in seq_along(pval_adj_fx)) {
      power_list[[S]][[j]] = setup_array(pl)
      for(k in seq_len(nrow(cases))) {
        i1 = which(cases$cc_adj[k] == pl$cc_adj); i2 = cases$thres_ignore[k] == pl$thres_ignore
        i3 = cases$alphas[k] == pl$alphas; i4 = cases$n_models[k] == pl$n_models
        
        test_stats = sim_subsetresults(sim_res, pl, cases$cc_adj[k], cases$thres_ignore[k], cases$alphas[k], cases$n_models[k])
        
        if (names(pval_adj_fx)[j] == "mult_bonferroni") {
          power_list[[S]][[j]][i1, i2, i3, i4] = mean(test_stats[,j] < 0.05 / cases$n_models[k])
        } else {
          power_list[[S]][[j]][i1, i2, i3, i4] = mean(test_stats[,j] > critvals[[j]][i1, i2, i3, i4])
        }
      }
    }    
  }
  
  
  return(power_list)
}



