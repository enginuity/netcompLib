


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (network_comparison)
#' <<BasicInfo>> 
#' 
#' @param first.network temp
#' @param second.network temp
#' @param verbose temp
#' @param iterations temp
#' @param new.replicates temp
#' @param run.id temp
#' 
#' @return temp
#' 
#' @export
#' 
network_comparison = function (first.network,
                               second.network,
                               verbose=TRUE,
                               iterations=NULL,
                               new.replicates=100,
                               run.id=NULL) {
  # Written by Andrew Thomas
  
  
  #blocks <- 3;  netdraw <- simple_block_model(block.matrix={out <- array(0.05, rep(blocks,2)); diag(out) <- 0.7; out}); first.network <- netdraw[,,1:5]; second.network <- netdraw[,,5+1:5]; verbose=TRUE;  iterations=NULL; new.replicates=100; run.id=NULL
  #source("allcode-general.R")
  
  if (is.null(run.id)) run.id <- as.character(round(1000000*runif(1)))
  
  ##### Make them both 3D arrays.
  if (is.na(dim(first.network)[3])) first.network <- array(c(first.network), c(dim(first.network)[1:2], 1))
  if (is.na(dim(second.network)[3])) second.network <- array(c(second.network), c(dim(second.network)[1:2], 1))
  dim.net.1 <- dim(first.network); dim.net.2 <- dim(second.network);
  
  if (any(dim.net.1[1:2] != dim.net.2[1:2])) stop ("Network node sizes do not match.")
  if (dim.net.1[1] != dim.net.1[2]) stop ("Input networks not square.")
  message(paste("Nodes: ",dim.net.1[1],", Total ties: ",dim.net.1[3]," ",dim.net.1[3],sep=""))
  
  if (is.null(iterations)) iterations <- round(dim.net.1[1]^2)
  
  combined.network <- array(c(c(first.network), c(second.network)), c(dim.net.1[1:2], dim.net.1[3]+dim.net.2[3]))
  
  if (verbose) message ("Begin solution, network 1")
  net.sol.1 <- cmn_mcmc (first.network, iterations=iterations, verbose=verbose)
  #|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
  if (verbose) message ("Begin solution, network 2")
  net.sol.2 <- cmn_mcmc (second.network, iterations=iterations, verbose=verbose)
  #|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
  
  true.dist <- tree_distance(net.sol.1$choice.tree.ret, net.sol.2$choice.tree.ret)
  
  
  if (verbose) message ("Begin solution, combined network")
  net.sol.comb <- cmn_mcmc (combined.network, iterations=iterations, verbose=verbose)
  #|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
  
  #simulation time!   new.replicates=15
  library(doMC)
  registerDoMC(15)
  
  done.sims <- foreach (kk=1:new.replicates) %dopar% {   #replicate (new.replicates, {
    
    sim.nets <- cmn_network (net.sol.comb$choice.tree.ret, dim(combined.network)[3])
    
    if (verbose) message ("Begin solution, sim network 1")
    net.sol.1 <- cmn_mcmc (sim.nets[,,1:dim(first.network)[3]], iterations=iterations, verbose=verbose)
    #|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
    if (verbose) message ("Begin solution, sim network 2")
    net.sol.2 <- cmn_mcmc (sim.nets[,,dim(first.network)[3]+1:dim(second.network)[3]],
                           #|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
                           iterations=iterations, verbose=verbose)
    
    return (list(net.sol.1=net.sol.1, net.sol.2=net.sol.2))
  }
  distance.distribution <- sapply(done.sims, function(item) tree_distance(item$net.sol.1$choice.tree.ret, item$net.sol.2$choice.tree.ret))
  
  save(net.sol.1, net.sol.2, net.sol.comb, done.sims, file=paste(run.id, "-netcomp-sims.RData", sep=""))
  
  
  
  
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (network_comparison)
#' <<BasicInfo>> 
#' 
#' @param first.network temp
#' @param second.network temp
#' @param verbose temp
#' @param iterations temp
#' @param new.replicates temp
#' @param run.id temp
#' @param prob_est temp
#' @param prob_smooth_param temp
#' @param smooth_root temp
#' 
#' @return temp
#' 
#' @export
#' 
network_comparison = function(first.network,
                                second.network,
                                verbose=TRUE,
                                iterations=NULL,
                                new.replicates=100,
                                run.id=NULL,
                                prob_est = c("default", "convex", "convex_scaled"),
                                prob_smooth_param = 0.5,
                                smooth_root = TRUE) {

  #blocks <- 3;  netdraw <- simple_block_model(block.matrix={out <- array(0.05, rep(blocks,2)); diag(out) <- 0.7; out}); first.network <- netdraw[,,1:5]; second.network <- netdraw[,,5+1:5]; verbose=TRUE;  iterations=NULL; new.replicates=100; run.id=NULL
  #source("allcode-general.R")
  
  if (is.null(run.id)) run.id <- as.character(round(1000000*runif(1)))
  
  ##### Make them both 3D arrays.
  if (is.na(dim(first.network)[3])) first.network <- array(c(first.network), c(dim(first.network)[1:2], 1))
  if (is.na(dim(second.network)[3])) second.network <- array(c(second.network), c(dim(second.network)[1:2], 1))
  dim.net.1 <- dim(first.network); dim.net.2 <- dim(second.network);

  if (any(dim.net.1[1:2] != dim.net.2[1:2])) stop ("Network node sizes do not match.")
  if (dim.net.1[1] != dim.net.1[2]) stop ("Input networks not square.")
  message(paste("Nodes: ",dim.net.1[1],", Total ties: ",dim.net.1[3]," ",dim.net.1[3],sep=""))

  if (is.null(iterations)) iterations <- round(dim.net.1[1]^2)
  
  combined.network <- array(c(c(first.network), c(second.network)), c(dim.net.1[1:2], dim.net.1[3]+dim.net.2[3]))

  if (verbose) message ("Begin solution, network 1")
  net.sol.1 <- cmn_mcmc (first.network, iterations=iterations, verbose=verbose, prob_est = prob_est, prob_smooth_param = prob_smooth_param, smooth_root = smooth_root)
#|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
  if (verbose) message ("Begin solution, network 2")
  net.sol.2 <- cmn_mcmc (second.network, iterations=iterations, verbose=verbose, prob_est = prob_est, prob_smooth_param = prob_smooth_param, smooth_root = smooth_root)
#|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--

  true.dist <- tree_distance(net.sol.1$choice.tree.ret, net.sol.2$choice.tree.ret)

  
  if (verbose) message ("Begin solution, combined network")
  net.sol.comb <- cmn_mcmc (combined.network, iterations=iterations, verbose=verbose, prob_est = prob_est, prob_smooth_param = prob_smooth_param, smooth_root = smooth_root)
#|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--

  #simulation time!   new.replicates=15
  library(doMC)
  
  registerDoMC(15)
  
  done.sims <- foreach (kk=1:new.replicates) %dopar% {   #replicate (new.replicates, {

    sim.nets <- cmn_network (net.sol.comb$choice.tree.ret, dim(combined.network)[3])
    
    if (verbose) message ("Begin solution, sim network 1")
    net.sol.1 <- cmn_mcmc (sim.nets[,,1:dim(first.network)[3]], iterations=iterations, verbose=verbose, prob_est = prob_est, prob_smooth_param = prob_smooth_param, smooth_root = smooth_root)
#|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
    if (verbose) message ("Begin solution, sim network 2")
    net.sol.2 <- cmn_mcmc (sim.nets[,,dim(first.network)[3]+1:dim(second.network)[3]],
#|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
                           iterations=iterations, verbose=verbose, prob_est = prob_est, prob_smooth_param = prob_smooth_param, smooth_root = smooth_root)

    return (list(net.sol.1=net.sol.1, net.sol.2=net.sol.2))
  }
  distance.distribution <- sapply(done.sims, function(item) tree_distance(item$net.sol.1$choice.tree.ret, item$net.sol.2$choice.tree.ret))

  save(net.sol.1, net.sol.2, net.sol.comb, done.sims, file=paste(run.id, "-netcomp-sims.RData", sep=""))

  return(list(true.dist, distance.distribution))

}












### Version for using C Code. 
## TODO: [Refactor] the version of the code using C can be combined with the raw version. 

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (network_comparison)
#' <<BasicInfo>> 
#' 
#' @param first.network temp
#' @param second.network temp
#' @param verbose temp
#' @param iterations temp
#' @param new.replicates temp
#' @param run.id temp
#' @param prob_est temp
#' @param prob_smooth_param temp
#' @param smooth_root temp
#' 
#' @return temp
#' 
#' @export
#' 
network_comparison_C = function(first.network,
                                second.network,
                                verbose=TRUE,
                                iterations=NULL,
                                new.replicates=100,
                                run.id=NULL,
                                prob_est = c("default", "convex", "convex_scaled"),
                                prob_smooth_param = 0.5,
                                smooth_root = TRUE) {

  #blocks <- 3;  netdraw <- simple_block_model(block.matrix={out <- array(0.05, rep(blocks,2)); diag(out) <- 0.7; out}); first.network <- netdraw[,,1:5]; second.network <- netdraw[,,5+1:5]; verbose=TRUE;  iterations=NULL; new.replicates=100; run.id=NULL
  #source("allcode-general.R")
  
  if (is.null(run.id)) run.id <- as.character(round(1000000*runif(1)))
  
  ##### Make them both 3D arrays.
  if (is.na(dim(first.network)[3])) first.network <- array(c(first.network), c(dim(first.network)[1:2], 1))
  if (is.na(dim(second.network)[3])) second.network <- array(c(second.network), c(dim(second.network)[1:2], 1))
  dim.net.1 <- dim(first.network); dim.net.2 <- dim(second.network);

  if (any(dim.net.1[1:2] != dim.net.2[1:2])) stop ("Network node sizes do not match.")
  if (dim.net.1[1] != dim.net.1[2]) stop ("Input networks not square.")
  message(paste("Nodes: ",dim.net.1[1],", Total ties: ",dim.net.1[3]," ",dim.net.1[3],sep=""))

  if (is.null(iterations)) iterations <- round(dim.net.1[1]^2)
  
  combined.network <- array(c(c(first.network), c(second.network)), c(dim.net.1[1:2], dim.net.1[3]+dim.net.2[3]))

  if (verbose) message ("Begin solution, network 1")
  ST = starter_tree(dim.net.1[1])
  iserror = TRUE
  while(iserror) {
    net_sol_1 = try(run_C_mcmc(ST$parents, adj_mat = first.network,
      num_iters = iterations, verbose = verbose, convex_wt = prob_smooth_param))
    if (class(net_sol_1) != "try-error") { iserror = FALSE }
  }

  if (verbose) message ("Begin solution, network 2")
  iserror = TRUE
  while(iserror) {
    net_sol_2 = try(run_C_mcmc(ST$parents, adj_mat = second.network,
      num_iters = iterations, verbose = verbose, convex_wt = prob_smooth_param))
    if (class(net_sol_2) != "try-error") { iserror = FALSE }
  }

  true.dist = tree_distance(net_sol_1$bestTree, net_sol_2$bestTree)

  
  if (verbose) message ("Begin solution, combined network")
  iserror = TRUE
  while(iserror) {
    net_sol_comb = try(run_C_mcmc(ST$parents, adj_mat = combined.network,
      num_iters = iterations, verbose = verbose, convex_wt = prob_smooth_param))
    if (class(net_sol_comb) != "try-error") { iserror = FALSE }
  }  
  
  dist_distrib = rep(-1, times = new.replicates)
  for(iter in 1:new.replicates) {
    if (iter %% 10 == 0) {cat(".")}
    sim.nets = cmn_network(net_sol_comb$bestTree, dim(combined.network)[3])
    
    ST = starter_tree(dim.net.1[1])
    iserror = TRUE
    while(iserror) {
      net_sol_1 = try(run_C_mcmc(ST$parents, adj_mat = sim.nets[,,1:dim(first.network)[3]],
        num_iters = iterations, verbose = verbose, convex_wt = prob_smooth_param))
      if (class(net_sol_1) != "try-error") { iserror = FALSE }
    }
    
    if (verbose) message ("Begin solution, network 2")
    iserror = TRUE
    while(iserror) {
      net_sol_2 = try(run_C_mcmc(ST$parents, adj_mat = sim.nets[,,dim(first.network)[3]+1:dim(second.network)[3]],
        num_iters = iterations, verbose = verbose, convex_wt = prob_smooth_param))
      if (class(net_sol_2) != "try-error") { iserror = FALSE }
    }
    
    dist_distrib[iter] = tree_distance(net_sol_1$bestTree, net_sol_2$bestTree)
  }

  return(list(true.dist, dist_distrib))
}



