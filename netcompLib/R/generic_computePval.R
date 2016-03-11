##@S Generic function that performs the hypothesis test (computes the p-value for the likelihood ratio test)

setGeneric("computePval", function(NetS, adja1, adja2, Nobs, pl, output_mode, model_type, verbose, vbset) standardGeneric("computePval"))


#' Computes single p-value for likelihood ratio test
#' 
#' This function takes the input network structure (or a list of network structures), and computes test statistics for the likelihood ratio test for the two input adjacency matrices (or arrays). 
#' 
#' This function computes a single p-value per network structure (it doesn't do any p-value aggregation). Thus, the only entries in 'pl' it considers are 'cc_adj' and 'thres_ignore'. 
#' 
#' @param NetS [\code{\link{NetworkStruct}}] :: Model to compute p-value with
#' @param adja1 [matrix/array] :: Adjacency matrix/array
#' @param adja2 [matrix/array] :: Adjacency matrix/array
#' @param Nobs [int] :: Number of network observations per class (default = 1)
#' @param pl [list] :: Simulation/Testing parameters, set by set_sim_param
#' @param output_mode [char] :: How to output results? Options are: 
#' \itemize{
#' \item 'chisq' :: Only gives the asymptotically chi-square test statistic
#' \item 'nodal' :: Outputs break down of test statistic into each vertex
#' \item 'pval' :: Outputs p-value based on chi-square distribution
#' }
#' @param model_type [char] :: What model type to fit? Options are: 
#' \itemize{
#' \item 'default' :: Uses the standard hypothesis test
#' \item 'densitydiff' :: Uses version where null allows for a global additive parameter in probability difference
#' \item 'correlated' :: Uses version where null allows for correlated models
#' }
#' @param verbose [logical] :: Prints progress if TRUE
#' @param vbset [vector-int] :: See \link{OutputSettings}
#' Regular-verbose -- outputs general information for processing the NetworkStructList case
#' High-verbose -- outputs dots for each struct processed
#' 
#' @return [] :: A matrix (or a list of matrices) of p-values (depending on the testing parameters)
#' 
#' @export
#' 
computePval = function(NetS, adja1, adja2, Nobs = 1, pl, output_mode = 'pval', model_type = "default", verbose = TRUE, vbset = c(1,0,0)) {
  stop("Placeholder for documentation purposes")
}


computePval.NetworkStruct = function(NetS, adja1, adja2, Nobs = 1, pl, output_mode = "pval", model_type = "default", verbose = TRUE, vbset = c(1,0,0)) {
  
  ## TODO: [Update documentation] to show this case
  if (is.list(adja1) & is.list(adja2) & (length(adja1) == length(adja2))) {
    return(mapply(function(x,y) { computePval(NetS, x, y, Nobs, pl, output_mode, model_type, verbose, vbset) }, adja1, adja2, SIMPLIFY = FALSE))
    ## TODO: update verbosity settings in this case?
  }
  
  ## TODO: Change code to work when Nobs != 1... 
  
  ## Compute everything by each node
  adj_abound = abind::abind(adja1, adja2, along = 3)
  by_group = (output_mode %in% c("nodal", "pval"))
  by_node = (output_mode == "nodal")
  
  ## Fit on appropriate model_type and compute log-likelihood contributions
  if (model_type %in% c("default", "densitydiff")) {
    fitN = fitModel(NetS, adj_abound, mode = model_type)
    fit1 = fitModel(NetS, adja1); fit2 = fitModel(NetS, adja2)
    
    llN = computeLik(fitN, adj_abound, by_node, by_group)
    llA = addComputeLik(computeLik(fit1, adja1, by_node, by_group), 
                        computeLik(fit2, adja2, by_node, by_group))
  } else if (model_type == "correlated") {
    fitN = fitModel(NetS, adj_abound, mode = "corr-global-null")
    fitA = fitModel(NetS, adj_abound, mode = "corr-global")
    
    llN = computeLik(fitN, adj_abound, by_node, by_group)
    llA = computeLik(fitA, adj_abound, by_node, by_group)
  }
  
  ## Compute chi-square test statistic
  chisq = -2 * (llN$sum - llA$sum)
  
  ## Compute p-value if desired 
  if (output_mode == "chisq") {
    return(chisq)
    
  } else if (output_mode %in% c("pval", "nodal")) {
    chisqByGroup = -2 * (llN$group_ll - llA$group_ll)
    
    ## Setup pval result matrix. 
    pvals = matrix(0, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
    dyad_counts = llA$group_size
    
    dyad_dfEst = computeEmpDfAdj(adja1, adja2, NetS)
    for (j in seq_along(pl$cc_adj)) {
      ## TODO: Add the computation to paper / thesis document
      SEs = 1/sqrt(dyad_counts)
      
      ## Adjust estimate by some number of SEs (given by cc_adj argument)
      dyad_dfAdj = sapply(dyad_dfEst + SEs * pl$cc_adj[j], function(x) {min(1, x)})
      
      for (k in seq_along(pl$thres_ignore)) {
        ## Only keepdyad groups with enough observations (given by thres_ignore argument)
        inds = which(dyad_counts >= pl$thres_ignore[k])
        
        pvals[j,k] = pchisq(q = sum(chisqByGroup[inds]), df = sum(dyad_dfAdj[inds]), lower.tail = FALSE)
      }
    }
    
    if (output_mode == "nodal") {
      chisqByNode = -2 * (llN$bynode - llA$bynode)
      
      return(list(chisq = chisq, pvals = pvals, nodecontrib = chisqByNode))
    } else {
      return(pvals)
    }
  } else {
    stop("Invalid output_mode")
  }
}


computePval.NetworkStructList = function(NetS, adja1, adja2, Nobs = 1, pl, output_mode = "chisq", model_type = "default", verbose = TRUE, vbset = c(1,0,0)) {
  if (verbose & vbset[1] > 0) { 
    cat("\n", stringr::str_pad(string = "", width = vbset[3], pad = "-"), date(), "-- Fitting on", length(NetS@models), "random structures")
  }
  
  ## Call appropriate computePval for each individual NetworkStruct
  res = lapply(NetS@models, function(x) { 
    if (verbose & vbset[1] > 0 & vbset[2] > 0) { cat(".") }
    computePval(x, adja1, adja2, Nobs, pl, output_mode, model_type, verbose, vs_new) 
  } )
  return(res)
}


# setMethod ---------------------------------------------------------------
setMethod("computePval", signature(NetS = "NetworkStruct"), computePval.NetworkStruct)
setMethod("computePval", signature(NetS = "NetworkStructList"), computePval.NetworkStructList)

