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
  
  ## TODO: Change code to work when Nobs != 1... 
  
  ## Compute everything by each node
  adj_abound = abind::abind(adja1, adja2, along = 3)
  
  ## Fit on appropriate model_type
  if (model_type %in% c("default", "densitydiff")) {
    fitN = fitModel(NetS, adj_abound, mode = model_type)
    fit1 = fitModel(NetS, adja1); fit2 = fitModel(NetS, adja2)
    
    llN = computeLik(fitN, adj_abound, by_node = TRUE)$bynode
    llA = computeLik(fit1, adja1, by_node = TRUE)$bynode + 
      computeLik(fit2, adja2, by_node = TRUE)$bynode
  } else if (model_type == "correlated") {
    fitN = fitModel(NetS, adj_abound, mode = "corr-global-null")
    fitA = fitModel(NetS, adj_abound, mode = "corr-global")
    
    llN = computeLik(fitN, adj_abound, by_node = TRUE)$bynode
    llA = computeLik(fitA, adj_abound, by_node = TRUE)$bynode
  }
  
  ## Compute chi-square test statstic
  chisqByNode = -2 * (llN - llA)
  chisq = sum(chisqByNode) 
  
  ## Compute p-value if desired 
  if (output_mode == "chisq") {
    return(chisq)
  } else if (output_mode %in% c("pval", "nodal")) {
    
    ## Setup pval result matrix. 
    pvals = matrix(0, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
    dyad_counts = computeLik(fit1, adja1, by_node = TRUE)$group_size
    
    dfadj_perdyad = computeEmpDfAdj(adja1[,,1], adja2[,,1], NetS)
    for(j in seq_along(pl$cc_adj)) {
      for(k in seq_along(pl$thres_ignore)) {
        ## TODO Fill in this, using dfadj_perdyad. 
        
      }
    }
    
    ## TODO: [Update] fix implmenetation of parameter list; since set_sim_param has been updated. 
    ## 'pl' not even used here... should it be? -- yes, when computing p-values. 
    ## New implmenetation -- this can be used for all NetworkStruct (as long as its not a list)
    ## look at cc_adj, thres_ignore. If these are non-vectors, pvals should be returned as a matrix. 
    
    if (output_mode == "nodal") {
      
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

