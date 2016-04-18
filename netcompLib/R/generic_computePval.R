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
#' \item 'default-slow' :: Uses the standard hypothesis test, but slower (but cleaner) R code
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
computePval = function(NetS, adja1, adja2, Nobs = 1, pl, output_mode = 'pval', model_type = "default-slow", verbose = TRUE, vbset = c(1,0,0)) {
  stop("Placeholder for documentation purposes")
}


computePval.NetworkStruct = function(NetS, adja1, adja2, Nobs = 1, pl, output_mode = "pval", model_type = "default-slow", verbose = TRUE, vbset = c(1,0,0)) {

  if (model_type != "default") {
    ## This only handles the main default case, with no frills - just outputs chi-square test statistic -- cannot do nodal contributions
    ## Special case - handle this faster with improper functions (uses variables in this function environment, instead of their own local function environment)
    ## This case handles multiple input networks gracefully / quickly... 
    if (getNetType(NetS) == "random") { stop("Cant use fast mode for 'random' structs")}
    if (output_mode == "nodal") { stop("Cannot use 'nodal' method in fast mode") }
    
    compute_one_pair = function() {
      ## Inputs are A1, A2, NetS; assume Nobs = 1; A1/A2 are matrices
      
      ## Count needed stuff 
      if (getNetType(NetS) == "tree") { correct = rep(1, times = getNnodes(NetS)) }
      if (getNetType(NetS) == "block") { correct = NetS@correct }
      
      n = rep(0, times = length(NetS@counts)); x1 = n; x2 = n
      for(j in seq_along(NetS@expand)) {
        x1[j] = sum(A1[NetS@expand[[j]][[1]],NetS@expand[[j]][[2]]], na.rm = TRUE)/correct[j]
        x2[j] = sum(A2[NetS@expand[[j]][[1]],NetS@expand[[j]][[2]]], na.rm = TRUE)/correct[j]
        xc[j] = sum((A1*A2)[NetS@expand[[j]][[1]],NetS@expand[[j]][[2]]], na.rm = TRUE)/correct[j]
        n[j] = sum(!is.na(A1[NetS@expand[[j]][[1]],NetS@expand[[j]][[2]]], na.rm = TRUE))/correct[j]
      }
      xc = x1 + x2
      
      ## Compute likelihoods by group
      llN = compute_loglik_fromPC(x1+x2, 2*n)
      llA = compute_loglik_fromPC(x1, n) + compute_loglik_fromPC(x2, n)
      
      chisqByGroup = -2 * (llN - llA)
      
      if (output_mode == "chisq") { return(sum(chisqByGroup)) }
      
      ## Fast compute df adjustment
      pvals = matrix(0, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
      dyad_dfEst = computeEmpDfAdj(A1, A2, NetS)
      
      px = x1/n; py = x2/n; pc = (x1+x2)/(2*n)
      nums = pc - px*py; dens = sqrt(px * (1-px) * py * (1-py))
      cor_by_dyadgroup = nums/dens
      cor_by_dyadgroup[nums == 0] = 0 ## zero out any 0/0. 
      
      ssadj = sapply(seq_along(cor_by_dyadgroup), function(j) {compute_small_samp_dfadj(fit1$n[j], pc[j], mode = "bound")})
      dfadj = data.frame(coradj = 1 - cor_by_dyadgroup, ssadj = ssadj)
      
      ## Setup pval result matrix. 
      for (j in seq_along(pl$cc_adj)) {
        SEs = 1/sqrt(n)
        
        ## Adjust estimate by some number of SEs (given by cc_adj argument)
        dyad_dfAdj = sapply(dfadj$coradj + SEs * pl$cc_adj[j], function(x) {min(1, x)}) * dfadj$ssadj
        
        for (k in seq_along(pl$thres_ignore)) {
          ## Only keepdyad groups with enough observations (given by thres_ignore argument)
          inds = which(dyad_counts >= pl$thres_ignore[k])
          pvals[j,k] = pchisq(q = sum(chisqByGroup[inds]), df = sum(dyad_dfAdj[inds]), lower.tail = FALSE)
        }
      }
      return(pvals)
    }
    
    ## Set up A1 and A2 and iterate...
    if (is.list(adja1)) {
      reslist = list()
      for (j in seq_along(adjda1)) {
        A1 = adja1[[j]][,,1]; A2 = adja2[[j]][,,1]
        reslist[[j]] = compute_one_pair()
      }
      return(reslist)
    } else { ## Only one adjm pair: 
      A1 = adja1[,,1]; A2 = adja2[,,1]
      return(compute_one_pair())
    }
    
  } else {
    ## TODO: [Update documentation] to show this case
    if (is.list(adja1) & is.list(adja2)) { # If inputs are both lists: 
      if (length(adja1) == length(adja2)) { # Make sure there is a proper amount of comparisons
        return(mapply(function(x,y) { computePval(NetS, x, y, Nobs, pl, output_mode, model_type, verbose, vbset) }, adja1, adja2, SIMPLIFY = FALSE))
        ## TODO: update verbosity settings in this case?
      } else {
        stop("Unable to compute p-values given an unequal amount of pairs of network data")
      }
    } else if (is.list(adja1) + is.list(adja2) == 1) {
      stop("Unable to compute p-values given an unequal amount of pairs of network data")
    }
    
    ## TODO: Change code to work when Nobs != 1... 
    
    ## Compute everything by each node
    adj_abound = abind::abind(adja1, adja2, along = 3)
    by_group = (output_mode %in% c("nodal", "pval"))
    by_node = (output_mode == "nodal")
    
    ## Fit on appropriate model_type and compute log-likelihood contributions
    if (model_type %in% c("default-slow", "densitydiff")) {
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
        dyad_dfAdj = sapply(dyad_dfEst$coradj + SEs * pl$cc_adj[j], function(x) {min(1, x)}) * dyad_dfEst$ssadj
        
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
}


computePval.NetworkStructList = function(NetS, adja1, adja2, Nobs = 1, pl, output_mode = "chisq", model_type = "default-slow", verbose = TRUE, vbset = c(1,0,0)) {
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

