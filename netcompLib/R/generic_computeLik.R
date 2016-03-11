##@S Generic function that performs the hypothesis test (computes the p-value for the likelihood ratio test)

setGeneric("computeLik", function(NetM, adja, by_node, by_group, na.rm) standardGeneric("computeLik"))


#' Compute Log-likelihood of a network on a specific model
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model to compute log-likelihood on
#' @param adja [matrix/array] :: Adjacency matrix or array
#' @param by_node [logical] :: If TRUE, returns by_node breakdown of log-likelihood
#' @param by_group [logical] :: If TRUE, returns breakdown of log-likelihood by edge groups
#' @param na.rm [logical] :: If TRUE, ignores NAs in the adjacency matrix
#' 
#' @return [list] :: Log-likelihood in various formats, stored in a list as follows:
#' \itemize{
#' \item sum -- [double] :: Entire log-likelihood
#' \item bynode -- [vector-double] :: Log-likelihood summed for each node
#' \item group_ll -- [vector-double] :: Log-likelihood summed in each dyad group
#' \item group_size -- [vector-int] :: Dyad group sizes
#' }
#' 
#' @export
#' 
computeLik = function(NetM, adja, by_node = FALSE, by_group = FALSE, na.rm = TRUE) {
  stop("Placeholder for generic function -- this call is meaningless for a generic NetworkModel")
}


computeLik.NetworkModel = function(NetM, adja, by_node = FALSE, by_group = FALSE, na.rm = TRUE) {
  res = list(sum = NULL, bynode = NULL, group_ll = NULL, group_size = NULL)
  
  if (length(dim(adja)) == 2) {
    adjm = adja; Nobs = 1
  } else if (length(dim(adja)) == 3) {
    adjm = apply(adja, c(1,2), sum); Nobs = dim(adja)[3]
  } else {
    stop("Invalid input 'adja' (not a 2D or 3D array)")
  }
  
  ## Compute log-likelihoods by dyad
  epmat = getEdgeProbMat(NetM)
  ll_dyad = adjm * log(epmat) + (Nobs - adjm) * log(1 - epmat)
  diag(ll_dyad) = 0
  
  ll_node = apply(ll_dyad, 1, sum, na.rm = na.rm)/2
  
  ## Process output, computing other quantities as needed
  res$sum = sum(ll_node, na.rm = na.rm)
  if (by_node) { res$bynode = ll_node }
  if (by_group) {
    bl_tri = lower.tri(ll_dyad, diag = FALSE)
    ll_dyad_lower = ll_dyad[bl_tri]
    egmat = getEdgeProbMat(NetM, 'group')[bl_tri]
    res$group_ll = unname(tapply(ll_dyad_lower, egmat, sum))
    res$group_size = unname(tapply(ll_dyad_lower, egmat, function(x) { sum(!is.na(x))}))
  }
  return(res) 
}

computeLik.NetworkModelPair = function(NetM, adja, by_node = FALSE, by_group = FALSE, na.rm = TRUE) {
  ## TODO: apply for general cases -- right now, ASSUMES adja is dimension 2. 
  ## TODO: Allow this to work for the case by_group = TRUE! 
  
  res = list(sum = NULL, bynode = NULL, group_ll = NULL, group_size = NULL)
  
  if (NetM@model_type == "correlated") {
    pt = hCorr_paramToProb(NetM@addl_param$c_param_corr, NetM@addl_param$c_param_a, NetM@addl_param$c_param_b)
    adjm = adja[,,1] * 10 + adja[,,2]
    matches_x = matrix(match(adjm, as.numeric(colnames(pt))), nrow = nrow(adjm))
    matches_group = matrix(match(getEdgeProbMat(NetM, "group")[[1]], as.numeric(NetM@addl_param$c_names)), nrow = nrow(adjm)) ## Fix c_names call here. 
    
    ll_node = sapply(seq_len(nrow(adjm)), function(i) { sum(log(mapply(function(x,y) {pt[x,y]}, y=matches_x[i,],x=matches_group[i,])[-i]))/2 })
    
    if (by_node) { res$bynode = ll_node }
    res$sum = sum(ll_node, na.rm = na.rm)
    
  } else { 
    ## Run computeLik separately and add results
    l1 = computeLik(NetM@m1, adja[,,1], loglik, by_node, na.rm)
    l2 = computeLik(NetM@m2, adja[,,2], loglik, by_node, na.rm)
    for(s in names(res)) { res[s] = l1[s] + l2[s] }
  }
  return(res)
}

# setMethod ---------------------------------------------------------------
setMethod("computeLik", signature(NetM = "NetworkModel"), computeLik.NetworkModel)
setMethod("computeLik", signature(NetM = "NetworkModelPair"), computeLik.NetworkModelPair)

