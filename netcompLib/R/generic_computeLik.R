##@S Generic function that performs the hypothesis test (computes the p-value for the likelihood ratio test)

setGeneric("computeLik", function(NetM, adja, loglik, by_node, na.rm) standardGeneric("computeLik"))

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computeLik)
#' Compute Likelihood of Network given a model
#' 
#' @param NetM Network Model 
#' @param adja Adjacency matrix/array
#' @param loglik if true -- give log-likelihood instead of likelihood
#' @param by_node temp
#' @param na.rm if true -- ignores NAs in the adjacency matrix
#' 
#' @return likelihood (or log-likelihood)
#' 
#' @export
#' 
computeLik = function(NetM, adja, loglik = TRUE, by_node = FALSE, na.rm = TRUE) {
 ## by node -- true -> return a vector of log likelihoods. (summed up by node)
  stop("Placeholder for generic function -- this call is meaningless for a generic NetworkModel")
}


computeLik.NetworkModel = function(NetM, adja, loglik = TRUE, by_node = FALSE, na.rm = TRUE) {
  if (length(dim(adja)) == 2) {
    adjm = adja; Nobs = 1
  } else if (length(dim(adja)) == 3) {
    adjm = apply(adja, c(1,2), sum); Nobs = dim(adja)[3]
  } else {
    stop("Invalid input 'adja' (not a 2D or 3D array)")
  }
  
  eps = getEdgeProbMat(NetM)
  resm = adjm * log(eps) + (Nobs - adjm) * log(1 - eps)
  diag(resm) = 0
  tempres = apply(resm, 1, sum, na.rm = na.rm)/2
  
  if (by_node) { res = tempres } else { res = sum(tempres, na.rm = na.rm) }
  # old cleaner code... 
  # res = sum(resm[upper.tri(x = resm, diag = FALSE)], na.rm = na.rm)
    
  if (loglik) {
    return(res) 
  } else {
    return(exp(res))
  }
}

computeLik.NetworkModelPair = function(NetM, adja, loglik = TRUE, by_node = FALSE, na.rm = TRUE) {
  ## TODO: apply for general cases -- right now, ASSUMES adja is dimension 2. 
  
  if (NetM@model_type == "correlated") {
    pt = hCorr_paramToProb(NetM@addl_param$c_param_corr, NetM@addl_param$c_param_a, NetM@addl_param$c_param_b)
    adjm = adja[,,1] * 10 + adja[,,2]
    matches_x = matrix(match(adjm, as.numeric(colnames(pt))), nrow = nrow(adjm))
    matches_group = matrix(match(getEdgeProbMat(NetM, "group")[[1]], as.numeric(NetM@addl_param$c_names)), nrow = nrow(adjm))
    
    llbynode = sapply(seq_len(nrow(adjm)), function(i) { sum(log(mapply(function(x,y) {pt[x,y]}, y=matches_x[i,],x=matches_group[i,])[-i]))/2 })
    
    if (by_node) { res = llbynode } else { res = sum(llbynode, na.rm = na.rm) }
    if (loglik) { return(res) } else { return(exp(res)) }
  }

  return(computeLik(NetM@m1, adja[,,1], loglik, by_node, na.rm) + 
           computeLik(NetM@m2, adja[,,2], loglik, by_node, na.rm))
}

# setMethod ---------------------------------------------------------------
setMethod("computeLik", signature(NetM = "NetworkModel"), computeLik.NetworkModel)
setMethod("computeLik", signature(NetM = "NetworkModelPair"), computeLik.NetworkModelPair)

