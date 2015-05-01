##@S Generic function that performs the hypothesis test (computes the p-value for the likelihood ratio test)

setGeneric("computeLik", function(NetM, adja, loglik, na.rm) standardGeneric("computeLik"))

#' Compute Likelihood of Network given a model
#' 
#' @param NetM Network Model 
#' @param adja Adjacency matrix/array
#' @param loglik if true -- give log-likelihood instead of likelihood
#' @param na.rm if true -- ignores NAs in the adjacency matrix
#' 
#' @return likelihood (or log-likelihood)
#' 
#' @export
#' 
computeLik = function(NetM, adja, loglik = TRUE, na.rm = TRUE) {
  stop("Placeholder for generic function -- this call is meaningless for a generic NetworkModel")
}


computeLik.NetworkModel = function(NetM, adja, loglik = TRUE, na.rm = TRUE) {
  ## TODO: [Improvement] Can handle fully-missing edges, but cannot handle partially-missing edges. 
  if (length(dim(adja)) == 2) {
    adjm = adja; Nobs = 1
  } else if (length(dim(adja)) == 3) {
    adjm = apply(adja, c(1,2), sum); Nobs = dim(adja)[3]
  } else {
    stop("Invalid input 'adja' (not a 2D or 3D array)")
  }
  
  eps = getEdgeProbMat(NetM)
  resm = adjm * log(eps) + (Nobs - adjm) * log(1 - eps)
  res = sum(resm[upper.tri(x = resm, diag = FALSE)], na.rm = na.rm)
    
  if (loglik) {
    return(res) 
  } else {
    return(exp(res))
  }
}

# setMethod ---------------------------------------------------------------
setMethod("computeLik", signature(NetM = "NetworkModel"), computeLik.NetworkModel)

