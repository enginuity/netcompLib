
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (fit_SBM)
#' <What does this function do>
#' 
#' @param adjm temp
#' @param Nobs temp
#' 
#' @return temp
#' 
#' @export
#' 
fit_SBM = function(adjm, Nobs) { 
  
  N = nrow(adjm) ## This is the number of nodes
  nodeps = rep(1/Nclass, length = Nclass)
  edgeps = symmetrize_mat(matrix(sum(adjm, na.rm = TRUE) / (Nobs * N * (N-1) / 2) * runif(n = Nclass * Nclass, min = 0.1, max = 0.9), nrow = Nclass))
  
  
  H = matrix(0, nrow = N, ncol = Nclass)
  PHI = matrix(runif(Nclass*N, min = 0.1, max = 0.9), nrow = N)
  PHI = PHI / rowSums(PHI)
  
  results = EM_SBM_mf(adjm = adjm, Nobs = Nobs, nodeps = nodeps, edgeps = edgeps, H = H, PHI = PHI, 
                      Niter = Niter, stop_thres = stop_thres, verbose = verbose)
  return(results)            
}


#' Fit SBM using mean field approx
#' 
#' @param adjm temp
#' @param Nobs temp
#' @param nodeps temp
#' @param edgeps temp
#' @param H temp
#' @param PHI temp
#' @param Niter temp
#' @param stop_thres temp
#' @param verbose temp
#' 
#' @return temp
#' 
#' @export
#' 
EM_SBM_mf = function(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, verbose) {
  
  N = nrow(adjm)
  Nclass = nrow(edgeps)
  
  for(I in 1:Niter) {
    ## E step
    H = H * 0
    for(r in 1:Nclass) { for(s in 1:Nclass) {
      h1 = adjm * log(edgeps[r,s]) + (Nobs - adjm) * log(1 - edgeps[r,s])
      H[,r] = H[,r] + sapply(1:N, function(x) {sum(h1[x,-x] * PHI[-x,s], na.rm = TRUE)})
    }}
    H = H + abs(max(H)) #rescale to prevent exponentiation errors
    peH = t(nodeps * t(exp(H)))
    PHI = peH / rowSums(peH)
    
    ## M step
    nodeps_new = apply(PHI, 2, sum) / N
    edgeps_new = edgeps * 0
    for(r in 1:Nclass) {
      for(s in r:Nclass) {
        Psq = PHI[,r,drop = FALSE] %*% t(PHI[,s,drop = FALSE])
        num = adjm * Psq
        den = matrix(as.numeric(!is.na(adjm)), nrow = N) * Nobs * Psq
        edgeps_new[r,s] = max(min(sum(num[lower.tri(num)], na.rm = TRUE) / sum(den[lower.tri(den)]), 0.999), 0.001) ## bound the prob by 0.999 and 0.001 to prevent weird stuff?
      }
    }
    edgeps_new = symmetrize_mat(edgeps_new)    
    
    ## Check for convergence
    ## Compute change in parameter estimates
    delta = sum(abs(edgeps_new - edgeps)) + sum(abs(nodeps_new - nodeps))
    if (verbose > 0) { cat("Iteration ", I, " ----- change in edgeps and nodeps = ", delta, "\n", sep = "") }
    if (verbose > 1) { cat("\t\t Log-likelihood: ", compute_sbm_loglik(class_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), adjm = adjm, nodeps = nodeps_new, edgeps = edgeps_new, Nobs = Nobs), "\n") }
    
    ## Update old parameters
    edgeps = edgeps_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    if (delta < stop_thres) { break }
  }
  
  return(list(nodeps = nodeps, edgeps = edgeps, PHI = PHI, 
              classes = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), nsteps = I))
}
