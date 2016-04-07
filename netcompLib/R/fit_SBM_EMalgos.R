
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_mf_chooser)
#' Fit SBM using mean field approx
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param Nobs [int] :: Number of observations per dyad
#' @param nodeps [vec-double] :: Initial node probabilities
#' @param edgeps [matrix-double] :: Initial between-block probabilities
#' @param H [matrix-double] :: Initial value of H
#' @param PHI [matrix-double] :: Initial value of PHI
#' @param Niter [int] :: Number of iterations max
#' @param stop_thres [double] :: Stop condition
#' @param mode temp
#' @param verbose [logical] :: Output stuff?
#' 
#' @return [list] :: EM algorithm output, with the following elements: 
#' \itemize{
#' \item nodeps -- [vec-double] :: Estimated class probabilities
#' \item edgeps -- [matrix-double] :: Estimated between-class edge probabilities
#' \item classes -- [vec-int] :: Estimated class assignments
#' \item model -- [\code{\link{NetworkModel}}] :: Estimated Network Model
#' }
#' 
#' @export
#' 
EM_SBM_mf_chooser = function(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, mode, verbose) {
  # mode = default -> fit using only one edgepm. 
  # mode = no-node-default -> fit only using one edgepm, but ignoring nodeps
  
  if (mode == 'default') {
    return(EM_SBM_mf_default(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, verbose))
  } else if (mode == "no-node-default") {
    return(EM_SBM_mf_default_nonode(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, verbose))
  }
}


  

#' Fit SBM using mean field approx
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param Nobs [int] :: Number of observations per dyad
#' @param nodeps [vec-double] :: Initial node probabilities
#' @param edgeps [matrix-double] :: Initial between-block probabilities
#' @param H [matrix-double] :: Initial value of H
#' @param PHI [matrix-double] :: Initial value of PHI
#' @param Niter [int] :: Number of iterations max
#' @param stop_thres [double] :: Stop condition
#' @param verbose [logical] :: Output stuff?
#' 
#' @return [list] :: EM algorithm output, with the following elements: 
#' \itemize{
#' \item nodeps -- [vec-double] :: Estimated class probabilities
#' \item edgeps -- [matrix-double] :: Estimated between-class edge probabilities
#' \item classes -- [vec-int] :: Estimated class assignments
#' \item model -- [\code{\link{NetworkModel}}] :: Estimated Network Model
#' }
#' 
#' @export
#' 
EM_SBM_mf_default = function(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, verbose) {
  
  N = nrow(adjm)
  Nclass = nrow(edgeps)
  
  compute_objfx = function() {
    runsum = 0
    for(r in 1:Nclass) { for (s in 1:Nclass) {
      s1 = adjm * log(edgeps[r,s]) + (Nobs - adjm) * log(1 - edgeps[r,s])
      s1 = s1 * (PHI[,r] %*% t(PHI[,s]))
      diag(s1) = 0
      runsum = runsum + sum(s1, na.rm = TRUE)
    }}
    for(r in 1:Nclass) {
      runsum = runsum + sum(PHI[,r] * (log(nodeps[r]) - log(PHI[,r])), na.rm = TRUE)
    }
    return(runsum)
  }
  
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
    if (verbose > 3) { print(compute_objfx()) }
    
    ## M step
    nodeps_new = apply(PHI, 2, sum, na.rm = TRUE) / N
    edgeps_new = edgeps * 0
    for(r in 1:Nclass) {
      for(s in r:Nclass) {
        Psq = PHI[,r,drop = FALSE] %*% t(PHI[,s,drop = FALSE])
        num = adjm * Psq
        den = matrix(as.numeric(!is.na(adjm)), nrow = N) * Nobs * Psq
        edgeps_new[r,s] = max(min(sum(num[lower.tri(num)], na.rm = TRUE) / sum(den[lower.tri(den)], na.rm = TRUE), 0.999), 0.001) ## bound the prob by 0.999 and 0.001 to prevent weird stuff?
        edgeps_new[!is.finite(edgeps_new)] = .999
      }
    }
    edgeps_new = symmetrize_mat(edgeps_new)    
    if (verbose > 2) { print(compute_objfx()) }
    
    ## Check for convergence
    ## Compute change in parameter estimates
    #print(edgeps_new)
    #print(nodeps_new)
    delta = sum(abs(edgeps_new - edgeps)) + sum(abs(nodeps_new - nodeps), na.rm = TRUE)
    if (verbose > 1) { cat("Iteration ", I, " ----- change in edgeps and nodeps = ", delta, "\n", sep = "") }
    #if (verbose > 1) { cat("\t\t Log-likelihood: ", compute_sbm_loglik(class_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), adjm = adjm, nodeps = nodeps_new, edgeps = edgeps_new, Nobs = Nobs), "\n") }
    
    ## Update old parameters
    edgeps = edgeps_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    # print(delta)
    if (delta < stop_thres) { break }
  }
  
  return(list(nodeps = nodeps, edgeps = edgeps, PHI = PHI, 
              classes = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), nsteps = I,
              model = NetworkModel(set_model_param(
                Nnodes = N, type = 'block', block_nclass = Nclass, 
                block_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }),
                block_probs = edgeps))))
}



#' Fit SBM using mean field approx
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param Nobs [int] :: Number of observations per dyad
#' @param nodeps [vec-double] :: Initial node probabilities
#' @param edgeps [matrix-double] :: Initial between-block probabilities
#' @param H [matrix-double] :: Initial value of H
#' @param PHI [matrix-double] :: Initial value of PHI
#' @param Niter [int] :: Number of iterations max
#' @param stop_thres [double] :: Stop condition
#' @param verbose [logical] :: Output stuff?
#' 
#' @return [list] :: EM algorithm output, with the following elements: 
#' \itemize{
#' \item nodeps -- [vec-double] :: Estimated class probabilities
#' \item edgeps -- [matrix-double] :: Estimated between-class edge probabilities
#' \item classes -- [vec-int] :: Estimated class assignments
#' \item model -- [\code{\link{NetworkModel}}] :: Estimated Network Model
#' }
#' 
#' @export
#' 
EM_SBM_mf_default_nonode = function(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, verbose) {
  
  N = nrow(adjm)
  Nclass = nrow(edgeps)
  
  compute_objfx = function() {
    runsum = 0
    for(r in 1:Nclass) { for (s in 1:Nclass) {
      s1 = adjm * log(edgeps[r,s]) + (Nobs - adjm) * log(1 - edgeps[r,s])
      s1 = s1 * (PHI[,r] %*% t(PHI[,s]))
      diag(s1) = 0
      runsum = runsum + sum(s1, na.rm = TRUE)
    }}
    return(runsum)
  }
  
  for(I in 1:Niter) {
    ## E step
    H = H * 0
    for(r in 1:Nclass) { for(s in 1:Nclass) {
      h1 = adjm * log(edgeps[r,s]) + (Nobs - adjm) * log(1 - edgeps[r,s])
      H[,r] = H[,r] + sapply(1:N, function(x) {sum(h1[x,-x] * PHI[-x,s], na.rm = TRUE)})
    }}
    H = H + abs(max(H)) #rescale to prevent exponentiation errors
    peH = exp(H)
    PHI = peH / rowSums(peH)
    if (verbose > 3) { print(compute_objfx()) }
    
    ## M step
    nodeps_new = apply(PHI, 2, sum, na.rm = TRUE) / N
    edgeps_new = edgeps * 0
    for(r in 1:Nclass) {
      for(s in r:Nclass) {
        Psq = PHI[,r,drop = FALSE] %*% t(PHI[,s,drop = FALSE])
        num = adjm * Psq
        den = matrix(as.numeric(!is.na(adjm)), nrow = N) * Nobs * Psq
        edgeps_new[r,s] = max(min(sum(num[lower.tri(num)], na.rm = TRUE) / sum(den[lower.tri(den)], na.rm = TRUE), 0.999), 0.001) ## bound the prob by 0.999 and 0.001 to prevent weird stuff?
        edgeps_new[!is.finite(edgeps_new)] = .999
      }
    }
    edgeps_new = symmetrize_mat(edgeps_new)    
    if (verbose > 2) { print(compute_objfx()) }
    
    ## Check for convergence
    ## Compute change in parameter estimates
    #print(edgeps_new)
    #print(nodeps_new)
    delta = sum(abs(edgeps_new - edgeps)) + sum(abs(nodeps_new - nodeps), na.rm = TRUE)
    if (verbose > 1) { cat("Iteration ", I, " ----- change in edgeps and nodeps = ", delta, "\n", sep = "") }
    #if (verbose > 1) { cat("\t\t Log-likelihood: ", compute_sbm_loglik(class_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), adjm = adjm, nodeps = nodeps_new, edgeps = edgeps_new, Nobs = Nobs), "\n") }
    
    ## Update old parameters
    edgeps = edgeps_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    # print(delta)
    if (delta < stop_thres) { break }
  }
  
  return(list(nodeps = nodeps, edgeps = edgeps, PHI = PHI, 
              classes = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), nsteps = I,
              model = NetworkModel(set_model_param(
                Nnodes = N, type = 'block', block_nclass = Nclass, 
                block_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }),
                block_probs = edgeps))))
}



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_mf_twoedgeps_nonode)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_mf_twoedgeps_nonode)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_mf_twoedgeps_nonode)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_mf_twoedgeps_nonode)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_mf_twoedgeps_nonode)
#' Fit SBM using mean field approx
#' 
#' @param edgeps [matrix-double] :: Initial between-block probabilities
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param adjm1 temp
#' @param adjm2 temp
#' @param Nobs [int] :: Number of observations per dyad
#' @param nodeps [vec-double] :: Initial node probabilities
#' @param edgeps1 temp
#' @param edgeps2 temp
#' @param H [matrix-double] :: Initial value of H
#' @param PHI [matrix-double] :: Initial value of PHI
#' @param Niter [int] :: Number of iterations max
#' @param stop_thres [double] :: Stop condition
#' @param verbose [logical] :: Output stuff?
#' 
#' @return [list] :: EM algorithm output, with the following elements: 
#' \itemize{
#' \item nodeps -- [vec-double] :: Estimated class probabilities
#' \item edgeps -- [matrix-double] :: Estimated between-class edge probabilities
#' \item classes -- [vec-int] :: Estimated class assignments
#' \item model -- [\code{\link{NetworkModel}}] :: Estimated Network Model
#' }
#' 
#' @export
#' 
EM_SBM_mf_twoedgeps_nonode = function(adjm1, adjm2, Nobs, nodeps, edgeps1, edgeps2, H, PHI, Niter, stop_thres, verbose) {
  
  N = nrow(adjm1)
  Nclass = nrow(edgeps1)
  
  compute_objfx = function() {
    runsum = 0
    for(r in 1:Nclass) { for (s in 1:Nclass) {
      s1 = adjm1 * log(edgeps1[r,s]) + adjm2 * log(edgeps2[r,s]) +
        (Nobs - adjm1) * log(1 - edgeps1[r,s]) + (Nobs - adjm2) * log(1 - edgeps2[r,s])
      s1 = s1 * (PHI[,r] %*% t(PHI[,s]))
      diag(s1) = 0
      runsum = runsum + sum(s1, na.rm = TRUE)
    }}
    return(runsum)
  }
  
  for(I in 1:Niter) {
    ## E step
    H = H * 0
    for(r in 1:Nclass) { for(s in 1:Nclass) {
      h1 = adjm1 * log(edgeps1[r,s]) + (Nobs - adjm1) * log(1 - edgeps1[r,s]) + 
        adjm2 * log(edgeps2[r,s]) + (Nobs - adjm2) * log(1 - edgeps2[r,s])
      H[,r] = H[,r] + sapply(1:N, function(x) {sum(h1[x,-x] * PHI[-x,s], na.rm = TRUE)})
    }}
    H = H + abs(max(H)) #rescale to prevent exponentiation errors
    peH = exp(H)
    PHI = peH / rowSums(peH)
    if (verbose > 3) { print(compute_objfx()) }
    
    ## M step
    nodeps_new = apply(PHI, 2, sum, na.rm = TRUE) / N
    edgeps1_new = edgeps1 * 0
    edgeps2_new = edgeps2 * 0
    for(r in 1:Nclass) {
      for(s in r:Nclass) {
        Psq = PHI[,r,drop = FALSE] %*% t(PHI[,s,drop = FALSE])
        num1 = adjm1 * Psq
        num2 = adjm2 * Psq
        den = matrix(as.numeric(!is.na(adjm1)), nrow = N) * Nobs * Psq
        edgeps1_new[r,s] = max(min(sum(num1[lower.tri(num1)], na.rm = TRUE) / sum(den[lower.tri(den)], na.rm = TRUE), 0.999), 0.001) ## bound the prob by 0.999 and 0.001 to prevent weird stuff?
        edgeps1_new[!is.finite(edgeps1_new)] = .999
        edgeps2_new[r,s] = max(min(sum(num2[lower.tri(num2)], na.rm = TRUE) / sum(den[lower.tri(den)], na.rm = TRUE), 0.999), 0.001) ## bound the prob by 0.999 and 0.001 to prevent weird stuff?
        edgeps2_new[!is.finite(edgeps2_new)] = .999
        
      }
    }
    edgeps1_new = symmetrize_mat(edgeps1_new)    
    edgeps2_new = symmetrize_mat(edgeps2_new)    
    
    if (verbose > 2) { print(compute_objfx()) }
    
    ## Check for convergence
    ## Compute change in parameter estimates
    #print(edgeps_new)
    #print(nodeps_new)
    delta = sum(abs(edgeps1_new - edgeps1)) + sum(abs(edgeps2_new - edgeps2)) + sum(abs(nodeps_new - nodeps), na.rm = TRUE)
    if (verbose > 1) { cat("Iteration ", I, " ----- change in edgeps and nodeps = ", delta, "\n", sep = "") }
    #if (verbose > 1) { cat("\t\t Log-likelihood: ", compute_sbm_loglik(class_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), adjm = adjm, nodeps = nodeps_new, edgeps = edgeps_new, Nobs = Nobs), "\n") }
    
    ## Update old parameters
    edgeps1 = edgeps1_new
    edgeps2 = edgeps2_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    # print(delta)
    if (delta < stop_thres) { break }
  }
  
  return(list(nodeps = nodeps, edgeps1 = edgeps1, edgeps2 = edgeps2, PHI = PHI, 
              classes = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), nsteps = I,
              model = NetworkModel(set_model_param(
                Nnodes = N, type = 'block', block_nclass = Nclass, 
                block_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }),
                block_probs = edgeps1))))
}






