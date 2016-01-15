

#' Fit SBM given adjacency matrix
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param Nobs [int] :: Number of observations
#' 
#' @return [\code{\link{NetworkModel}}] :: Best found fitted model
#' 
#' @export
#' 
fit_SBM = function(adjm, Nobs = 1, control_list = set_fit_param()) { 
  cl = control_list ## Shorten control_list variable name
  
  ## TODO: Make this take multiple inputs? adjacency arrays, multiple observations? Need to alter further arguments to make this work! 
  
  ## method = 'spectral' or 'mf'
  ## method 'mf' for mean field EM approach, 'spectral' for just doing spectral clustering
  if (cl$SBM_method == "spectral") {
    ## For spectral clustering, simply apply spectral clustering 
    return(specClust(adjm, cl$SBM_Nclass, cl$Ntries))
    
    
  } else if (cl$SBM_method == "mf") {
    ## For EM algo, use 'start' to decide how to initialize! If using spectral clustering to start, don't need to rerun the eigenvector computation part -- only need to rerun the k-means part! 
    ## start can take on values of 'spectral' or 'random'
    
    for(j in 1:cl$Ntries) {
      N = nrow(adjm) ## This is the number of nodes
      
      best_loglik = -Inf
      best_res = NULL
      if (cl$SBM_start == 'random') {
        ## Random start
        nodeps = rep(1/cl$SBM_Nclass, length = cl$SBM_Nclass)
        edgeps = symmetrize_mat(matrix(sum(adjm, na.rm = TRUE) / (Nobs * N * (N-1) / 2) * runif(n = cl$SBM_Nclass * cl$SBM_Nclass, min = 0.1, max = 0.9), nrow = cl$SBM_Nclass))
        
        H = matrix(0, nrow = N, ncol = cl$SBM_Nclass)
        PHI = matrix(runif(cl$SBM_Nclass*N, min = 0.1, max = 0.9), nrow = N)
        PHI = PHI / rowSums(PHI)
      } else if (cl$SBM_start == 'spectral') {
        
        ## Do spectral clustering once to give a rough start
        res = specClust(adjm, cl$SBM_Nclass, 2)
        clusts = res@groups
        
        nodeps = sapply(1:cl$SBM_Nclass, function(x) {mean(clusts == x)})
        edgeps = res@probmat
        
        H = matrix(0, nrow = N, ncol = cl$SBM_Nclass)
        PHI = matrix(runif(cl$SBM_Nclass*N, min = 0.1, max = 0.9), nrow = N)
        for(j in 1:cl$SBM_Nclass) {
          PHI[which(j == clusts),j] = PHI[which(j == clusts),j] + 2
        }
        PHI = PHI / rowSums(PHI)
      }
      
      edgeps[edgeps < 0.01] = .01; edgeps[edgeps > 0.99] = .99
      results = EM_SBM_mf(adjm = adjm, Nobs = Nobs, nodeps = nodeps, edgeps = edgeps, H = H, PHI = PHI, 
                          Niter = cl$SBM_EM_Niter, stop_thres = cl$SBM_EM_stopthres, verbose = cl$verbose)
      # newm = NetworkModel(set_model_param(Nnodes = 30, block_assign = clusts, block_probs = edgeps))
      # computeLik(newm, adjm)$sum
      loglik = computeLik(results$model, adja = adjm)$sum
      if (cl$verbose > 0) { print(loglik) }
      
      if (loglik > best_loglik) {
        best_loglik = loglik
        best_res = results
      }
    }
    
    return(best_res)
  }
}


#' Runs spectral clustering on an input matrix
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param Nclass [int] :: Number of clusters
#' @param Ntries [int] :: Number of attempts
#' 
#' @return [\code{\link{NetworkModel}}] :: Output bestfit model
#' 
#' @export
#' 
specClust = function(adjm, Nclass, Ntries) {
  diag(adjm) = 0 # Should be 0 already, but this won't hurt. 
  
  ## IF there is missing data, use cheap algorithm to do matrix completing
  if (any(is.na(adjm))) {
    ## Fill with row average if available, if not, fill with global average. 
    N = nrow(adjm)
    below_thres = NULL
    for(j in seq_len(N)) {
      r = which(is.na(adjm[j,]))
      if (length(r) > N * .8) {
        below_thres = c(below_thres, r)
      } else {
        adjm[j,r] = sum(adjm[j,], na.rm = TRUE)/(N-1)
      }
    } 
    
    overall_density = sum(adjm[-r,-r], na.rm = TRUE) / ((N-r) * (N-r-1))
    for(j in below_thres) {
      adjm[j,] = overall_density
    }
    diag(adjm) = 0 # Should be 0 already, but this won't hurt. 
  }
  
  best_res = NULL
  best_loglik = -Inf
  
  for(j in 1:Ntries) {
    
    clusts = kmeans(scale(eigen(adjm)$vectors[,1:Nclass]), centers = Nclass, nstart = 10)$cluster
    nodeps = sapply(1:Nclass, function(x) { sum(x == clusts) }) / length(clusts)
    nets = NetworkStruct(set_model_param(Nnodes = length(clusts), type = 'block', block_nclass = Nclass, block_assign = clusts))
    model = fitModel(nets, adjm)
    loglik = computeLik(model, adja = adjm)$sum
    if (loglik > best_loglik) {
      best_loglik = loglik
      best_res = model
    }
  }
  return(best_res)
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
#' \item model -- [\link{\code{NetworkModel}}] :: Estimated Network Model
#' }
#' 
#' @export
#' 
EM_SBM_mf = function(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, verbose) {
  
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

