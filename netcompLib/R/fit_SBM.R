

#' Fit a SBM given an input adjacency matrix (pair) and possibly starting model
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param Nobs [int] :: Number of observations
#' @param control_list [list] :: List of control parameters -- see \code{\link{set_fit_param}}
#' 
#' @return [\code{\link{NetworkModel}}] :: Best found fitted model
#' 
#' @export
#' 
fit_SBM = function(adjm, Nobs = 1, control_list = set_fit_param()) { 
  
  ## TODO: Implement CV to help choose parameters... 
  ## TODO: Make this take multiple inputs? adjacency arrays, multiple observations? Need to alter further arguments to make this work! 
  
  # TODO: Implement bic?
  #   best_loglik = best_loglik - (Nclass*(Nclass-1)/2 + Nclass)*log(sum(!is.na(adjm))) ## (n + n(n-1)/2) [parameters] * log(n) [from BIC]
  #   best_k = which.max(best_loglik)
  
  ## 1-2: If desired
  ## 1. Do matrix completion
  ## -- This will only ever be done once, as the results are deterministic
  ## 2. Do spectral clustering
  ## -- can rerun this (since it gives a different warm start)
  ## 3. Throw stuff into EM algorithm
  
  ## For EM algo, use 'start' to decide how to initialize! If using spectral clustering to start, don't need to rerun the eigenvector computation part -- only need to rerun the k-means part! 
  ## start can take on values of 'spectral' or 'random'
  
  cl = control_list ## Shorten control_list variable name
  if (cl$SBM_start != "random") {
    ## (1) Matrix Completion (if needed)
    ## Some error checking / warnings
    
    if (cl$SBM_MC_Neigenvecs > 0 & cl$SBM_MC_laplacian) { print("WARNING: Using Laplacian and largest eigenvectors... is this intended?")}
    if (cl$SBM_MC_Neigenvecs < 0 & !cl$SBM_MC_laplacian) { print("WARNING: Using adjacnecy matrices and smallest eigenvectors... is this intended?")}
    if (cl$SBM_MC_Neigenvecs == 0) { stop("Must use a nonzero number of eigenvectors") }
    
    if (cl$SBM_start == "spectral-mean") {
      cmMethod = "rcmeans"
    } else if (cl$SBM_start == "spectral-complete") {
      cmMethod = "complete"
    }
    compEVS = completeMatrix(adjm = adjm, method = cmMethod, laplacian = cl$SBM_MC_laplacian, eigenvecs = cl$SBM_MC_Neigenvecs, softImpute_maxit = cl$SBM_MC_softImpute_maxit, softImpute_thresh = cl$SBM_MC_softImpute_thresh, softImpute_rankmax = cl$SBM_MC_softImpute_rankmax)
  }
  
  ## (2) Spectral Clustering and (3) EM algorithm
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
    } else if (length(grep("spectral", cl$SBM_start)) > 0) {
      
      ## Do spectral clustering once to give a rough start
      res = specClust(compEVS, Nclass = cl$SBM_SC_Nclass, NStart = cl$SBM_SC_Nstart)
      res = fitModel(res, adjm)
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
    results = EM_SBM_mf_chooser(adjm = adjm, Nobs = Nobs, nodeps = nodeps, edgeps = edgeps, H = H, PHI = PHI, 
                        Niter = cl$SBM_EM_Niter, stop_thres = cl$SBM_EM_stopthres, mode = cl$SBM_EM_mode, verbose = cl$verbose)
    
    loglik = computeLik(results$model, adja = adjm)$sum
    if (cl$verbose > 0) { print(loglik) }
    
    if (loglik > best_loglik) {
      best_loglik = loglik
      best_res = results
    }
  }
  
  return(best_res)
  
}



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (fit_SBM_pair)
#' Fit a SBM given an input adjacency matrix (pair) and possibly starting model
#' 
#' @param adjm [matrix-int] :: Input adjacency matrix
#' @param adjm1 temp
#' @param adjm2 temp
#' @param Nobs [int] :: Number of observations
#' @param control_list [list] :: List of control parameters -- see \code{\link{set_fit_param}}
#' 
#' @return [\code{\link{NetworkModel}}] :: Best found fitted model
#' 
#' @export
#' 
fit_SBM_pair = function(adjm1, adjm2, Nobs = 1, control_list = set_fit_param()) { 
  
  ## TODO: Implement CV to help choose parameters... 
  ## TODO: Make this take multiple inputs? adjacency arrays, multiple observations? Need to alter further arguments to make this work! 
  
  # TODO: Implement bic?
  #   best_loglik = best_loglik - (Nclass*(Nclass-1)/2 + Nclass)*log(sum(!is.na(adjm))) ## (n + n(n-1)/2) [parameters] * log(n) [from BIC]
  #   best_k = which.max(best_loglik)
  
  ## 1-2: If desired
  ## 1. Do matrix completion
  ## -- This will only ever be done once, as the results are deterministic
  ## 2. Do spectral clustering
  ## -- can rerun this (since it gives a different warm start)
  ## 3. Throw stuff into EM algorithm
  
  ## For EM algo, use 'start' to decide how to initialize! If using spectral clustering to start, don't need to rerun the eigenvector computation part -- only need to rerun the k-means part! 
  ## start can take on values of 'spectral' or 'random'
  
  cl = control_list ## Shorten control_list variable name
  if (cl$SBM_start != "random") {
    ## (1) Matrix Completion (if needed)
    ## Some error checking / warnings
    
    if (cl$SBM_MC_Neigenvecs > 0 & cl$SBM_MC_laplacian) { print("WARNING: Using Laplacian and largest eigenvectors... is this intended?")}
    if (cl$SBM_MC_Neigenvecs < 0 & !cl$SBM_MC_laplacian) { print("WARNING: Using adjacnecy matrices and smallest eigenvectors... is this intended?")}
    if (cl$SBM_MC_Neigenvecs == 0) { stop("Must use a nonzero number of eigenvectors") }
    
    if (cl$SBM_start == "spectral-mean") {
      cmMethod = "rcmeans"
    } else if (cl$SBM_start == "spectral-complete") {
      cmMethod = "complete"
    }
    compEVS1 = completeMatrix(adjm = adjm1, method = cmMethod, laplacian = cl$SBM_MC_laplacian, eigenvecs = cl$SBM_MC_Neigenvecs, softImpute_maxit = cl$SBM_MC_softImpute_maxit, softImpute_thresh = cl$SBM_MC_softImpute_thresh, softImpute_rankmax = cl$SBM_MC_softImpute_rankmax)
    compEVS2 = completeMatrix(adjm = adjm2, method = cmMethod, laplacian = cl$SBM_MC_laplacian, eigenvecs = cl$SBM_MC_Neigenvecs, softImpute_maxit = cl$SBM_MC_softImpute_maxit, softImpute_thresh = cl$SBM_MC_softImpute_thresh, softImpute_rankmax = cl$SBM_MC_softImpute_rankmax)
  }
  
  ## (2) Spectral Clustering and (3) EM algorithm
  for(j in 1:cl$Ntries) {
    N = nrow(adjm1) ## This is the number of nodes
    
    best_loglik = -Inf
    best_res = NULL
    if (cl$SBM_start == 'random') {
      ## Random start
      nodeps = rep(1/cl$SBM_Nclass, length = cl$SBM_Nclass)
      edgeps1 = symmetrize_mat(matrix(sum(adjm, na.rm = TRUE) / (Nobs * N * (N-1) / 2) * runif(n = cl$SBM_Nclass * cl$SBM_Nclass, min = 0.1, max = 0.9), nrow = cl$SBM_Nclass))
      edgeps2 = symmetrize_mat(matrix(sum(adjm, na.rm = TRUE) / (Nobs * N * (N-1) / 2) * runif(n = cl$SBM_Nclass * cl$SBM_Nclass, min = 0.1, max = 0.9), nrow = cl$SBM_Nclass))
      
      H = matrix(0, nrow = N, ncol = cl$SBM_Nclass)
      PHI = matrix(runif(cl$SBM_Nclass*N, min = 0.1, max = 0.9), nrow = N)
      PHI = PHI / rowSums(PHI)
    } else if (length(grep("spectral", cl$SBM_start)) > 0) {
      
      ## Do spectral clustering once to give a rough start
      res = specClust(cbind(compEVS1, compEVS2), Nclass = cl$SBM_SC_Nclass, NStart = cl$SBM_SC_Nstart)
      res1 = fitModel(res, adjm1)
      clusts = res1@groups
      
      nodeps = sapply(1:cl$SBM_Nclass, function(x) {mean(clusts == x)})
      edgeps1 = res1@probmat
      
      res2 = fitModel(res, adjm2)
      edgeps2 = res2@probmat
      
      H = matrix(0, nrow = N, ncol = cl$SBM_Nclass)
      PHI = matrix(runif(cl$SBM_Nclass*N, min = 0.1, max = 0.9), nrow = N)
      for(j in 1:cl$SBM_Nclass) {
        PHI[which(j == clusts),j] = PHI[which(j == clusts),j] + 2
      }
      PHI = PHI / rowSums(PHI)
    }
    
    edgeps1[edgeps1 < 0.01] = .01; edgeps1[edgeps1 > 0.99] = .99
    
    results = EM_SBM_mf_twoedgeps_nonode(adjm1, adjm2, nodeps = nodeps, edgeps1 = edgeps1, edgeps2 = edgeps2, H = H, PHI = PHI, Niter = cl$SBM_EM_Niter, stop_thres = cl$SBM_EM_stopthres, verbose = cl$verbose, Nobs = 1)
    
    m1 = results$model; m2 = results$model
    m1@probmat = edgeps1; m2@probmat = edgeps2
    
    loglik = computeLik(m1, adja = adjm1)$sum + computeLik(m2, adja = adjm2)$sum 
    if (cl$verbose > 0) { print(loglik) }
    
    if (loglik > best_loglik) {
      best_loglik = loglik
      best_res = results
    }
  }
  
  return(best_res)
}


#' Runs spectral clustering
#' 
#' This requires input of eigenvectors. The eigenvectors can come from the Laplacian matrix, however. 
#' 
#' @param evs [matrix-numeric] :: Matrix of eigenvectors, with the columns as individual eigenvectors
#' @param Nclass [int] :: Number of classes to return for spectral clustering. 
#' @param NStart [int] :: Number of random starts for the k-means algorithm
#' 
#' @return [\code{\link{NetworkStruct}}] :: Output clustering's network structure
#' 
#' @export
#' 
specClust = function(evs, Nclass, NStart = 3) {
  ## evs should be a matrix of eigenvectors, with the columns as individual eigenvectors
  evs = scale(evs)
  
  clusts = kmeans(evs, centers = Nclass, nstart = NStart)$cluster
  NetS = NetworkStruct(set_model_param(Nnodes = length(clusts), type = 'block', block_nclass = Nclass, block_assign = clusts))
  
  return(NetS)
}


#' Estimates missing entries of an input matrix
#' 
#' This function applies one of two matrix-completion techniques. The cheap version ('rcmeans') simply fills each missing value with the average of all cells in its row/column. This results in a matrix that is still a properly weighted undirected network. The other version ('softImpute') calls the function in the corresponding R package 'softImpute'. Essentially, it looks for entries to fill that matrix such that it preserves low rank. I think the resulting matrix is symmetric if the input is symmetric, but there's no guarantee that it's a proper adjacency matrix (in that the values can possibly be negative...)
#' 
#' There is an option to return the filled adjacency matrix, the matrix converted to Laplacian form, or even simply the first few eigenvectors of the matrix (or the Laplacian)
#' 
#' This function requires existence of the 'rARPACK' package, since that is used to perfrom eigenvector/eigenvalue computations efficiently. Also, if method is set to 'softImpute', the 'softImpute' package is required. 
#' 
#' @param softImpute_rank [int] :: If using the 'softImpute method, this sets the maximum rank of the low-rank approximation. 
#' @param adjm [matrix-numeric] :: Input matrix
#' @param method [character] :: Two possible values: 'rcmeans' or 'softImpute'
#' @param laplacian [logical] :: If TRUE, the results are based on the Laplacian matrix
#' @param eigenvecs [int] :: If nonzero, this function returns eigenvectors instead of a matrix. If this is positive, it returns the 'eigenvecs' largest eigenvectors; if negative, it returns the 'eigenvecs' smallest eigenvectors. 
#' @param softImpute_rankmax [int] :: The max rank of the svd approximation, when using softImpute
#' @param softImpute_thresh [numeric] :: Convergence threshold, when using softImpute
#' @param softImpute_maxit [int] :: Maximum number of iterations, when using softImpute
#' 
#' @return [matrix] :: The returned value is either a square matrix, or a matrix containing the first few eigenvectors. 
#' 
#' @export
#' 
completeMatrix = function(adjm, method = "rcmeans", laplacian = FALSE, eigenvecs = 0, softImpute_rankmax = 5, softImpute_thresh = 1e-05, softImpute_maxit = 100) {
  require(rARPACK)
  
  ## method -- values = "rcmeans", "softImpute"
  ## softImpute requires the softImpute package
  ## softImpute_rank = number of SVD vectors to use in the approximation for imputation (be precise in final documentation)
  
  ## Note, for a real symmetric matrix with positive eigenvectors, the SVD vectors and eigenvectors are the same! Use this!! 
  
  ## If eigenvecs != 0, then return the components of the first 'eigenvecs' eigenvectors instead of returning an adjacency matrix. If positive, it returns the largest eigenvectors; if negatie, it returns the smallest eigenvectors -- explain this better! 
  
  
  ## TODO: Allow more general adjacency arrays / lists?
  N = nrow(adjm)
  diag(adjm) = 0
  
  if (method == "rcmeans") {
    ## Compute all row-means (thesee are the same as column means)
    rowEdgeCount = rowSums(adjm, na.rm = TRUE)
    rowDyadCount = N - 1 - rowSums(is.na(adjm)) 
    
    ## Fill in all missing cells
    for(j in seq_len(N)) { 
      for (k in seq_len(N)) {
        if ((j > k) & is.na(adjm[j,k])) { 
          val = (rowEdgeCount[j] + rowEdgeCount[k]) / (rowDyadCount[j] + rowDyadCount[k])
          adjm[j,k] = val; adjm[k,j] = val
        }
      }
    }
  } else if (method == "softImpute") {
    require(softImpute)
    est_svd = softImpute::softImpute(adjm, rank.max = softImpute_rankmax, thresh = softImpute_thresh, maxit = softImpute_maxit)
    if (!laplacian & (eigenvecs > 0) & (eigenvecs <= softImpute_rankmax)) {
      ## Special case for speed -- If we want eigenvectors, we can just output results from the softImpute stage, instead of estimating the adjacency matrix and then computing the eigenvectors (SVD) again. 
      return(est_svd$u[,seq_len(eigenvecs)])
    }
    adjm = softImpute::complete(adjm, est_svd)
  } else {
    stop("Invalid method choice (function = completeMatrix)")
  }
  
  ## Return the Laplacian instead of so desired, after completing the adjacency matrix
  if (laplacian) { res = graph_laplacian(adjm) } else { res = adjm }
  if (eigenvecs != 0) {
    require(rARPACK)
    ## requires rARPACK
    if (eigenvecs > 0) {
      return(eigs_sym(res, k = eigenvecs, which = "LM")$vectors)
    } else {
      return(eigs_sym(res, k = abs(eigenvecs), which = "SM")$vectors)
    } 
  } else {
    return(res) 
  }
}

