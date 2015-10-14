library(netcompLib)
library(igraph)


## For running fit_SBM
if (FALSE) {
  nm = NetworkModel(set_model_param(Nnodes = 100, block_assign = rep(c(1,2), each = 50)))
  adjm = sampleNetwork(nm)[,,1]
  Nobs = 1
  Nclass = 3
  Niter = 100
  Ntries = 10
} 

## For running EM mean field version
if (FALSE) {
  N = nrow(adjm) ## This is the number of nodes
  nodeps = rep(1/Nclass, length = Nclass)
  edgeps = symmetrize_mat(matrix(sum(adjm, na.rm = TRUE) / (Nobs * N * (N-1) / 2) * runif(n = Nclass * Nclass, min = 0.1, max = 0.9), nrow = Nclass))
    
  H = matrix(0, nrow = N, ncol = Nclass)
  PHI = matrix(runif(Nclass*N, min = 0.1, max = 0.9), nrow = N)
  PHI = PHI / rowSums(PHI)
  
  stop_thres = 0.00005
  verbose = TRUE
}


nm = NetworkModel(set_model_param(Nnodes = 100, block_assign = rep(c(1,2), each = 50)))
adjm = sampleNetwork(nm)[,,1]
z = fit_SBM(adjm)
computeLik(z$model, adjm)$sum




