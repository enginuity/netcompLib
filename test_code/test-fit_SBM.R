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



nm1 = NetworkModel(set_model_param(Nnodes = 50, block_assign = rep(c(1,2), each = 25), block_probs = matrix(c(.5, .1, .1, .5), nrow = 2) ))
nm2 = NetworkModel(set_model_param(Nnodes = 50, block_assign = rep(c(1,2), each = 25), block_probs = matrix(c(.1, .5, .5, .1), nrow = 2) ))


## nulldist
adjm = sampleNetwork(nm1)[,,1]
fit_SBM(adjm, start = 'spectral', Nclass = 2)

plot_probmatrix(pm = adjm, add_legend = FALSE, colors = c(0,1))



NETSIZE = 30
GL = NetworkModelPair(
  m1 = set_model_param(type = "block", Nnodes = NETSIZE, 
                       block_assign = rep(1:3, each = ceiling(NETSIZE/3))[1:NETSIZE], 
                       block_probs = matrix(c(.5, .1, .1,
                                              .1, .5, .2,
                                              .1, .2, .5), nrow = 3)), 
  m2 = set_model_param(type = 'block', Nnodes = NETSIZE, 
                       block_assign = rep(1:3, each = ceiling(NETSIZE/3))[1:NETSIZE], 
                       block_probs = matrix(c(.1, .5, .1,
                                              .5, .1, .2,
                                              .1, .2, .5), nrow = 3)), 
  is_null = FALSE)
adjm = sampleNetwork(GL@m1)[,,1]
