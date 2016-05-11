library(netcompLib)
library(igraph)


## For running fit_SBM
if (FALSE) {
  nm = NetworkModel(set_model_param(Nnodes = 60, block_assign = rep(c(1,2), each = 30), block_probs = matrix(c(.4,.2,.2,.4), nrow = 2)))
  adjm = sampleNetwork(nm)[,,1]
  Nobs = 1
  adjm = hide_edges(adjm, frac = .5)
  control_list = set_fit_param(SBM_start = "spectral-mean", SBM_Nclass = 2)
  control_list = set_fit_param(SBM_start = "spectral-mean", SBM_Nclass = 2, SBM_MC_Neigenvecs = 4, SBM_MC_softImpute_maxit = 500, SBM_MC_softImpute_rankmax = 5, SBM_EM_mode = "no-node-default")
  
  adjm1 = sampleNetwork(nm)[,,1]
  adjm2 = sampleNetwork(nm)[,,1]
  
  NN = 60
  genSBM = NetworkModelSBM(
    model_params = set_model_param(
      Nnodes = NN, block_assign = rep(1:3, times = NN/3)[seq_len(NN)],
      block_probs = matrix(c(.8,.3,.1, .3,.6,.1, .1,.1,.5), nrow = 3)))
  adjm = sampleNetwork(genSBM)[,,1]
  control_list = set_fit_param(SBM_start = "spectral-mean", SBM_Nclass = 3, SBM_MC_Neigenvecs = 4)
  
  bestmod = fitModel(NetworkStruct(set_model_param(Nnodes = 30, block_assign = rep(c(1,2), each = 15))), adjm)
  computeLik(bestmod, adjm)$sum
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



adjm = hide_edges(adjm, frac = .5)


nm = NetworkModel(set_model_param(Nnodes = 40, block_assign = rep(c(1,2), each = 20)))
adjm = sampleNetwork(nm)[,,1]
z = fit_SBM(adjm, control_list = set_fit_param(SBM_EM_mode = "no-node-default"))
#|----##If calling the version from netcompLib, the input parameters have been changed! --Thu Jan 14 20:48:05 2016--
computeLik(z$model, adjm)$sum

EM_SBM_mf_default_nonode

nm1 = NetworkModel(set_model_param(Nnodes = 50, block_assign = rep(c(1,2), each = 25), block_probs = matrix(c(.5, .1, .1, .5), nrow = 2) ))
nm2 = NetworkModel(set_model_param(Nnodes = 50, block_assign = rep(c(1,2), each = 25), block_probs = matrix(c(.1, .5, .5, .1), nrow = 2) ))


## nulldist
adjm = sampleNetwork(nm1)[,,1]
fit_SBM(adjm, start = 'spectral', Nclass = 2)
#|----##If calling the version from netcompLib, the input parameters have been changed! --Thu Jan 14 20:48:05 2016--

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


specClust(hide_edges(adjm, frac = .1), 3)

adjm = hide_edges(adjm, frac = .2)




nm = NetworkModel(set_model_param(Nnodes = 30, block_assign = rep(c(1,2), each = 15), block_probs = matrix(c(.4,.2,.2,.4), nrow = 2)))
adjm = sampleNetwork(nm)[,,1]
Nobs = 1
adjm = hide_edges(adjm, frac = .5)
test = 
  fit_SBM(adjm, control_list = set_fit_param(SBM_start = "random", SBM_Nclass = 2, SBM_MC_Neigenvecs = 4, SBM_EM_mode = "no-node-default", SBM_SC_Nstart = 10))
