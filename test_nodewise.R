# testing new functionality

library(netcompLib)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")


n1 = NetworkModel(Nnodes = 30, type = "block", model_param = set_model_param(block_assign = rep(c(1,2,3), each = 10), block_probs = matrix(c(.4, .1, .6, .1, .7, .4, .6, .4, .5), nrow = 3)))
n2 = NetworkModel(Nnodes = 30, type = "block", model_param = set_model_param(block_assign = rep(c(1,2,3), each = 10), block_probs = matrix(c(.4, .1, .6, .1, .3, .2, .6, .2, .5), nrow = 3)))

m1 = sampleNetwork(n1, Nobs = 1)
m2 = sampleNetwork(n2, Nobs = 1)

computePval(NetworkStructSBM(Nnodes = 30), adja1 = m1, adja2 = m2, Nobs = 1, pl = set_sim_param(n_models = 1), mode = "nodewise")

netsl = NetworkStructList(Nnodes = 30, type = "block", Nmodels = 10)

computePval(netsl, m1, m2, 1, pl = set_sim_param(n_models = c(1, 10)), mode = "nodewise")


nmp = NetworkModelPair(m1 = NetworkModel(type = "block"), is_null = TRUE)
sampleNetwork(NetM = nmp)

