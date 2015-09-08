# testing new functionality

library(netcompLib)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")


test1 = NetworkModel(Nnodes = 10)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
getNnodes(test1)
getNetType(test1)

test2 = NetworkModel(Nnodes = 10, type = "block", model_param = set_model_param(block_assign = rep(c(1,2), each = 5), block_probs = matrix(c(0.1, 0.7, 0.7, 0.5), nrow = 2)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
getNnodes(test2)
getNetType(test2)
getEdgeProbMat(test2)
sampleNetwork(test2)
extractStruct(test2)

test2@assign
test2@probmat

test3 = NetworkModel(Nnodes = 10, type = "tree")
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
getNnodes(test3)
getNetType(test3)
test3
getEdgeProbMat(test3)
extractStruct(test3)

test4 = NetworkModel(Nnodes = 10, type = "latent")
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
getNnodes(test4)
getNetType(test4)
round(getEdgeProbMat(test4),2)

test5 = NetworkModel(Nnodes = 5, type = "random", model_param = set_model_param(random_ngroups = 3))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
getNnodes(test5)
getNetType(test5)
round(getEdgeProbMat(test5),2)
extractStruct(test5)

NetworkStruct(Nnodes = 15, type = "block")
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:31 2015--
NetworkStruct(Nnodes = 15, type = "tree")
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:31 2015--
NetworkStruct(Nnodes = 15, type = "random")
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:31 2015--


pl = list(cc_adj = c(1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_models = c(1,20))
#|----##In set_sim_param, n_models is renamed into n_structs --Tue Sep  8 02:23:23 2015--
computePval(extractStruct(test2), sampleNetwork(test2), sampleNetwork(test2), 1, pl)
computePval(extractStruct(test2), sampleNetwork(test2), sampleNetwork(test2), 1, pl)
computePval(extractStruct(test2), sampleNetwork(test2), sampleNetwork(test2), 1, pl)
computePval(extractStruct(test2), sampleNetwork(test2), sampleNetwork(test2), 1, pl)

netsl = NetworkStructList(type = "block")
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:32 2015--
computePval(netsl, sampleNetwork(test2), sampleNetwork(test2), 1, pl)
getNetType(netsl)

NetM = NetworkModel(Nnodes = 30, type = "block")
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
adja1 = sampleNetwork(NetM)
adja2 = sampleNetwork(NetM)
Nobs = 1

NetS = NetworkStructSBM(Nnodes = 30, model_param = set_model_param(block_nclass = 3))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:33 2015--
getNetType(NetS)

m1 = NetworkModel(Nnodes = 50, type = "block", model_param = set_model_param(block_probs = matrix(0.45, nrow = 1, ncol = 1), block_assign = rep(1, times = 50)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
m2 = NetworkModel(Nnodes = 50, type = "block", model_param = set_model_param(block_probs = matrix(0.5, nrow = 1, ncol = 1), block_assign = rep(1, times = 50)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--

computePval(NetworkStructList(Nnodes = 50, type = "tree"), sampleNetwork(m1), sampleNetwork(m2), 1, pl)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:32 2015--
computePval(NetworkStructList(Nnodes = 50, type = "random"), sampleNetwork(m1), sampleNetwork(m2), 1, pl)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:32 2015--


nmp = NetworkModelPair(m1 = NetworkModel(type = "block"), is_null = TRUE)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:28 2015--
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
sampleNetwork(NetM = nmp)


