## For testing fitModel
library(netcompLib)
library(abind)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")


##
NetM = NetworkModel(set_model_param(Nnodes = 6))
adja1 = sampleNetwork(NetM)
adja2 = sampleNetwork(NetM)

getEdgeProbMat(NetM, "group")
getEdgeProbMat(NetM, "prob")

NetM@probmat

aggstat_dendiff(NetM, adja1, adja2)
aggstat_corr(NetM, adja1, adja2)

aggstat_single(NetworkModelHRG(set_model_param(Nnodes = 6)), adja1)
aggstat_single(NetworkModelRND(set_model_param(Nnodes = 6)), adja1)


aggstat_single(NetM, getEdgeProbMat(NetM))
getEdgeProbMat(NetM)
