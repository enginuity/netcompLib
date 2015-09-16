## For testing fitModel
library(netcompLib)
library(abind)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")


##
NetM = NetworkModel()
adja1 = sampleNetwork(NetM)
adja2 = sampleNetwork(NetM)


