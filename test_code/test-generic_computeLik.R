## For testing fitModel
library(netcompLib)
library(abind)
library(faraway)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")

if (FALSE) {
  NetM = NetworkModel(set_model_param(Nnodes = 50))
  adja = sampleNetwork(NetM)
  by_node = TRUE
  by_group = TRUE
  na.rm = TRUE
}



