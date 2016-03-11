## For testing fitModel
library(netcompLib)
library(abind)
library(faraway)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")

if (FALSE) { ## Case -- single network
  NetM = NetworkModel(set_model_param(Nnodes = 10))
  adja = sampleNetwork(NetM)
  by_node = TRUE
  by_group = TRUE
  na.rm = TRUE
}

if (FALSE) { ## Case -- correlated network
  NetM = NetworkModelPair(set_model_param(Nnodes = 50, pairtype = "correlated-null"))
}

## Test Code
NetS = NetworkStruct(set_model_param(Nnodes = 20, block_nclass = 1))
adjm = matrix(0, nrow = 20, ncol = 20)
computeLik(fitModel(NetS, adjm), adjm)
