## For testing fitModel
library(netcompLib)
library(abind)
library(faraway)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")


## For SBM version
NetS = NetworkStruct(set_model_param(Nnodes = 100))
genmod = NetworkModelSBM(set_model_param(Nnodes = 100))
getNnodes(genmod)
adja = sampleNetwork(genmod, Nobs = 2)
optim_tries = 10
mode = 'densitydiff'


getEdgeProbMat(genmod, 'group')



