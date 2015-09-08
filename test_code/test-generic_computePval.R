## For testing computePval 
library(netcompLib)

## For setting parameters
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
NetM = NetworkModel(set_model_param(Nnodes = 30, type = "block"))
adja1 = sampleNetwork(NetM)
adja2 = sampleNetwork(NetM)
Nobs = 1
pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_structs = c(1,20))

## Pick the appropriate parameter
NetS = NetworkStructRND(set_model_param(Nnodes = 30, block_nclass = 3))
NetS = NetworkStructHRG(set_model_param(Nnodes = 30, block_nclass = 3))
NetS = NetworkStructSBM(set_model_param(Nnodes = 30, block_nclass = 3))
