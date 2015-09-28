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

NetM = fitModel(NetS, adja, "corr-global-null")





C = rbind(c(50, 30, 30, 20), c(60, 40, 40, 20))
n = c(130, 160)
llFx_cnull(c(0,0,0), C, n)
llFx_calt(c(0,0,0,0,0), C, n)

optim(par = c(0,0,0,0,0), fn = llFx_calt, gr = llGrFx_calt, method = "BFGS", control = list(fnscale = "-1"), C = C, n = n)
optim(par = c(0,0,1), fn = llFx_cnull, gr = llGrFx_cnull, method = "BFGS", control = list(fnscale = "-1"), C = C, n = n)




