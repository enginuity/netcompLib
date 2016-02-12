## For testing computePval 
library(netcompLib)
library(abind)
library(faraway)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")

## For setting parameters
NetM = NetworkModel(set_model_param(Nnodes = 30, type = "block"))
adja1 = sampleNetwork(NetM)
adja2 = sampleNetwork(NetM)
Nobs = 1
pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_structs = c(1,20))


te = hide_edges(adja1[,,1], frac = 0.5)
ta = adja2[,,1]
ta[is.na(te)] = NA
## Pick the appropriate parameter
NetS = NetworkStructRND(set_model_param(Nnodes = 30, block_nclass = 3))
NetS = NetworkStructHRG(set_model_param(Nnodes = 30, block_nclass = 3))
NetS = NetworkStructSBM(set_model_param(Nnodes = 30, block_nclass = 3))


computePval(NetS, ta, te, pl = pl, mode = 'chisq')
#|----##Changed parameter 'mode' to 'output_mode' --Fri Feb 12 15:17:37 2016--
computePval(NetworkStructList(Nmodels = 10, set_model_param()), ta, te, pl = pl, mode = 'chisq')
#|----##Changed parameter 'mode' to 'output_mode' --Fri Feb 12 15:17:37 2016--
