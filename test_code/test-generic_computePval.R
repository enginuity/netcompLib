## For testing computePval 

library(netcompLib)
library(abind)
library(faraway)
library(microbenchmark)

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

computePval(NetS = NetS, adja1, adja2, pl = pl, output_mode = 'chisq', model_type = 'default')
computePval(NetS = NetS, adja1, adja2, pl = pl, output_mode = 'chisq', model_type = 'default-slow')
computePval(NetS = NetS, ta, te, pl = pl, output_mode = 'chisq', model_type = 'default')
computePval(NetS = NetS, ta, te, pl = pl, output_mode = 'chisq', model_type = 'default-slow')

microbenchmark(computePval(NetS = NetS, ta, te, pl = pl, output_mode = 'chisq', model_type = 'default'), 
               computePval(NetS = NetS, ta, te, pl = pl, output_mode = 'chisq', model_type = 'default-slow'))



