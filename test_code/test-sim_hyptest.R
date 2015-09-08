## For testing sim_hyptest functions
library(netcompLib)
library(abind)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")


## For sim_hyptest
gen_NetMPair = NetworkModelPair(set_model_param(), set_model_param(), is_null = FALSE)
fit_NetSList = NetworkStructList(Nmodels = 100)
Nobs = 1; Nsim = 100; verbose = TRUE;
pl = set_sim_param(cc_adj = c(0,1,2), thres_ignore = c(5, 10), n_structs = c(1, 25, 50, 100))


## For sim_power_rpart
GL = list(NetworkModelPair(m1 = NetworkModelSBM(set_model_param(Nnodes = 30)), is_null = TRUE), 
          NetworkModelPair(m1 = NetworkModelSBM(set_model_param(Nnodes = 30)), is_null = TRUE))
NL = GL
FL = list(set_model_param(Nnodes = 30), set_model_param(Nnodes = 30))
Nsim = 20
Nsim_crit = 20
Nobs = 1
verbose = 5
pl = set_sim_param()
