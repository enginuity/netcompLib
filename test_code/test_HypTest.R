
library(netcompLib)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
library(abind)

test9 = sim_hyptest(gen_NetMPair = NetworkModelPair(NetworkModelSBM(), NetworkModelSBM(), is_null = FALSE), fitm_params = set_model_param(30, 'block'), Nobs = 1, Nsim = 100, pl = set_sim_param(n_structs = c(1, 10, 50, 100), fitstruct_method = 'random'), verbose = TRUE, vbset = c(1,1,0))
#|----##Parameters changed --Thu Sep 17 05:27:59 2015--

