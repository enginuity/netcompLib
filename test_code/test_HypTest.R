
library(netcompLib)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
library(abind)
system.time(test9 <- sim_hyptest(gen_NetMPair = NetworkModelPair(NetworkModelSBM(100), NetworkModelSBM(100), is_null = FALSE), fit_NetSList = NetworkStructList(100, 100, type = "block"), Nobs = 1, Nsim = 100, param_list = set_sim_param(n_models = c(1, 10, 100)), pval_adj_fx = list(mult_bonferroni, mult_pearson, mult_highcrit), verbose = TRUE)) ## not much slower! 

system.time(test8 <- sim_test_v2(Nnodes = 100, Nsim = 100, Nobs = 1, gen_null = FALSE, gen_mode = "block", fit_mode = "block", param_list = set_sim_param(n_models = c(1, 10, 100)), complete = FALSE, pval_adj_fx = list(mult_bonferroni, mult_pearson, mult_highcrit), verbose = TRUE))