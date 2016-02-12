# testing new functionality

library(netcompLib)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")


n1 = NetworkModel(Nnodes = 30, type = "block", model_param = set_model_param(block_assign = rep(c(1,2,3), each = 10), block_probs = matrix(c(.4, .1, .6, .1, .7, .4, .6, .4, .5), nrow = 3)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
n2 = NetworkModel(Nnodes = 30, type = "block", model_param = set_model_param(block_assign = rep(c(1,2,3), each = 10), block_probs = matrix(c(.4, .1, .6, .1, .3, .2, .6, .2, .5), nrow = 3)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--

m1 = sampleNetwork(n1, Nobs = 1)
m2 = sampleNetwork(n2, Nobs = 1)

computePval(NetworkStructSBM(Nnodes = 30), adja1 = m1, adja2 = m2, Nobs = 1, pl = set_sim_param(n_models = 1), mode = "nodewise")
#|----##Changed parameter 'mode' to 'output_mode' --Fri Feb 12 15:17:37 2016--
#|----##In set_sim_param, n_models is renamed into n_structs --Tue Sep  8 02:23:23 2015--
#|----##Function modified -- all calls need to be updated.. --Wed Sep  2 20:40:58 2015--
#|----##Function modified -- all calls need to be updated.. --Tue Aug  4 00:51:48 2015--
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:33 2015--

netsl = NetworkStructList(Nnodes = 30, type = "block", Nmodels = 100)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:32 2015--
res = computePval(netsl, m1, m2, 1, pl = set_sim_param(n_models = c(100)), mode = "nodewise")
#|----##Changed parameter 'mode' to 'output_mode' --Fri Feb 12 15:17:37 2016--
#|----##In set_sim_param, n_models is renamed into n_structs --Tue Sep  8 02:23:23 2015--
#|----##Function modified -- all calls need to be updated.. --Wed Sep  2 20:40:58 2015--
#|----##Function modified -- all calls need to be updated.. --Tue Aug  4 00:51:48 2015--

par(mfrow = c(2,2))
temp = sapply(res, function(x) x$node)
plot(1:30, apply(temp, 1, function(x) { mean(x[is.finite(x)])}), ylab = "Mean Nodal Contribution", xlab = "Node", main = "Network Size = 30", ylim = c(-5, 10))
abline(v = c(10.5, 20.5))

y = sapply(list(1:10, 11:20, 21:30), function(m) {mean(apply(temp, 1, function(x) { mean(x[is.finite(x)])})[m])})
segments(c(1,11,21), y, c(10, 20, 30), y, lty = 2)


plot(1:30, apply(m1 - m2 != 0, 1, sum), ylab = "Number of Inconsistent Edges", xlab = "Nodes", main = "Network Size = 30")
abline(v = c(10.5, 20.5))
y = sapply(list(1:10, 11:20, 21:30), function(m) {mean(apply(m1 - m2 != 0, 1, sum)[m])})
segments(c(1,11,21), y, c(10, 20, 30), y, lty = 2)


n1 = NetworkModel(Nnodes = 60, type = "block", model_param = set_model_param(block_assign = rep(c(1,2,3), each = 20), block_probs = matrix(c(.4, .1, .6, .1, .7, .4, .6, .4, .5), nrow = 3)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--
n2 = NetworkModel(Nnodes = 60, type = "block", model_param = set_model_param(block_assign = rep(c(1,2,3), each = 20), block_probs = matrix(c(.4, .1, .6, .1, .3, .2, .6, .2, .5), nrow = 3)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:00:51 2015--

m1 = sampleNetwork(n1, Nobs = 1)
m2 = sampleNetwork(n2, Nobs = 1)

computePval(NetworkStructSBM(Nnodes = 60), adja1 = m1, adja2 = m2, Nobs = 1, pl = set_sim_param(n_models = 1), mode = "nodewise")
#|----##Changed parameter 'mode' to 'output_mode' --Fri Feb 12 15:17:37 2016--
#|----##In set_sim_param, n_models is renamed into n_structs --Tue Sep  8 02:23:23 2015--
#|----##Function modified -- all calls need to be updated.. --Wed Sep  2 20:40:58 2015--
#|----##Function modified -- all calls need to be updated.. --Tue Aug  4 00:51:48 2015--
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:33 2015--

netsl = NetworkStructList(Nnodes = 60, type = "block", Nmodels = 100)
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:32 2015--
res = computePval(netsl, m1, m2, 1, pl = set_sim_param(n_models = c(100)), mode = "nodewise")
#|----##Changed parameter 'mode' to 'output_mode' --Fri Feb 12 15:17:37 2016--
#|----##In set_sim_param, n_models is renamed into n_structs --Tue Sep  8 02:23:23 2015--
#|----##Function modified -- all calls need to be updated.. --Wed Sep  2 20:40:58 2015--
#|----##Function modified -- all calls need to be updated.. --Tue Aug  4 00:51:48 2015--

temp = sapply(res, function(x) x$node)
plot(1:60, apply(temp, 1, function(x) { mean(x[is.finite(x)])}), ylab = "Mean Nodal Contribution", xlab = "Node", main = "Network Size = 60", ylim = c(-5, 10))
abline(v = c(20.5, 40.5))
y = sapply(list(1:20, 31:40, 41:60), function(m) {mean(apply(temp, 1, function(x) { mean(x[is.finite(x)])})[m])})
segments(c(1,21,41), y, c(20, 40, 60), y, lty = 2)


plot(1:60, apply(m1 - m2 != 0, 1, sum), ylab = "Number of Inconsistent Edges", xlab = "Nodes", main = "Network Size = 60")
abline(v = c(20.5, 40.5))
y = sapply(list(1:20, 31:40, 41:60), function(m) {mean(apply(m1 - m2 != 0, 1, sum)[m])})
segments(c(1,21,41), y, c(20, 40, 60), y, lty = 2)

