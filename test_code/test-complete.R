## Code to test the simple hypothesis testing, but check all intermediate outputs. 
library(netcompLib)
library(abind)

## 
adj1 = matrix(c(0,1,1,0,0,
                1,0,0,1,1,
                1,0,0,1,0,
                0,1,1,0,1,
                0,1,0,1,0), nrow = 5)
adj2 = matrix(c(0,0,1,0,0,
                0,0,0,0,1,
                1,0,0,1,1,
                0,0,1,0,1,
                0,1,1,1,0), nrow = 5)

NetS = NetworkStruct(set_model_param(Nnodes = 5, block_assign = c(1,1,1,2,2)))
m1 = fitModel(NetS, adj1)
m2 = fitModel(NetS, adj2)
mc = fitModel(NetS, abind(adj1, adj2, along = 3))
computePval(NetS, adj1, adj2, pl = set_sim_param(cc_adj = 0, thres_ignore = 0))
