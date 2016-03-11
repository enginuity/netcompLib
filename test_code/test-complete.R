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
                1,0,0,1,0,
                0,0,1,0,1,
                0,1,0,1,0), nrow = 5)

NetS = NetworkStruct(set_model_param(Nnodes = 5, block_assign = c(1,1,1,2,2)))
m1 = fitModel(NetS, adj1)
m2 = fitModel(NetS, adj2)

abs(m1@probmat - matrix(c(2/3, 3/6, 3/6, 1/1), nrow = 2))
abs(m2@probmat - matrix(c(1/3, 2/6, 2/6, 1/1), nrow = 2))
mc = fitModel(NetS, abind(adj1, adj2, along = 3))
abs(mc@probmat - matrix(c(3/6, 5/12, 5/12, 2/2), nrow = 2))





l1 = computeLik(m1, adj1, by_node = TRUE, by_group = TRUE)
computePval(NetS, adj1, adj2, pl = set_sim_param(cc_adj = 0, thres_ignore = 0), output_mode = "chisq")




