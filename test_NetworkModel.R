# testing new functionality

library(netcompLib)

test1 = NetworkModel(Nnodes = 10)
getNnodes(test1)
getNetType(test1)

test2 = NetworkModel(Nnodes = 10, type = "block", model_param = set_model_param(block_assign = rep(c(1,2), each = 5), block_probs = matrix(c(0.1, 0.7, 0.7, 0.5), nrow = 2)))
getNnodes(test2)
getNetType(test2)
getEdgeProbMat(test2)
sampleNetwork(test2)
extractStruct(test2)

test2@assign
test2@probmat

test3 = NetworkModel(Nnodes = 10, type = "tree")
getNnodes(test3)
getNetType(test3)
test3
getEdgeProbMat(test3)
extractStruct(test3)

test4 = NetworkModel(Nnodes = 10, type = "latent")
getNnodes(test4)
getNetType(test4)
round(getEdgeProbMat(test4),2)

test5 = NetworkModel(Nnodes = 5, type = "random", model_param = set_model_param(random_ngroups = 3))
getNnodes(test5)
getNetType(test5)
round(getEdgeProbMat(test5),2)
extractStruct(test5)

NetworkStruct(Nnodes = 15, type = "block")
NetworkStruct(Nnodes = 15, type = "tree")
NetworkStruct(Nnodes = 15, type = "random")
