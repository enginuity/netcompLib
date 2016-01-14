
library(netcompLib)

net1 = NetworkModel(set_model_param(type = 'block', Nnodes = 10, block_nclass = 4))
getEdgeProbMat(net1)
getEdgeProbMat(net1, mode = 'group')

net2 = NetworkModel(set_model_param(type = 'tree', Nnodes = 10))
getEdgeProbMat(net2)
getEdgeProbMat(net2, mode = 'group')

net3 = NetworkModel(set_model_param(type = 'random', Nnodes = 10))
getEdgeProbMat(net3)
getEdgeProbMat(net3, mode = 'group')
