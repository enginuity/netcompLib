##@S Testing all the functions that work with HRGs. 

library(netcompLib)
net1 = NetworkModel(set_model_param(Nnodes = 10, type = "tree"))

depth_from_parents(net1@parents)
