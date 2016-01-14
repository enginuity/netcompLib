##@S Testing all the functions that work with HRGs. 

library(microbenchmark)
library(netcompLib)
net1 = NetworkModel(set_model_param(Nnodes = 10, type = "tree"))

HRG_depthFromParents(net1@parents)
HRG_treeFromParents(net1@parents)

## Old tree format... 
oldtree = list(nodes = net1@Nnodes, parents = net1@parents, children = net1@children)

HRG_expandedChildren(net1)
HRG_closestAncestor(net1)$anc_table
