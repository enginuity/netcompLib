##@S Testing all the functions that work with HRGs. 

library(microbenchmark)
library(netcompLib)
net1 = NetworkModel(set_model_param(Nnodes = 10, type = "tree"))

depth_from_parents(net1@parents)
tree_from_parents(net1@parents)


oldtree = list(nodes = net1@Nnodes, parents = net1@parents, children = net1@children)
expanded_children_from_tree(net1)

closest_ancestor_new(net1)$anc_table
