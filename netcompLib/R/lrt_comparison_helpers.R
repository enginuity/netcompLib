##@S This is a collection of functions that aid with network comparisons.
##@S These are NOT generic functions, but they DO use NetM/NetS's as arguments. 

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computeDfAdj)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computeDfAdj)
#' Function to compute the theoretical df adjustment
#' 
#' @param NetM Given NetworkModel object
#' @param NetS Given NetworkStruct object
#' @param hidden_nodes ewefwf
#' 
#' @return temp
#' 
#' @export
#' 
computeDfAdj = function(NetM, NetS, hidden_nodes = NULL) {
  # hidden_nodes should be a logical matrix, with TRUE identifying edges that are hidden.
  ## TODO: [Dumb Variable Name] -- rename hidden_nodes with hidden_edges !?!??!?!
  if (FALSE) {
    ## Test code
    NetM = NetworkModel(Nnodes = 20, type = "block", model_param = set_model_param(block_assign = rep(c(1,2), each = 10), block_probs = matrix(c(.3, .3, .3, .7), nrow = 2)))
    NetS = NetworkStructSBM(Nnodes = 20, model_param = set_model_param(block_assign = rep(c(1,2), times = 10)))
  }
  
  if (getNnodes(NetM) != getNnodes(NetS)) { stop("Number of nodes is not consistent between input NetM and NetS") }
  test_mod = fitModel(NetS, adja = NULL)
  test_nodeids = getEdgeProbMat(NetM = test_mod, mode = "group")
  fit_nodeprobs = getEdgeProbMat(NetM, mode = "prob")
  test_nodeids[lower.tri(test_nodeids, diag = TRUE)] = 0
  
  if (!is.null(hidden_nodes)) { test_nodeids[hidden_nodes] = 0 }
  
  ids = unique(as.vector(test_nodeids))
  ids = ids[ids > 0]
  
  cors = 0 * seq_along(ids)
  for(j in seq_along(ids)) {
    matches = (test_nodeids == ids[j])
    probs = fit_nodeprobs[matches]
    avprob = mean(probs)
    cors[j] = (mean(probs^2) - avprob^2)  /  (avprob - avprob^2)
  }
  return(1 - cors)
}

