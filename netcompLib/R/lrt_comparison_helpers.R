##@S This is a collection of functions that aid with network comparisons.
##@S These are NOT generic functions, but they DO use NetM/NetS's as arguments. 

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computeDfAdj)
#' Function to compute the theoretical df adjustment
#' 
#' @param NetM temp
#' @param NetS temp
#' 
#' @return temp
#' 
#' @export
#' 
computeDfAdj = function(NetM, NetS) {
  if (FALSE) {
    ## Test code
    NetM = NetworkModel(Nnodes = 50, type = "block")
    NetS = NetworkStructSBM(Nnodes = 50)
  }
  
  if (getNnodes(NetM) != getNnodes(NetS)) { stop("Number of nodes is not consistent between input NetM and NetS") }
  test_mod = fitModel(NetS, adja = NULL)
  test_nodeids = getEdgeProbMat(NetM = test_mod, mode = "group")
  fit_nodeprobs = getEdgeProbMat(NetM, mode = "prob")
  test_nodeids[lower.tri(test_nodeids, diag = TRUE)] = 0
  
  ids = unique(as.vector(test_nodeids))
  ids = ids[ids > 0]
  
  cors = 0 * seq_along(ids)
  for(j in seq_along(ids)) {
    matches = test_nodeids == ids[j]
    probs = fit_nodeprobs[matches]
    avprob = mean(probs)
    count = length(probs)
    cors[j] = (1/count * sum(probs^2) - avprob^2)  /  (avprob - avprob^2)
  }
  return(1 - cors)
}

