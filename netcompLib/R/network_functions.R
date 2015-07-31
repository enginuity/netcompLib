##@S This is a collection of functions that aid with network comparisons.
##@S These are NOT generic functions, but they DO use NetM/NetS's as arguments. 


#' Compute Theoretical df Adjustment
#' 
#' For a fixed null hypothesis model, and for a fixed network structure, compute the true theoretical df adjustment
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Underlying null model
#' @param NetS [\code{\link{NetworkStruct}}] :: Testing model structure
#' @param hidden_edges [matrix-logical] :: Entries of matrix are 'true' if they are hidden in the inference process. (so default would be a matrix of falses)
#' 
#' @return A vector of theoretical df adjustments (one for each edge group)
#' 
#' @export
#' 
computeDfAdj = function(NetM, NetS, hidden_edges = NULL) {
  if (FALSE) {
    ## Test code
    NetM = NetworkModel(set_model_param(Nnodes = 20, type = "block", block_assign = rep(c(1,2), each = 10), block_probs = matrix(c(.3, .3, .3, .7), nrow = 2)))
    NetS = NetworkStructSBM(Nnodes = 20, model_param = set_model_param(block_assign = rep(c(1,2), times = 10)))
  }
  
  if (getNnodes(NetM) != getNnodes(NetS)) { stop("Number of nodes is not consistent between input NetM and NetS") }
  test_mod = fitModel(NetS, adja = NULL)
  test_nodeids = getEdgeProbMat(NetM = test_mod, mode = "group")
  fit_nodeprobs = getEdgeProbMat(NetM, mode = "prob")
  test_nodeids[lower.tri(test_nodeids, diag = TRUE)] = 0
  
  if (!is.null(hidden_edges)) { test_nodeids[hidden_edges] = 0 }
  
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



#' Computes a distance between two network models
#' 
#' @param x [\code{\link{NetworkModel}}] :: First network model
#' @param y [\code{\link{NetworkModel}}] :: Second network model
#' @param type [char; ALLOWED = c("KL", "KLsym")] :: What distance measure to use? \cr
#' \itemize{
#'   \item KL -- KL distance 
#'   \item KLsym -- symmetric KL distance (average of KL distances in both ways)
#' }
#' 
#' @return [double] :: distance value
#' 
#' @export
#' 
computeDist = function(x, y, type = "KLsym") {
  p1 = getEdgeProbMat(x); p2 = getEdgeProbMat(y)
  
  h_kldist = function(p1, p2) {
    diag(p1) <- 0.5; diag(p2) <- 0.5
    q1 = 1 - p1; q2 = 1 - p2
    res = log(q1) - log(q2) + p1 * (log( (p1 / q1) / (p2 / q2 ) ) )
    diag(res) <- 0
    return(sum(res))
  }
  
  if (type == "KL") { 
    return(h_kldist(p1, p2))
  } else if (type == "KLsym") {
    return( (h_kldist(p1,p2)+h_kl_dist(p2,p1)) / 2)
  }
  stop("Invalid type of distance measure inputted")
  return(NULL)
}
