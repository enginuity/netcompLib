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
#' @return [vector-double] :: Theoretical df adjustments for each edge group
#' 
#' @export
#' 
computeTrueDfAdj = function(NetM, NetS, hidden_edges = NULL) {
  if (getNnodes(NetM) != getNnodes(NetS)) { stop("Number of nodes is not consistent between input NetM and NetS") }
  
  ## TODO: [Issue #16] Need to handle case when it's not standard hypothesis test...
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


#' Estimates the degrees of freedom empirically
#' 
#' @param adjm1 [matrix-int] :: Input adjacency matrix
#' @param adjm2 [matrix-int] :: Input adjacency matrix
#' @param NetS [\code{\link{NetworkStruct}}] :: Structure used for testing
#' @param model_type [char] :: Model type to compute degrees of freedom for. See ## FIX THIS## for description of model_type. 
#' 
#' @return [vector-double] :: Degrees of freedom for each parameter in NetS
#' 
#' @export
#' 
computeEmpDfAdj = function(adjm1, adjm2, NetS, model_type = "default") {
  ## Computes empirical df with respect to a specific network structure
  ## TODO: This only works when the input is a single pair of matrices, and also with no missing values. Need to check what happens when there are missing values. 
  
  ## Input adjm1,2 should be matrices. 
  
  ##  TODO: Implement (and also calculate) what to do for non-deault model types. 
  if (model_type == "default") {
    fit1 = aggstat_single(fitModel(NetS, adjm1), adjm1)
    fit2 = aggstat_single(fitModel(NetS, adjm2), adjm2)
    fitc = aggstat_single(fitModel(NetS, adjm1*adjm2), adjm1*adjm2)
    
    px = fit1$x / fit1$n
    py = fit2$x / fit2$n
    pc = fitc$x / fitc$n
    
    nums = pc - px*py
    dens = sqrt(px * (1-px) * py * (1-py))
    cor_by_dyadgroup = nums/dens
    cor_by_dyadgroup[nums == 0] = 0 ## zero out any 0/0. 
    
    ## Returns a per-dyad-group adjustment (returns the estimated 'df' for that dyad group)
    return(1 - cor_by_dyadgroup) ## If correlation = 0, then the df = 1. 
  }
}

#' Computes a distance between two network models
#' 
#' @param NetM1 [\code{\link{NetworkModel}}] :: First network model
#' @param NetM2 [\code{\link{NetworkModel}}] :: Second network model
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
computeDist = function(NetM1, NetM2, type = "KLsym") {
  p1 = getEdgeProbMat(NetM1); p2 = getEdgeProbMat(NetM2)
  
  h_kldist = function(p1, p2) {
    diag(p1) <- 0.5; diag(p2) <- 0.5
    q1 = 1 - p1; q2 = 1 - p2
    ## TODO: Check wtf this distance is. 
    res = log(q1) - log(q2) + p1 * (log( (p1 / q1) / (p2 / q2 ) ) )
    diag(res) <- 0
    return(sum(res))
  }
  
  if (type == "KL") { 
    return(h_kldist(p1, p2))
  } else if (type == "KLsym") {
    return( (h_kldist(p1,p2)+h_kldist(p2,p1)) / 2)
  }
  stop("Invalid type of distance measure inputted")
  return(NULL)
}
