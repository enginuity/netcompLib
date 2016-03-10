## This is a collection of functions relating to the likelihood function for the extended models (density difference, correlation)


# Functions to aggregate over edge groups ---------------------------------


#' Compute sufficient statstics for a single network
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model with respect to which structure would be computed
#' @param adja [matrix/array-int] :: Adjacency matrix or arrays
#' 
#' @return [list] :: List with following entries
#' \itemize{
#'  \item n -- [vector-int] :: Edge group sizes
#'  \item x -- [vector-int] :: Dyad group edge presence counts
#'  \item names -- [vector-char] :: Dyad group IDs
#' }
#' 
#' @export
#' 
aggstat_single = function(NetM, adja) {
  if (length(dim(adja)) == 2) { 
    Nobs = 1
    adjm = adja
  } else {
    Nobs = dim(adja)[3]
    adjm = apply(adja, c(1,2), sum)
  }
  
  m = getEdgeProbMat(NetM, 'group')
  subset = lower.tri(m)
  inds = m[subset]
  
  xs = adjm[subset]
  
  return(lapply(list(
    n = Nobs * tapply(xs, inds, function(x) { sum(!is.na(x)) }), 
    x = tapply(xs, inds, sum, na.rm = TRUE), 
    names = names(tapply(xs, inds, sum))
  ), unname))
}


#' Compute sufficient statistics for density-difference estimation
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model with respect to which structure would be computed
#' @param adja1 [array-int] :: Adjacency matrix or arrays 1
#' @param adja2 [array-int] :: Adjacency matrix or arrays 2
#' 
#' @return [list] :: List with following entries
#' \itemize{
#'  \item n -- [vector-int] :: Edge group sizes
#'  \item x -- [vector-int] :: Dyad group edge presence counts in input 1
#'  \item y -- [vector-int] :: Dyad group edge presence counts in input 2
#'  \item names -- [vector-char] :: Dyad group IDs
#' }
#' 
#' @export
#' 
aggstat_dendiff = function(NetM, adja1, adja2) {
  m = getEdgeProbMat(NetM, 'group')
  subset = lower.tri(m)
  inds = m[subset]
  if (length(dim(adja1)) == 2) { xraw = adja1 } else { xraw = apply(adja1, c(1,2), sum) }
  if (length(dim(adja2)) == 2) { yraw = adja2 } else { yraw = apply(adja2, c(1,2), sum) }
  xs = xraw[subset]; ys = yraw[subset]
  
  return(lapply(list(
    n = tapply(xs, inds, function(x) { sum(!is.na(x)) }), 
    x = tapply(xs, inds, sum, na.rm = TRUE), 
    y = tapply(ys, inds, sum, na.rm = TRUE), 
    names = names(tapply(ys, inds, sum))
  ), unname))
}


#' Compute sufficient statistics for correlated-model estimation
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model with respect to which structure would be computed
#' @param adja1 [array-int] :: Adjacency matrix or arrays 1
#' @param adja2 [array-int] :: Adjacency matrix or arrays 2
#' 
#' @return [list] :: List with following entries
#' \itemize{
#'  \item n -- [vector-int] :: Edge group sizes
#'  \item names -- [vector-char] :: Dyad group IDs
#'  \item C - [matrix-int] :: columns are counts of 11, 10, 01, 00 cases respectively
#' }
#' 
#' @export
#' 
aggstat_corr = function(NetM, adja1, adja2) {
  m = getEdgeProbMat(NetM, 'group')
  subset = lower.tri(m)
  inds = m[subset]
  if (length(dim(adja1)) == 2) { xraw = adja1 } else { xraw = apply(adja1, c(1,2), sum) }
  if (length(dim(adja2)) == 2) { yraw = adja2 } else { yraw = apply(adja2, c(1,2), sum) }
  xs = xraw[subset]; ys = yraw[subset]
  xy = xs+ys
  xyd = xs-ys
  
  ## C matrix returned hsould be c11, c10, c01, c00
  return(c(lapply(list(
    n = tapply(xs, inds, function(x) { sum(!is.na(x)) }), 
    names = names(tapply(ys, inds, sum))), unname), 
    list(C = cbind(c11 = tapply(xy == 2, inds, sum, na.rm = TRUE), 
                   c10 = tapply(xyd == 1, inds, sum, na.rm = TRUE),
                   c01 = tapply(xyd == -1, inds, sum, na.rm = TRUE), 
                   c00 = tapply(xy == 0, inds, sum, na.rm = TRUE)))
  ))
}


#' Get dyad group probabilities
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model to get group probabilities of
#' 
#' @return [list] :: Two elements: 
#' \itemize{
#' \item probs -- [vector-double] :: Probabilities for each dyad group
#' \item names -- [vector-char] :: Names for each dyad group for ID purposes
#' }
#' 
#' @export
#' 
get_dyadgroup_prob = function(NetM) {
  ## TODO: Rename output 'names' by 'ids' in this and all edge group aggregation functions
  
  temp = aggstat_single(NetM, getEdgeProbMat(NetM))
  return(list(probs = temp$x / temp$n, names = temp$names))
}


#' Set dyad group probabilities
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model to reassign probabilities of
#' @param ids [vector-char OR vector-int] :: Dyad group IDs
#' @param probs [vector-double] :: Probabilities to assign to dyad groups
#' 
#' @return [\code{\link{NetworkModel}}] :: Adjusted network model with reassigned probabilities
#' 
#' @export
#' 
set_dyadgroup_prob = function(NetM, ids, probs) {
  ids = as.numeric(ids)
  ## TODO: Make this repsect standard order (based on ids?), since can be buggy if out of order! 
  
  if (inherits(NetM, "NetworkModelSBM")) {
    for(j in seq_along(ids)) {
      NN = NetM@Nnodes
      s = ids[j] %/% NN; r = ids[j] %% NN
      NetM@probmat[r,s] = probs[j]; NetM@probmat[s,r] = probs[j]
    }
    
  } else if (inherits(NetM, "NetworkModelHRG")) {
    NetM@prob = probs
    
  } else if (inherits(NetM, "NetworkModelRND")) {
    NetM@prob = probs
    
  } else {
    stop("Invalid NetM object")
  } 
  
  return(NetM)
}


# Likelihood Functions -- based on edge groups ----------------------------


#' Optimize a log-likelihood - helper
#' 
#' @param nparam [int] :: number of parameters needing randomization
#' @param optim_tries [int] :: number of runs of optimization
#' @param fn [function] :: log-likelihood function
#' @param gn [function] :: gradient function
#' @param ... [] :: other parameters to pass onto \code{\link{optim}}
#' 
#' @return [list] :: best output of optim
#' 
#' @export
#' 
h_optim = function(nparam, optim_tries, fn, gn, ...) {
  bestval = -Inf
  for(i in 1:optim_tries) { 
    temp = optim(rnorm(n = nparam), fn = fn, gr = gn, control = list(fnscale = -1), method = "BFGS", ...)
    if (temp$value < 0) {
      
      if (temp$value > bestval) {
        bestval = temp$value; best = temp
      }
    } else {
      cat("badrun.")
    }
  }
  return(best)
}


#' Converts parameters to probabilities in correlation model
#' 
#' @param rho [double] :: correlation parameter
#' @param a [vector-double] :: network1 parameters
#' @param b [vector-double] :: network2 parameters
#' 
#' @return [matrix-double] :: table with 4 columns -- probs of 11, 10, 01, 00
#' 
#' @export
#' 
hCorr_paramToProb = function(rho, a, b = a) {
  lambda = -log(exp(a+b+rho) + exp(a) + exp(b) + 1)
  df = cbind(exp(a+b+rho+lambda), exp(a+lambda), exp(b+lambda), exp(lambda))
  colnames(df) = c("11", "10", "01", "00")
  return(df)
}


#' Converts probabilities to parameters in correlation model
#' 
#' @param pt [matrix-double] :: table with 4 columns -- probs of 11, 10, 01, 00
#' 
#' @return [list] :: Named list of parameters in correlation model -- names are: 'a', 'b', 'lambda', 'rho'
#' 
#' @export
#' 
hCorr_probToParam = function(pt) {
  lambda = log(pt[,4])
  b = log(pt[,3]) - lambda
  a = log(pt[,2]) - lambda
  rho = log(pt[,1]) - lambda - a - b
  return(list(a = a, b = b, lambda = lambda, rho = rho))
}


#' Compute log-likelihood function for density-difference
#' 
#' @param t [vector-double] :: c(b,a) - b is single value; a is vector
#' @param x [vector-int] :: edge counts in net 1
#' @param y [vector-int] :: edge counts in net 2
#' @param n [vector-int] :: total edge group size
#' 
#' @return [double] :: value of log-likelihood
#' 
#' @export
#' 
llFx_dendiff = function(t,x,y,n) {
  b = t[1]; a = t[-1]
  return(sum(x*a - n * log(1 + exp(a)) + y * (a+b) - n * log(1 + exp(a+b))))
}


#' Compute gradient vector of log-likelihood in density-difference model
#' 
#' @param t [vector-double] :: c(b,a) - b is single value; a is vector
#' @param x [vector-int] :: edge counts in net 1
#' @param y [vector-int] :: edge counts in net 2
#' @param n [vector-int] :: total edge group size
#' 
#' @return [vector-double] :: gradient vector
#' 
#' @export
#' 
llGrFx_dendiff = function(t,x,y,n) {
  oe = function(x) { exp(x) / (1+exp(x)) }
  b = t[1]; a = t[-1]
  return(c(sum(y - n*oe(a+b)), x + y - n * (oe(a) + oe(a+b))))
}


#' Computes Log-likelihood for null-hypothesis correlation case
#' 
#' @param t [vector-double] :: c(rho, a_v)
#' @param C [matrix-int] :: Dyad group counts, in format: cbind(c11, c10, c01, c00)
#' @param n [vector-int] :: Dyad group sizes
#' 
#' @return [double] :: Log-likelihood at current point
#' 
#' @export
#' 
llFx_cnull = function(t, C, n) {
  N = nrow(C)
  
  rho = t[1]
  a = t[-1]
  
  lambda = - log( exp(a+a+rho) + exp(a) + exp(a) + 1)
  return(sum(C[,1]*(a+a+rho) + C[,2]*a + C[,3]*a + n*lambda))
}


#' Computes gradient of log-likelihood for null-hypothesis correlation case
#' 
#' @param t [vector-double] :: c(rho, a_v)
#' @param C [matrix-int] :: Dyad group counts, in format: cbind(c11, c10, c01, c00)
#' @param n [vector-int] :: Dyad group sizes
#' 
#' @return [vector-double] :: Gradient at current point
#' 
#' @export
#' 
llGrFx_cnull = function(t,C,n) {
  rho = t[1]; a = t[-1]
  # term1 = exp(2 * theta_k + rho); term2 = 2 * exp(theta_k)
  tr1 = exp(a + a + rho)
  tr2 = exp(a)
  den = tr1 + tr2 + tr2 + 1
  lambda = -log(den)
  
  grad = t * 0
  grad[1] = sum(C[,1] - n * tr1/den)
  grad[1+seq_along(a)] = 2 * C[,1] + C[,2] + C[,3] - 2 * n * ((tr1+tr2)/den)
  
  return(grad)
}


#' Computes Log-likelihood for alternate-hypothesis correlation case
#' 
#' @param t [vector-double] :: c(rho, a_v, b_v)
#' @param C [matrix-int] :: Dyad group counts, in format: cbind(c11, c10, c01, c00)
#' @param n [vector-int] :: Dyad group sizes
#' 
#' @return [double] :: Log-likelihood at current point
#' 
#' @export
#' 
llFx_calt = function(t,C,n) {
  N = nrow(C)
  rho = t[1]; rt = t[-1]
  a = rt[seq_len(N)]; rt = rt[-seq_len(N)]
  b = rt
  lambda = - log(exp(a + b + rho) + exp(a) + exp(b) + 1)
  
  return(sum(C[,1]*(a+b+rho) + C[,2]*a + C[,3]*b + n*lambda))
}


#' Computes gradient of log-likelihood for alternate-hypothesis correlation case
#' 
#' @param t [vector-double] :: c(rho, a_v, b_v)
#' @param C [matrix-int] :: Dyad group counts, in format: cbind(c11, c10, c01, c00)
#' @param n [vector-int] :: Dyad group sizes
#' 
#' @return [vector-double] :: Gradient at current point
#' 
#' @export
#' 
llGrFx_calt = function(t, C, n) {
  N = nrow(C)
  rho = t[1]; rt = t[-1]
  a = rt[seq_len(N)]; rt = rt[-seq_len(N)]
  b = rt
  
  ## term1 = exp(a + b + rho)
  tr1 = exp(a + b + rho); tr2 = exp(a); tr3 = exp(b); 
  den = tr1 + tr2 + tr3 + 1
  lambda = - log(den)
  
  grad = t * 0
  grad[1] = sum(C[,1] - n * tr1/den)
  grad[seq_len(N) + 1] = C[,1] + C[,2] - n * (tr1 + tr2) / den
  grad[seq_len(N) + 1 + N] = C[,1] + C[,3] - n * (tr1 + tr3) / den
  return(grad)
}



# Function to add computeLik outputs --------------------------------------

addComputeLik = function(cl1, cl2) {
  ## Adds the lists for outputs of computeLik. This is needed for alternate hypotheses... 
  cl1$sum = cl1$sum + cl2$sum
  cl1$bynode = cl1$bynode + cl2$bynode
  cl1$group_ll = cl1$group_ll + cl2$group_ll
  return(cl1)
}


