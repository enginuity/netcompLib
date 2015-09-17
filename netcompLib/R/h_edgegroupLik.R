## This is a collection of functions relating to the likelihood function for the extended models (density difference, correlation)


# Functions to aggregate over edge groups ---------------------------------

#' Aggregate statistics for density-difference estimation
#' 
#' @param NetM [] :: model
#' @param adja1 [] :: adjm 1
#' @param adja2 [] :: adjm 2
#' 
#' @return [] :: list of statistics: n, x, y -- n = total edges in edge group, x,y = edge appearance count in each edge group
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
  
  return(lapply(list(n = tapply(inds, inds, function(x) { sum(x > 0) }), x = tapply(xs, inds, sum), y = tapply(ys, inds, sum), names = names(tapply(ys, inds, sum))), unname))
}

reassign_edgegroup_prob = function(NetM, ids, probs) {
  if (inherits(NetM, "NetworkModelSBM")) {
    ids = as.numeric(ids)
    for(j in seq_along(ids)) {
      NN = NetM@Nnodes
      s = ids[j] %/% NN; r = ids[j] %% NN
      NetM@probmat[r,s] = probs[j]; NetM@probmat[s,r] = probs[j]
    }
    return(NetM)
  } else {
    ## Implement cases for HRG, RND
    stop("Case not implemented")
  } 
}
# Likelihood Functions -- based on edge groups ----------------------------


#' Compute log-likelihood function for density-difference
#' 
#' @param t [] :: c(b,a) - b is single value; a is vector
#' @param x [] :: edge counts in net 1
#' @param y [] :: edge counts in net 2
#' @param n [] :: total edge group size
#' 
#' @return [] :: value of log-likelihood
#' 
#' @export
#' 
llFx_dendiff = function(t,x,y,n) {
  b = t[1]; a = t[-1]
  return(sum(x*a - n * log(1 + exp(a)) + y * (a+b) - n * log(1 + exp(a+b))))
}


#' Compute gradient vector of log-likelihood in density-difference model
#' 
#' @param t [] :: c(b,a) - b is single value; a is vector
#' @param x [] :: edge counts in net 1
#' @param y [] :: edge counts in net 2
#' @param n [] :: total edge group size
#' 
#' @return [] :: gradient vector
#' 
#' @export
#' 
llGrFx_dendiff = function(t,x,y,n) {
  oe = function(x) { exp(x) / (1+exp(x)) }
  b = t[1]; a = t[-1]
  return(c(sum(y - n*oe(a+b)), x + y - n * (oe(a) + oe(a+b))))
}



