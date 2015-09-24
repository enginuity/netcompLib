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

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (aggstat_corr)
#' <What does this function do>
#' 
#' @param NetM temp
#' @param adja1 temp
#' @param adja2 temp
#' 
#' @return temp
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
  
  ## C matrix returned hsould be c11, c10, c01, c00
  return(c(lapply(list(n = tapply(inds, inds, function(x) { sum(x > 0) }), names = names(tapply(ys, inds, sum))), unname), 
         list(C = cbind(c11 = tapply(xy == 2, inds, sum), c10 = tapply(xy == xs, inds, sum),
                        c01 = tapply(xy == ys, inds, sum), c00 = tapply(xy == 0, inds, sum)))))
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (reassign_edgegroup_prob)
#' <What does this function do>
#' 
#' @param NetM temp
#' @param ids temp
#' @param probs temp
#' 
#' @return temp
#' 
#' @export
#' 
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


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (llFx_cnull)
#' <What does this function do>
#' 
#' @param t temp
#' @param C temp
#' @param n temp
#' 
#' @return temp
#' 
#' @export
#' 
llFx_cnull = function(t, C, n) {
  # C = cbind(c11, c10, c01, c00); n = vector of total counts; t = c(rho, theta_v)
  rho = t[1]
  theta = t[-1]
  
  lambda = -log(exp(2 * theta + rho) + 2 * exp(theta) + 1)
  return(sum(C[,1] * (2 * theta + rho) + (C[,2] + C[,3]) * theta + n * lambda))
}

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (llGrFx_cnull)
#' <What does this function do>
#' 
#' @param t temp
#' @param C temp
#' @param n temp
#' 
#' @return temp
#' 
#' @export
#' 
llGrFx_cnull = function(t,C,n) {
  rho = t[1]; theta = t[-1]
  # term1 = exp(2 * theta_k + rho); term2 = 2 * exp(theta_k)
  tr1 = exp(2 * theta + rho)
  tr2 = 2 * exp(theta)
  den = tr1 + tr2 + 1
  lambda = -log(den)
  
  grad = t * 0
  grad[1] = sum(C[,1] - n * tr1/den)
  grad[-1] = 2 * C[,1] + C[,2] + C[,3] - 2 * n * (1 - 1/den)
  
  return(grad)
}

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (llFx_calt)
#' <What does this function do>
#' 
#' @param t temp
#' @param C temp
#' @param n temp
#' 
#' @return temp
#' 
#' @export
#' 
llFx_calt = function(t,C,n) {
  # t  = c(rho, a_v, b_v)
  N = nrow(C)
  rho = t[-1]; rt = t[-1]
  a = rt[seq_len(N)]; rt = rt[-seq_len(N)]
  b = rt
  lambda = - log(exp(a + b + rho) + exp(a) + exp(b) + 1)
  
  return(sum(C[,1]*(a+b+rho) + C[,2]*a + C[,3]*b + n*lambda))
}

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (llGrFx_calt)
#' <What does this function do>
#' 
#' @param t temp
#' @param C temp
#' @param n temp
#' 
#' @return temp
#' 
#' @export
#' 
llGrFx_calt = function(t, C, n) {
  N = nrow(C)
  rho = t[-1]; rt = t[-1]
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





