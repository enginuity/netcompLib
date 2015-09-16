## This is a collection of functions relating to the likelihood function for the extended models (density difference, correlation)


# Functions to aggregate over edge groups ---------------------------------

aggstat_dendiff = function(NetM, adja1, adja2) {
  m = getEdgeProbMat(NetM, 'group')
  subset = lower.tri(m)
  inds = m[subset]
  if (length(dim(adja1)) == 2) { xraw = adja1 } else { xraw = apply(adja1, c(1,2), sum) }
  if (length(dim(adja2)) == 2) { yraw = adja2 } else { yraw = apply(adja2, c(1,2), sum) }
  xs = xraw[subset]; ys = yraw[subset]
  
  return(lapply(list(n = tapply(inds, inds, function(x) { sum(x > 0) }), x = tapply(xs, inds, sum), y = tapply(ys, inds, sum)), unname))
}


# Likelihood Functions -- based on edge groups ----------------------------


llFx_dendiff = function(t,x,y,n) {
  ## x,y,n,a are vectorized; b should be constant; t = c(b,a)
  b = t[1]; a = t[-1]
  return(sum(x*a - n * log(1 + exp(a)) + y * (a+b) - n * log(1 + exp(a+b))))
}


llGrFx_dendiff = function(t,x,y,n) {
  ## x,y,n,a are vectorized; b should be constant; t = c(b,a)
  oe = function(x) { exp(x) / (1+exp(x)) }
  b = t[1]; a = t[-1]
  return(c(sum(y - n*oe(a+b)), x + y - n * (oe(a) + oe(a+b))))
}






