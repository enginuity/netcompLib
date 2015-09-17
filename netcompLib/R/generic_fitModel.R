##@S Function to create a NetworkModel object from a NetworkStruct object

setGeneric("fitModel", function(NetS, adja, mode, optim_tries) standardGeneric("fitModel"))


#' Estimates the best-fit model for given structure
#' 
#' Note: if adja is NULL, will fill in random probabilities drawn from an uniform distribution
#' 
#' @param NetS [\code{\link{NetworkStruct}}] :: Edge partition to fit model on
#' @param adja [matrix/array] :: Adjacency array or matrix
#' @param mode [char] :: Method to fit model -- choices are: "default", "densitydiff", "corr-global-null", "corr-global"
#' @param optim_tries [int] :: Number of attempts at optimization of likelihood function
#' 
#' @return [\code{\link{NetworkModel}}] :: Returns the best-fit model for each structure. If input is a \code{\link{NetworkStructList}}, output will actually be a list of \code{\link{NetworkModel}}s. 
#' 
#' @export
#' 
fitModel = function(NetS, adja, mode = "default", optim_tries = 10) {
  stop("No implementation in template")
  return(NULL)
}


fitModel.NetworkStruct = function(NetS, adja, mode = "default", optim_tries = 10) { 
  stop("No implmentation for template NetworkStruct")
}


fitModel.NetworkStructList = function(NetS, adja, mode = "default", optim_tries = 10) {
  return(lapply(NetS@models, function(x) { fitModel(x, adja) } ))
}


fitModel.NetworkStructSBM = function(NetS, adja, mode = "default", optim_tries = 10) {
  ## mode -- default
  ## if mode = densitydiff, then, assume global density difference parameter
  ## TODO: Re-implement this?!?! the modes may not be the best way to deal with this... 
  
  res = NetworkModel(set_model_param(Nnodes = getNnodes(NetS), type = "block", block_assign = NetS@groups))
  Nobs = dim(adja)[3]
  
  NG = max(res@assign)
  gr = res@assign
  res@probmat = matrix(NA, NG, NG)
  if (mode == "default") {
    if (Nobs == 1) { adjm = adja[,,1] }
    if (Nobs > 1) { adjm = apply(adja, c(1,2), sum) }
    for(i in 1:NG) { for(j in 1:NG) {
      ii = which(gr == i); ij = which(gr == j)
      ci = sum(gr == i); cj = sum(gr == j)
      if (i != j) { res@probmat[i,j] = sum(adjm[ii,ij]) / (Nobs * ci * cj) }
      if (i == j) { res@probmat[i,j] = sum(adjm[ii,ij]) / (Nobs * ci * (ci - 1)) }
    }}
    return(res) 
  }
  
  if (mode == "densitydiff") {
    ## TODO: Allow fitting when inputting longer adjacency arrays? (or perhaps a pair of them)
    
    temp = aggstat_dendiff(res, adja[,,1], adja[,,2]) ## should work for all struct types
    x = temp$x; y = temp$y; n = temp$n
    
    bestval = -Inf
    for(i in 1:optim_tries) { 
      temp = optim(rnorm(n = length(n)+1), fn = llFx_dendiff, gr = llGrFx_dendiff, x=x, y=y, n=n, control = list(fnscale = -1), method = "BFGS")
      if (temp$value > bestval) {
        bestval = temp$value; best = temp
      }
    }
    
    m1 = 3
    NetworkModelPair(m1, m2, is_null = FALSE)
    
    return(list(res, temp)) ## TODO: Fix output
    ## should output networkmodel pair! 
  }
  
  
  
  
  
  ## TO REDO. 
  if (mode == "corr-global-null") {
    sample_start = function(n) { 
      t = rnorm(n)
      p = 1/(1 + exp(-t[-1]))
      rmin = max(-c(p/(1-p), (1-p)/p))
      r = runif(1, min = rmin, max = 1) * .8
      t[1] = log((1+r)/(1-r))
      return(t)
    }
    
    a1 = adja[,,1]; a2 = adja[,,2]; diag(a1) = -1; diag(a2) = -1
    llfx = function(t,C,n) {
      # C = cbind(c11, c10, c01, c00); n = vector of total counts; t = c(rho, p)
      # p, r in logit form
      p = 1/(1 + exp(-t[-1]))
      r = 2 / (1 + exp(-t[1])) - 1 
      q = 1- p; pq = p * q 
      
      p11 = p^2 + r * pq
      p00 = q^2 + r * pq
      p10 = (1-r) * pq
      if (!all(is_prob(c(p11, p10, p00)))) return(10^-239)
      
      c11 = C[,1]; c10 = C[,2]; c01 = C[,3]; c00 = C[,4]
      res = c11 * log(p11) + c00 * log(p00) + (c10 + c01) * log(p10)
      return(sum(res))
    }
    grfx = function(t, C, n) {
      p = 1/(1 + exp(-t[-1]))
      r = 2 / (1 + exp(-t[1])) - 1 
      q = 1 - p; pq = p * q
      
      p11 = p^2 + r * pq
      p00 = q^2 + r * pq
      p10 = (1-r) * pq
      
      c11 = C[,1]; c10 = C[,2]; c01 = C[,3]; c00 = C[,4]
      res = t * 0
      
      res[1] = sum(c11*pq/p11 + c00*pq/p00 - (c10+c01)*pq/p10)
      res[-1] = c11*(2*p+r*(1-2*p))/p11 + c00*(-2*(1-p)+r*(1-2*p))/p00 + (c10*c01)*(1-2*p)/(p*q)
      return(res/1000)
    }
    c11m = matrix(0, NG, NG); c00m = matrix(0, NG, NG)
    c10m = matrix(0, NG, NG); c01m = matrix(0, NG, NG)
    tm = matrix(0, NG, NG)
    for(i in 1:NG) { for(j in 1:NG) {
      ii = which(gr == i); ij = which(gr == j)
      ci = sum(gr == i); cj = sum(gr == j)
      if (i == j) { 
        c11m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 1, na.rm = TRUE)/2
        c10m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 0, na.rm = TRUE)/2
        c01m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 1, na.rm = TRUE)/2
        c00m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 0, na.rm = TRUE)/2
        tm[i,j] = ci * (ci -1)/2
      }
      if (i != j) { 
        c11m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 1, na.rm = TRUE)
        c10m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 0, na.rm = TRUE)
        c01m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 1, na.rm = TRUE)
        c00m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 0, na.rm = TRUE)
        tm[i,j] = ci*cj 
      }
    }}
    
    C = cbind(c11m[lower.tri(c11m, diag = TRUE)], c10m[lower.tri(c11m, diag = TRUE)],
              c01m[lower.tri(c11m, diag = TRUE)], c00m[lower.tri(c11m, diag = TRUE)])
    n = tm[lower.tri(tm, diag = TRUE)]
    
    bestval = -Inf
    v = rep(1:optim_tries)
    for(i in 1:optim_tries) { 
      print(i)
      temp = optim(sample_start(length(n) + 1), fn = llfx, gr = grfx, C = C, n=n, control = list(fnscale = -1), method = "BFGS")
      if (temp$value > bestval) {
        bestval = temp$value; best = temp
      }
      v[i] = temp$value
    }
    
    return(list(res, best, v))
  }
  
  
  if (mode == "corr-global") {
    sample_start = function(n) { 
      t = rnorm(n)
      p = 1/(1 + exp(-t[-1]))
      r = 0
      t[1] = log((1+r)/(1-r))
      return(t)
    }
    
    a1 = adja[,,1]; a2 = adja[,,2]; diag(a1) = -1; diag(a2) = -1
    llfx = function(t,D,n) {
      # D = cbind(c11, c10, c01, c00); n = vector of total counts; t = c(rho, p)
      # p1,p2,r in logit form
      N = nrow(D)
      r = 2 / (1 + exp(-t[1])) - 1; rt = t[-1]
      p1 = 1/(1 + exp(-t[seq_len(N)])); rt = rt[-seq_len(N)]
      p2 = 1/(1 + exp(-rt))
      c11 = D[,1]; c10 = D[,2]; c01 = D[,3]; c00 = D[,4]
      A = r * sqrt(p1 * p2 * (1 - p1) * (1 - p2))
      
      p11 = p1 * p2 + A
      p10 = p1 * (1-p2) - A
      p01 = (1-p1) * p2 - A
      p00 = (1-p1)*(1-p2) + A
      if (!all(is_prob(c(p11, p10, p01, p00)))) return(10^-239)
      res = c11 * log(p11) + c00 * log(p00) + c10 * log(p01) + c01 * log(p00)
      return(sum(res))
    }
    grfx = function(t, D, n) {
      # D = cbind(c11, c10, c01, c00); n = vector of total counts; t = c(rho, p)
      # p1,p2,r in logit form
      N = nrow(D)
      r = 2 / (1 + exp(-t[1])) - 1; rt = t[-1]
      p1 = 1/(1 + exp(-t[seq_len(N)])); rt = rt[-seq_len(N)]
      p2 = 1/(1 + exp(-rt))
      q1 = 1 - p1; q2 = 1 - p2
      c11 = D[,1]; c10 = D[,2]; c01 = D[,3]; c00 = D[,4]
      s = sqrt(p1 * p2 * q1 * q2)
      A = r * s
      dp1 = (1-2*p1)*p2*q2/(2*s)
      dp2 = (1-2*p2)*p1*q1/(2*s)
      
      p11 = p1 * p2 + A
      p10 = p1 * q2 - A
      p01 = q1 * p2 - A
      p00 = q1 * q2 + A
      
      res = t * 0
      res[1] = sum(c11*s/p11 + c00*s/p00 - c10*s/p10 - c01*s/p01)
      res[1 + seq_len(N)] = c11*(p2+dp1)/p11 + c10*(q2-dp1)/p10 + c01*(-p2-dp1)/p01 + c00*(-q2+dp1)/p00
      res[1 + N + seq_len(N)] = c11*(p1+dp2)/p11 + c10*(q1-dp2)/p10 + c01*(-p1-dp2)/p01 + c00*(-q1+dp2)/p00
      return(res/1000)
    }
    
    c11m = matrix(0, NG, NG); c00m = matrix(0, NG, NG)
    c10m = matrix(0, NG, NG); c01m = matrix(0, NG, NG)
    tm = matrix(0, NG, NG)
    for(i in 1:NG) { for(j in 1:NG) {
      ii = which(gr == i); ij = which(gr == j)
      ci = sum(gr == i); cj = sum(gr == j)
      if (i == j) { 
        c11m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 1, na.rm = TRUE)/2
        c10m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 0, na.rm = TRUE)/2
        c01m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 1, na.rm = TRUE)/2
        c00m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 0, na.rm = TRUE)/2
        tm[i,j] = ci * (ci -1)/2
      }
      if (i != j) { 
        c11m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 1, na.rm = TRUE)
        c10m[i,j] = sum(a1[ii,ij] == 1 & a2[ii,ij] == 0, na.rm = TRUE)
        c01m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 1, na.rm = TRUE)
        c00m[i,j] = sum(a1[ii,ij] == 0 & a2[ii,ij] == 0, na.rm = TRUE)
        tm[i,j] = ci*cj 
      }
    }}
    
    D = cbind(c11m[lower.tri(c11m, diag = TRUE)], c10m[lower.tri(c11m, diag = TRUE)],
              c01m[lower.tri(c11m, diag = TRUE)], c00m[lower.tri(c11m, diag = TRUE)])
    n = tm[lower.tri(tm, diag = TRUE)]
    
    bestval = -Inf
    v = rep(1:optim_tries)
    for(i in 1:optim_tries) { 
      temp = optim(sample_start(2*length(n) + 1), fn = llfx, D = D, n=n, control = list(fnscale = -1))
      if (temp$value > bestval) {
        bestval = temp$value; best = temp
      }
      v[i] = temp$value
    }
    
    return(list(res, best, v))
  }
  
  stop("Invalid 'mode'")
}


fitModel.NetworkStructRND = function(NetS, adja, mode = "default", optim_tries = 10) {
  if (is.null(adja)) {
    stop("Not implemented")
    ## TODO: Fill in code for conversion
  } else {
    ## TODO: Fill in code for estimating probabilities. 
    stop("Not implemented for non-null adjacency array")
  }
}


fitModel.NetworkStructHRG = function(NetS, adja, mode = "default", optim_tries = 10) {
  if (is.null(adja)) {
    stop("Not implemented")
    ## TODO: Fill in code for conversion
  } else {
    ## TODO: Fill in code for estimating probabilities. 
    stop("Not implemented for non-null adjacency array")
  }
}


# setMethod ---------------------------------------------------------------
setMethod("fitModel", signature(NetS = "NetworkStruct"), fitModel.NetworkStruct)
setMethod("fitModel", signature(NetS = "NetworkStructList"), fitModel.NetworkStructList)
setMethod("fitModel", signature(NetS = "NetworkStructSBM"), fitModel.NetworkStructSBM)
setMethod("fitModel", signature(NetS = "NetworkStructRND"), fitModel.NetworkStructRND)
setMethod("fitModel", signature(NetS = "NetworkStructHRG"), fitModel.NetworkStructHRG)
