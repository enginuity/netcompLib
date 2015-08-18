##@S Function to create a NetworkModel object from a NetworkStruct object

setGeneric("fitModel", function(NetS, adja, mode) standardGeneric("fitModel"))


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (fitModel)
#' Estimates the best-fit model for given structure
#' 
#' Note: if adja is NULL, will fill in random probabilities drawn from an uniform distribution
#' 
#' @param NetS NetworkStruct object
#' @param adja adjacency array/ matrix
#' @param mode temp
#' 
#' @return NetworkModel object containing best-fit model per structure. 
#' 
#' @export
#' 
fitModel = function(NetS, adja, mode = "default") {
  stop("No implementation in template")
  return(NULL)
}


fitModel.NetworkStruct = function(NetS, adja, mode = "default") { 
  stop("No implmentation for template NetworkStruct")
}


fitModel.NetworkStructList = function(NetS, adja, mode = "default") {
  return(lapply(NetS@models, function(x) { fitModel(x, adja) } ))
}


fitModel.NetworkStructSBM = function(NetS, adja, mode = "default") {
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
    llfx = function(t,x,y,n) {
      ## x,y,n,a are vectorized; b should be constant; t = c(b,a)
      b = t[1]; a = t[-1]
      return(sum(x*a - n * log(1 + exp(a)) + y * (a+b) - n * log(1 + exp(a+b))))
    }
    
    xm = matrix(0, NG, NG); ym = xm; tm = matrix(0, NG, NG)
    for(i in 1:NG) { for(j in 1:NG) {
      ii = which(gr == i); ij = which(gr == j)
      ci = sum(gr == i); cj = sum(gr == j)
      if (i == j) { xm[i,j] = sum(adja[ii,ij,1])/2 ; ym[i,j] = sum(adja[ii,ij,2])/2; tm[i,j] = ci * (ci -1)/2}
      if (i != j) { xm[i,j] = sum(adja[ii,ij,1]) ; ym[i,j] = sum(adja[ii,ij,2]); tm[i,j] = ci*cj }
    }}
    
    x = xm[lower.tri(xm, diag = TRUE)]; y = ym[lower.tri(ym, diag = TRUE)]
    n = tm[lower.tri(tm, diag = TRUE)]
    temp = optim(rnorm(n = length(n)+1), fn = llfx, x=x, y=y, n=n, control = list(fnscale = -1))
    return(list(res, temp))
  }
  
  if (mode == "corr-global-null") {
    a1 = adja[,,1]; a2 = adja[,,2]; diag(a1) = -1; diag(a2) = -1
    llfx = function(t,C,n) {
      # C = cbind(c11, c10, c01, c00); n = vector of total counts; t = c(rho, p)
      # p, r in logit form
      p = 1/(1 + exp(-t[-1]))
      r = 2 / (1 + exp(-t[1])) - 1 
      q = 1- p; pq = p * q 
      c11 = C[,1]; c10 = C[,2]; c01 = C[,3]; c00 = C[,4]
      res = c11 * log(p^2 + r * pq) + c00 * log(q^2 + r * pq) + (c10 + c01) * log((1-r) * pq)
      if (any(is.nan(res))) return(-100000)
      return(sum(res))
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
    for(i in 1:10) { 
      temp = optim(rnorm(length(n) + 1), fn = llfx, C = C, n=n, control = list(fnscale = -1))
      if (temp$value > bestval) {
        bestval = temp$value; best = temp
      }
    }
    
    return(list(res, temp))
  }
  
  
  if (mode == "corr-global") {
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
      res = c11 * log(p1 * p2 + A) + c00 * log((1-p1)*(1-p2) + A) + c10 * log(p1 * (1-p2) - A) + c01 * log((1-p1) * p2 - A)
      if (any(is.nan(res))) return(-100000)
      return(sum(res))
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
    for(i in 1:10) { 
      temp = optim(rnorm(2*length(n) + 1), fn = llfx, D = D, n=n, control = list(fnscale = -1))
      if (temp$value > bestval) {
        bestval = temp$value; best = temp
      }
    }
    return(list(res, best))
  }
  
  stop("Invalid 'mode'")
}


fitModel.NetworkStructRND = function(NetS, adja, mode = "default") {
  if (is.null(adja)) {
    stop("Not implemented")
    ## TODO: Fill in code for conversion
  } else {
    ## TODO: Fill in code for estimating probabilities. 
    stop("Not implemented for non-null adjacency array")
  }
}


fitModel.NetworkStructHRG = function(NetS, adja, mode = "default") {
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
