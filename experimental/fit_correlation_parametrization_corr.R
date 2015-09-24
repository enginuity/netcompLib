## Code used for fitting models for global correlated -- this is the true correlation parametrization. 

# TODO: Save this code? This is the alternate parametrization... somehow allow this parametrization to be used? 

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
    p1 = 1/(1 + exp(-rt[seq_len(N)])); rt = rt[-seq_len(N)]
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
