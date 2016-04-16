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
  ## mode -- default
  ## if mode = densitydiff, then, assume global density difference parameter
  ## TODO: [Issue #33] Re-implement this?!?! the modes may not be the best way to deal with this... 
  
  res = extractModel(NetS)
  
  ## If adja is a list -- call fitModel separately on each entry
  ## TODO: [Update documentation] document this as an allowable case
  if (is.list(adja)) {
    return(lapply(adja, function(x) { fitModel(NetS, x, mode, optim_tries) }))
  }
  
  
  ## If adja is NULL, just fill in random probabilities, which extractModel does. 
  if (is.null(adja)) {
    # if (mode == "default") { 
    #   return()
    # }
    return(res)
    ## TODO: [Issue #33] This doesn't work when mode is NOT default. Need to handle this properly eventually. 
  }
  
  ## When the adja is not NULL:
  Nobs = dim(adja)[3]
  if (length(dim(adja)) == 2) { Nobs = 1 }
  
  if (mode == "default") {
    aggstat = aggstat_single(res, adja)
    
    return(set_dyadgroup_prob(res, aggstat$names, probs = aggstat$x/aggstat$n))
  } else if (mode == "densitydiff") {
    ## TODO: Allow fitting when inputting longer adjacency arrays? (or perhaps a pair of them)
    
    aggstat = aggstat_dendiff(res, adja[,,1], adja[,,2]) ## should work for all struct types
    x = aggstat$x; y = aggstat$y; n = aggstat$n
    
    best = h_optim(nparam = length(n) + 1, optim_tries = optim_tries, fn = llFx_dendiff, gn = llGrFx_dendiff, x = x, y = y, n = n)
    
    m1 = set_dyadgroup_prob(res, aggstat$names, probs = faraway::ilogit(best$par[-1]))
    m2 = set_dyadgroup_prob(res, aggstat$names, probs = faraway::ilogit(best$par[-1] + best$par[1]))
    return(NetworkModelPair(m1 = m1, m2 = m2, is_null = FALSE, model_type = "densitydiff", addl_param = list(dd_param_add = best$par[1])))
  }
  
  if (mode == "corr-global-null") {
    aggstat = aggstat_corr(res, adja[,,1], adja[,,2])
    n = aggstat$n; C = aggstat$C
    
    best = h_optim(nparam = length(n)+1, optim_tries = optim_tries, fn = llFx_cnull, gn = llGrFx_cnull, n = n, C = C)
    rho = best$par[1]; a = best$par[-1]
    pt = hCorr_paramToProb(rho, a)
    
    m1 = set_dyadgroup_prob(res, aggstat$names, probs = pt[,1] + pt[,2])
    m2 = set_dyadgroup_prob(res, aggstat$names, probs = pt[,1] + pt[,3])
    return(NetworkModelPair(m1 = m1, m2 = m2, is_null = FALSE, model_type = "correlated", addl_param = list(c_param_corr = rho, c_param_a = a, c_param_b = a, c_names = aggstat$names)))
  }
  
  if (mode == "corr-global") {
    aggstat = aggstat_corr(res, adja[,,1], adja[,,2])
    n = aggstat$n; C = aggstat$C
    
    best = h_optim(nparam = 2*length(n)+1, optim_tries = optim_tries, fn = llFx_calt, gn = llGrFx_calt, n = n, C = C)
    rho = best$par[1]; a = best$par[1 + seq_along(n)]; b = best$par[1+length(n)+seq_along(n)]
    pt = hCorr_paramToProb(rho, a, b)
    
    m1 = set_dyadgroup_prob(res, aggstat$names, probs = pt[,1] + pt[,2])
    m2 = set_dyadgroup_prob(res, aggstat$names, probs = pt[,1] + pt[,3])
    return(NetworkModelPair(m1 = m1, m2 = m2, is_null = FALSE, model_type = "correlated", addl_param = list(c_param_corr = rho, c_param_a = a, c_param_b = b, c_names = aggstat$names)))
  }
  
  stop("Invalid 'mode'")
}


fitModel.NetworkStructList = function(NetS, adja, mode = "default", optim_tries = 10) {
  return(lapply(NetS@models, function(x) { fitModel(x, adja) } ))
}


# setMethod ---------------------------------------------------------------
setMethod("fitModel", signature(NetS = "NetworkStruct"), fitModel.NetworkStruct)
setMethod("fitModel", signature(NetS = "NetworkStructList"), fitModel.NetworkStructList)

