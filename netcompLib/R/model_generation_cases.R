## Contains code to generate and simulate from each of the network models. 

## TODO: [Obselete] 2015.03.12 -- all of this code has been copied into NetworkModel***.R's. 

## TODO: [Cleanup] Some variables are placeholders (that stood for previous parameters). They can safely be removed if it doesn't help the code... 

# SBM ---------------------------------------------------------------------

#' Generate block model parameters
#' 
#' @param Nnodes Number of nodes
#' @param model_param List of model parameters
#' 
#' @return List with block model information
#' 
#' @export
#' 
sample_model_block = function(Nnodes, model_param = set_model_param()) {
  K = model_param$block_nclass
  
  ## Function for adjusting block model probabilities  
  adjust_blockprobs = function(mod, avgden = 0.4, plimit = c(0.05, 0.95)) {
    classct = table(mod$assign)
    params = length(classct) * (length(classct) - 1 ) / 2 + length(classct)
    coordmat = matrix(NA, nrow = params, ncol = 4)
    colnames(coordmat) = c("row", "col", "ID", "count")
    
    ct = 1
    for(j in 1:length(classct)) { for(k in j:length(classct)) {
      coordmat[ct,] = c(j,k,ct, ifelse(j == k, classct[j] * (classct[j]-1)/2, classct[j] * classct[k]))
      ct = ct + 1
    }}
    
    samp_probvec = function(counts, avgden, plimit) {
      notvalid = TRUE
      while(notvalid) {
        test = runif(n = length(counts), min = plimit[1], max = plimit[2])
        test[1] = (sum(counts) * avgden - sum(counts[-1] * test[-1]))/counts[1]
        if (test[1] > plimit[1] & test[1] < plimit[2]) { notvalid = FALSE }
      }
      return(test)
    }
    probvec = samp_probvec(counts = coordmat[,4], avgden = avgden, plimit = plimit)
    
    pmat = matrix(NA, ncol = length(classct), nrow = length(classct))
    ct = 1
    for(j in 1:length(classct)) { for(k in j:length(classct)) {
      pmat[j,k] = probvec[ct]
      if (j != k) { pmat[k,j] = probvec[ct] }
      ct = ct + 1
    }}
    return(pmat)
  }
  
  group_assign = sample(1:K, size = Nnodes, replace = TRUE)
  prob_matrix = matrix(0, nrow = K, ncol = K)
  for(j in 1:K) { for(i in 1:K) {
    if (i <= j) {
      prob_matrix[i,j] = runif(1, model_param$pmin, model_param$pmax)
      if (i != j) {prob_matrix[j,i] = prob_matrix[i,j]}
    }
  }}
  
  res = list(assign = group_assign, probmat = prob_matrix)
  if (!is.null(model_param$block_avgdensity)) { 
    ## Set appropriate probability matrix
    newpmat = adjust_blockprobs(mod = res, avgden = model_param$block_avgdensity, plimit = c(model_param$pmin, model_param$pmax))
    res$probmat = newpmat
  }
  
  return(res)
}


#' Generate data from a block model
#' 
#' @param bm Block model information
#' @param Nobs Number of observations
#' 
#' @return Array of network observations
#' 
#' @export
#' 
sample_network_block = function(bm, Nobs = 1) {
  Nnodes= length(bm$assign)
  res = array(0, dim = c(Nnodes, Nnodes, Nobs))
  for(d in 1:Nobs) {
    for(k in 1:Nnodes) {
      for(j in 1:Nnodes) {
        if (k < j) {
          v = rbinom(n = 1, size = 1, prob = bm$probmat[bm$assign[k], bm$assign[j]])
          res[j,k,d] = v
          res[k,j,d] = v
        }
      }}
  }
  return(res)
}




# HRG ---------------------------------------------------------------------


## Code to generate HRG models

#' Generate a tree: fixed structure, node probabilities at random (uniform)
#' 
#' @param Nnodes Number of nodes in network
#' 
#' @return Tree object (list)
#' 
#' @export
#' 
starter_tree = function(Nnodes = 10) {
  # Written by Andrew Thomas
  
  #format: id, left child, right child, value.
  id1 <- Nnodes+1:(Nnodes-1)  #internal node ids.
  lcrc <- array(NA, c(2,Nnodes-1)); lcrc[1:(Nnodes-2)] <- 2:(Nnodes-1)+Nnodes; lcrc[(Nnodes-1):(2*(Nnodes-1))] <- 1:Nnodes
  children <- list(); for (kk in 1:dim(lcrc)[2]) children[[kk]] <- lcrc[,kk]
  
  val <- runif(Nnodes-1)
  
  parents <- rep(0, 2*Nnodes - 1)
  for (kk in 1:length(children)) parents[children[[kk]]] <- kk+Nnodes
  
  output <- list(prob=val,
                 children=children,
                 parents=parents,
                 nodes=Nnodes)
  
  return(output)
}

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sample_model_tree)
#' Generate a random tree structure / model
#' 
#' @param Nnodes Number of nodes in the network
#' @param model_param temp
#' 
#' @return HRG tree object
#' 
#' @export
#' 
sample_model_tree = function(Nnodes, model_param = set_model_param()) {
  
  ## TODO: [TEMP] Using these two variables shouldn't be necessary. 
  random_type = model_param$tree_type
  random_plimit = c(model_param$pmin, model_param$pmax)
  
  ## Helper functions
  node_to_add = function(inp_clist, cur_node, NN) {
    lr = 1
    if (runif(1) > 0.5) { lr = 2 }
    
    if (inp_clist[[cur_node - NN]][lr] == -1) {
      return(c(cur_node, lr))
    } else {
      return(node_to_add(inp_clist, inp_clist[[cur_node - NN]][lr], NN))
    }
  }
  
  find_first_empty = function(inp_clist) {
    nores = TRUE
    j = 0
    while(nores) {
      j = j + 1
      z = which(inp_clist[[j]] == -1)
      nores = (length(z) == 0)
    }
    return(c(j,z[1]))
  }
  
  
  ## Note that random_plimit only applies in the case of "random" random_type. (and for "left" random_type also). 
  res_tree = list()
  res_tree$nodes = Nnodes
  res_tree$prob = runif(Nnodes-1, min = random_plimit[1], max = random_plimit[2])
  
  if (random_type == "original") {
    return(starter_tree(Nnodes))
    
  } else if (random_type == "random") {
    clist = list()
    for(j in 1:(Nnodes-1)) {
      clist[[j]] = c(-1,-1)
    }
    
    for(j in (Nnodes+2):(2*Nnodes-1)) {
      z = node_to_add(clist, cur_node = (Nnodes+1), Nnodes)
      clist[[z[1] - Nnodes]][z[2]] = j
    }
    
    ## fill remaining children :
    for(j in sample(1:Nnodes, size = Nnodes)) {
      z = find_first_empty(clist)
      clist[[z[1]]][z[2]] = j
    }
    
    pars = rep(0, times = (2*Nnodes - 1))
    for(j in (Nnodes+1):(2*Nnodes-1)) {
      pars[clist[[j-Nnodes]]] = j
    }
    
  } else if (random_type == "left") {
    pars = rep(0, times = 2*Nnodes - 1)
    pars[(Nnodes+2):(2*Nnodes-1)] = (Nnodes+1):(2*Nnodes-2)
    pars[1:Nnodes] = sample(c( (Nnodes+1):(2*Nnodes-1), (2*Nnodes-1)), size = Nnodes)
    
  } else {
    stop(paste("ERROR: no such type", random_type))
  }
  
  res_tree$parents = pars
  res_tree$children = tree_from_parents(pars)
  
  return(res_tree)
}



#' Generates observations from a CMN network model
#' 
#' @param tm CMN network object
#' @param Nobs Number of edge matrices to generate
#' 
#' @return Array of generated edge matrices
#' 
#' @export
#' 
sample_network_tree = function(tm, Nobs = 1) {
  # Written by Andrew Thomas
  
  nn <- tm$nodes
  
  clo.anc <- closest_ancestor(tm)$anc.table
  out <- array(0, c(nn, nn, Nobs))
  series <- lower_diag(nn)
  for (kk in 1:Nobs) {
    out[series+nn^2*(kk-1)] <- rbinom(nn*(nn-1)/2, 1, tm$prob[clo.anc[series]-nn])
    out[,,kk] <- out[,,kk]+t(out[,,kk])
  }
  return(out)
}



# Latent Space Models -----------------------------------------------------


#' Generate model parameters for a latent space model. 
#' 
#' @param Nnodes Number of nodes
#' @param model_param List of model parameters
#' 
#' @return Latent space model object
#' 
#' @export
#' 
sample_model_latent = function(Nnodes, model_param = set_model_param()) {
  K = model_param$latent_nclass
  D = model_param$latent_dim
  gen_norm = model_param$latent_isgennorm
  sd_center = model_param$latent_sdcenter
  ## TODO: [TEMP] REMOVE these variables. 
  
  # K,gen_norm = TRUE, D = 2, sd_center = 5) {
  
  group_assign = sample(1:K, size = Nnodes, replace = TRUE)
  centers = list()
  
  for(j in 1:K) {
    if (gen_norm) {
      centers[[j]] = rnorm(n = D, mean = 0, sd = sd_center)
    } else {
      centers[[j]] = runif(n = D, min = -2 * sd_center, max = 2 * sd_center)
    }
  }
  
  locs = matrix(nrow= Nnodes, ncol = D)
  for(j in 1:Nnodes) {
    locs[j,] = centers[[group_assign[j]]] + rnorm(n = D, mean = 0, sd = 1)
  }
  return(list(locs = locs, alpha = runif(1, 0, 10)))
}


#' Simulate networks from a latent space model object. 
#' 
#' @param lm Latent space model object
#' @param Nobs Number of networks drawn
#' 
#' @return Array of networks
#' 
#' @export
#' 
sample_network_latent = function(lm, Nobs = 1) {
  Nnodes = nrow(lm$locs)
  odds = exp(lm$alpha - dist(lm$locs))
  prob_mat = matrix(nrow = Nnodes, ncol = Nnodes)
  prob_mat[lower.tri(prob_mat)] = odds / (1 + odds)
  
  res = array(0, dim = c(Nnodes, Nnodes, Nobs))
  for(d in 1:Nobs) {
    for(k in 1:Nnodes) {
      for(j in 1:Nnodes) {
        if (k < j) {
          v = rbinom(n = 1, size = 1, prob = prob_mat[j,k])
          res[j,k,d] = v
          res[k,j,d] = v
        }
      }}
  }
  return(res) 
}

