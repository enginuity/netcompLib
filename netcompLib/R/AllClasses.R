# Class Definitions -------------------------------------------------------

setClass("NetworkModel", representation(Nnodes = "numeric"))
setClass("NetworkModelSBM", representation(groups = "numeric", probmat = "matrix"), contains = "NetworkModel")
setClass("NetworkModelHRG", representation(parents = "numeric", children = "list", prob = "numeric"), contains = "NetworkModel")
setClass("NetworkModelLSM", representation(locs = "matrix", alpha = "numeric"), contains = "NetworkModel")
setClass("NetworkModelRND", representation(counts = "numeric", prob = "numeric", ids = "list"), contains = "NetworkModel")
setClass("NetworkModelPair", representation(m1 = "NetworkModel", m2 = "NetworkModel", is_null = "logical", model_type = "character", addl_param = "list"), contains = "NetworkModel")

setClass("NetworkStruct", representation(Nnodes = "numeric"))
setClass("NetworkStructSBM", representation(groups = "numeric", counts = "numeric", expand = "list", correct = "numeric"), contains = "NetworkStruct")
setClass("NetworkStructRND", representation(counts = "numeric", ids = "list"), contains = "NetworkStruct")
setClass("NetworkStructHRG", representation(parents = "numeric", children = "list", expand = "list", counts = "numeric"), contains = "NetworkStruct")
setClass("NetworkStructList", representation(models = "list"), contains = "NetworkStruct")



# Constructors -- NetworkModel & subclasses -------------------------------

#' Instantiates an object of class NetworkModel
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' 
#' @return [NetworkModel] :: A representation of a model generated as specified
#' 
#' @export
#' 
NetworkModel = function(model_params = set_model_param()) {
  type = model_params$type
  if (type == "none") { return(new("NetworkModel")) }
  if (type == "block") { return(NetworkModelSBM(model_params)) }
  if (type == "tree") { return(NetworkModelHRG(model_params)) }
  if (type == "latent") { return(NetworkModelLSM(model_params)) }
  if (type == "random") { return(NetworkModelRND(model_params)) }
  stop("Invalid 'type' specified")
}



#' Constructor for SBM network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' 
#' @return [NetworkModelSBM] :: A representation of the generated model
#' 
#' @export
#' 
NetworkModelSBM = function(model_params = set_model_param()) {
  ########## Helper Functions ##########
  
  ## Helper function -- Reindexes group assignments to remove empty groups! 
  doGroupReassignment = function(groups, pmat) {
    ## Given input group IDs and probability matrix, reindex group IDs if necessary. 
    
    ## Compute new group assignments
    ordered = sort(unique(groups), decreasing = FALSE)
    reassignment = match(seq_len(max(groups)), ordered)
    
    ## Reassign group IDs
    new_groups = reassignment[groups]
    rc_keep = which(!is.na(reassignment))
    new_pmat = pmat[rc_keep, rc_keep, drop = FALSE]
    
    return(list(groups = new_groups, pmat = new_pmat))
  }
  
  ## Helper function -- Generates/adjusts block model probabilities  
  generateBlockProbs = function(groups, avgden = NULL, plims = c(0,1)) {
    ## This function, when given class assignments in 'groups', forces the probability matrix to take on values within 'plims' AND keep an average overall density given by 'avgden'. 
    
    classct = table(groups)
    NC = max(groups) # Number of classes (NC)
    if (NC < 1) { stop("Bad group assigments -- max group assignment < 1") }
    Nparams = NC * (NC - 1) / 2 + NC
    
    ## Stores assigned probabilities
    coordmat = matrix(NA, nrow = Nparams, ncol = 4)
    colnames(coordmat) = c("row", "col", "count", "prob")
    
    ## Fill out dyad group counts
    II = 1
    for(j in 1:NC) { for(k in j:NC) {
      coordmat[II,1:3] = c(j, k, ifelse(j==k, classct[j] * (classct[j]-1)/2, classct[j] * classct[k]))
      II = II + 1
    }}
    
    ## Assign probabilities
    sampProbVec = function(counts, avgden, plims) {
      ## samples vectors so that the average density matches 'avgden'. 
      valid = FALSE
      while(!valid) {
        testvec = runif(length(counts), plims[1], plims[2])
        j = sample(seq_along(counts), size = 1) ## Pick random index of dyad group to change probability of
        testvec[j] = (sum(counts) * avgden - sum(counts[-j] * testvec[-j]))/counts[j]
        if (testvec[j] > plims[1] & testvec[j] < plims[2]) { valid = TRUE }
      }
      return(testvec)
    }
    
    if (is.null(avgden)) {
      coordmat[,4] = runif(n = Nparams, min = plims[1], max = plims[2])
    } else {
      coordmat[,4] = sampProbVec(coordmat[,3], avgden, plims)
    }
    
    ## Fill out in matrix form
    pmat = matrix(NA, ncol = NC, nrow = NC)
    for(i in 1:Nparams) {
      pmat[coordmat[i,1], coordmat[i,2]] = coordmat[i,4]
      pmat[coordmat[i,2], coordmat[i,1]] = coordmat[i,4]
    }
    
    return(pmat)
  }
  
  ########## Beginning of constructor function ##########
  Nnodes = model_params$Nnodes
  
  ## Figure out group assignment - uniform sample or use input assignment
  if (is.null(model_params$block_assign)) {
    group_assign = sample(1:model_params$block_nclass, size = Nnodes, replace = TRUE)
  } else {
    group_assign = model_params$block_assign
  }
  
  ## Specify block probability matrix
  if (is.null(model_params$block_probs)) {
    prob_matrix = generateBlockProbs(group_assign, model_params$block_avgdensity, c(model_params$pmin, model_params$pmax))
  } else {
    prob_matrix = model_params$block_probs
  }
  
  ## Remove any empty group assignments and relabel groups as necessary
  temp = doGroupReassignment(group_assign, prob_matrix)
  group_assign = temp$groups; prob_matrix = temp$pmat
  
  ## Create network model object
  netm = new("NetworkModelSBM", Nnodes = Nnodes, groups = group_assign, probmat = prob_matrix)
  return(netm)
}


#' Constructor for HRG network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' 
#' @return [NetworkModelHRG] :: A representation of the generated model
#' 
#' @export
#' 
NetworkModelHRG = function(model_params = set_model_param()) {
  ## TODO: [Rewrite] this
  
  Nnodes = model_params$Nnodes
  
  # helper function that generates a fixed structure tree (as close to binary tree as possible)
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
                   Nnodes=Nnodes)
    
    return(output)
  }
  
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
  
  ## TODO: [TEMP] Using these two variables shouldn't be necessary. 
  random_type = model_params$tree_type
  random_plimit = c(model_params$pmin, model_params$pmax)
  
  res_tree = list()
  
  if (random_type == "original") {
    res_tree = starter_tree(Nnodes)  
    pars = res_tree$parents
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
  res_tree$children = HRG_treeFromParents(pars)
  res_tree$prob = runif(Nnodes-1, min = random_plimit[1], max = random_plimit[2])
  
  netm = new("NetworkModelHRG", Nnodes = Nnodes, parents = res_tree$parents, 
             children = res_tree$children, prob = res_tree$prob)
  return(netm)
}


#' Constructor for LSM network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' 
#' @return [NetworkModelLSM] :: A representation of the generated model
#' 
#' @export
#' 
NetworkModelLSM = function(model_params = set_model_param()) {
  ## TODO: [Rewrite] this
  
  Nnodes = model_params$Nnodes
  
  K = model_params$latent_nclass
  D = model_params$latent_dim
  gen_norm = model_params$latent_isgennorm
  sd_center = model_params$latent_sdcenter
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
  
  netm = new("NetworkModelLSM", Nnodes = Nnodes, locs = locs, alpha = runif(1,0,10))
  return(netm)
}



#' Constructor for RND network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' 
#' @return [NetworkModelRND] :: A representation of the generated model
#' 
#' @export
#' 
NetworkModelRND = function(model_params = set_model_param()) {
  ## TODO: [Rewrite] this
  Nnodes = model_params$Nnodes
  
  rnd_Ngroups = model_params$random_ngroups
  
  # Compute index numbers of the lower-diagonal portion of the matrix
  idm = matrix(1:(Nnodes^2), nrow = Nnodes)
  idm[lower.tri(x = idm, diag = TRUE)] = 0
  idm = as.vector(idm)
  good_ids = idm[idm != 0]
  if (length(good_ids) < rnd_Ngroups) { stop("Too many groups for too few edges") }
  
  # Sample ids to assign into groups
  rand_order = sample(good_ids, size = length(good_ids), replace = FALSE)
  sizes = rep(floor(length(good_ids) / rnd_Ngroups), times = rnd_Ngroups)
  m = length(good_ids) %% rnd_Ngroups
  if (m > 0) {
    sizes[1:m] = sizes[1:m] + 1
  }
  
  counted = 0
  id_list = list()
  for(k in 1:rnd_Ngroups) {
    id_list[[k]] = sort(rand_order[seq(from = counted+1, length.out = sizes[k])])
    counted = counted + sizes[k]
  }
  
  netm = new("NetworkModelRND", Nnodes = Nnodes, counts = sizes, ids = id_list, 
             prob = runif(n = rnd_Ngroups, min = model_params$pmin, max = model_params$pmax))
  return(netm)
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (NetworkModelPair)
#' Instantiates an object of class NetworkModelPair
#' 
#' This is done by providing both network models. If m2 is not given and is_null is FALSE, there is an error. Otherwise, m2 will be ignored if is_null is TRUE. (ie null hypothesis means both models are the same, so m1)
#' There is also no check that m1 and m2 are on the same network size, but they should be. -- add this in...
#' 
#' @param m1 [\code{\link{NetworkModel}} OR list] :: This can either be a model object, or a list containing model parameters
#' @param m2 [\code{\link{NetworkModel}} OR list] :: This can either be a model object, or a list containing model parameters
#' @param is_null [logical] :: If TRUE, m2 is ignored (since the null hypothesis means that both models are identical)
#' @param model_type [char] :: Can be 'default', 'correlated', or 'densitydiff'
#' @param addl_param [list] :: Additional parameters for correlated or densitydiff models
#' 
#' @return [NetworkModelPair] :: An object containing information about two network models. 
#' 
#' @export
#' 
NetworkModelPair = function(m1, m2 = NULL, is_null = FALSE, model_type = "default", addl_param = list()) {
  if (!(is(m1, "NetworkModel"))) { m1 = NetworkModel(m1) }
  if (!(is(m2, "NetworkModel")) & !is.null(m2)) { m2 = NetworkModel(m2) }
  
  if (model_type == "default") {
    if (is_null) {
      netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m1, is_null = TRUE, model_type = model_type, addl_param = addl_param)
    } else {
      if (is.null(m2)) { stop("---You must provide the second model if the null hypothesis is FALSE.---") }
      netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m2, is_null = FALSE, model_type = model_type, addl_param = addl_param)
    }
    
  } else if (model_type == "densitydiff") {
    ## In this case, null hypothesis must be true, as the second model must have same structure AND be exactly a fixed density parameter away. 
    
    if (!("dd_param_add" %in% names(addl_param))) { addl_param$dd_param_add = rnorm(1) }
    m2 = m1
    origprobs = get_dyadgroup_prob(m1)
    m2 = set_dyadgroup_prob(
      m2, origprobs$names, faraway::ilogit(faraway::logit(origprobs$probs) + addl_param$dd_param_add))
    
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m2, is_null = FALSE, model_type = model_type, addl_param = addl_param)
  } else if (model_type == "correlated") {
    
    ## In this case, the structure should be the same, but can either be null or non-null. 
    ## If null -> model must be the same, so m2 is ignored. 
    if (is_null) { 
      m2 = m1 
    } else if (any(getEdgeProbMat(m1, 'group') != getEdgeProbMat(m2, 'group'))) {
      stop("Invalid second model -- dyad partition is not the same.")
    } 
    
    ## Assign correlation parameter if not inputted
    if (!("c_param_corr" %in% names(addl_param))) { addl_param$c_param_corr = rnorm(1) }
    
    ## Check if remaining parameters are already assigned; assign if not. 
    if (!("c_param_a" %in% names(addl_param))) { 
      temp = aggstat_single(m1, getEdgeProbMat(m1)) 
      addl_param$c_param_a = temp$x / temp$n
    }
    
    if (!("c_param_b" %in% names(addl_param))) { 
      temp = aggstat_single(m2, getEdgeProbMat(m2)) 
      addl_param$c_param_b = temp$x / temp$n
    }
    
    if (!("c_names" %in% names(addl_param))) {
      addl_param$c_names = aggstat_single(m1, getEdgeProbMat(m1))$names
    }
    
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m2, is_null = FALSE, model_type = model_type, addl_param = addl_param)
  } else {
    stop("Invalid model_type")
  }
  
  return(netmp)
}




# Constructors -- NetworkStruct & subclasses ------------------------------

#' Instantiates an object of class NetworkStruct
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' @param NetM [\code{\link{NetworkModel}}] :: If not NULL, this model is used, and its structure extracted
#' 
#' @return [NetworkStruct] :: A representation of the generated model structure
#' 
#' @export
#' 
NetworkStruct = function(model_params = set_model_param(), NetM = NULL) {
  type = model_params$type
  
  if (type == "none") { return(new("NetworkStruct")) }
  if (type == "block") { return(NetworkStructSBM(model_params, NetM)) }
  if (type == "tree") { return(NetworkStructHRG(model_params, NetM)) }
  if (type == "latent") { stop("Not valid with latent space models") }
  if (type == "random") { return(NetworkStructRND(model_params, NetM)) }
  stop("Invalid 'type' specified")
}


#' Constructor for class NetworkStructList -- a list of network model structure information. 
#' 
#' These are the random edge partitions used in the network hypothesis testing. 
#' 
#' @param Nmodels [int] :: Number of random edge partitions to generate
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' @param NetMPair [\code{\link{NetworkModelPair}}] :: If not NULL, these models are used, and its structures extracted
#' 
#' @return [NetworkStructList] :: A representation of the generated model structure
#' 
#' @export
#' 
NetworkStructList = function(Nmodels = 10, model_params = set_model_param(), NetMPair = NULL) {
  
  if (is.null(NetMPair)) {
    res = replicate(n = Nmodels, expr = NetworkStruct(model_params))
  } else {
    res = replicate(n = 2, expr = NetworkStruct(model_params))
    res@models[[1]] = NetworkStruct(NetM = NetM@m1)
    res@models[[2]] = extractStruct(NetM = NetM@m2)
  }
  
  return(new("NetworkStructList", Nnodes = model_params$Nnodes, models = res))
}


#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' @param NetM [\code{\link{NetworkModelRND}}] :: If not NULL, this model is used, and its structure extracted
#' 
#' @return [NetworkStructRND] :: A representation of the generated model structure
#' 
#' @export
#' 
NetworkStructRND = function(model_params = set_model_param(), NetM = NULL) {
  if (is.null(NetM)) { NetM = NetworkModelRND(model_params) }
  
  nets = new("NetworkStructRND", Nnodes = getNnodes(NetM), counts = NetM@counts, ids = NetM@ids)
  return(nets)
}


#' Constructor for HRG network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' @param NetM [\code{\link{NetworkModelHRG}}] :: If not NULL, this model is used, and its structure extracted
#' 
#' @return [NetworkStructHRG] :: A representation of the generated model structure
#' 
#' @export
#' 
NetworkStructHRG = function(model_params = set_model_param(), NetM = NULL) {
  if (is.null(NetM)) { NetM = NetworkModelHRG(model_params) }
  expc = HRG_expandedChildren(NetM)
  
  nets = new("NetworkStructHRG", Nnodes = getNnodes(NetM), parents = NetM@parents, children = NetM@children, 
             expand = expc, counts = sapply(expc, function(x) {length(x[[1]]) * length(x[[2]]) }))
  return(nets)
}


#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param model_params [list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters
#' @param NetM [\code{\link{NetworkModelSBM}}] :: If not NULL, this model is used, and its structure extracted
#' 
#' @return [NetworkStructSBM] :: A representation of the generated model structure
#' 
#' @export
#' 
NetworkStructSBM = function(model_params = set_model_param(), NetM = NULL) {
  if (is.null(NetM)) { NetM = NetworkModelSBM(model_params) }
  
  # group assignments
  ga = NetM@groups
  NClass = length(unique(ga))
  
  counts = rep(0, times = NClass + NClass * (NClass - 1) / 2)
  correction = rep(0, times = NClass + NClass * (NClass - 1) / 2)
  
  expanded = list()
  cur = 1
  for(k in 1:NClass) { for(m in 1:NClass) {
    if (k >= m) {
      if (k == m) { 
        counts[cur] = sum(ga == k) * (sum(ga == k) - 1) / 2 
        correction[cur] = 2
      } else {
        counts[cur] = sum(ga == k) * sum(ga == m)
        correction[cur] = 1
      }
      expanded[[cur]] = list(which(ga == k), which(ga == m))
      cur = cur + 1
    }
  }}
  
  nets = new("NetworkStructSBM", Nnodes = getNnodes(NetM), groups = ga, counts = counts, expand = expanded, correct = correction)
  return(nets)
}

