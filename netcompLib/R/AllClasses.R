# Class Definitions -------------------------------------------------------
## TODO: Get default print/summary methods

setClass("NetworkModel", representation(Nnodes = "numeric"))
setClass("NetworkModelSBM", representation(assign = "numeric", probmat = "matrix"), contains = "NetworkModel")
setClass("NetworkModelHRG", representation(parents = "numeric", children = "list", prob = "numeric"), contains = "NetworkModel")
setClass("NetworkModelLSM", representation(locs = "matrix", alpha = "numeric"), contains = "NetworkModel")
setClass("NetworkModelRND", representation(counts = "numeric", prob = "numeric", ids = "list"), contains = "NetworkModel")
setClass("NetworkModelPair", representation(m1 = "NetworkModel", m2 = "NetworkModel", is_null = "logical"), contains = "NetworkModel")

setClass("NetworkStruct", representation(Nnodes = "numeric"))
setClass("NetworkStructSBM", representation(groups = "numeric", counts = "numeric", expand = "list", correct = "numeric"), contains = "NetworkStruct")
setClass("NetworkStructRND", representation(counts = "numeric", ids = "list"), contains = "NetworkStruct")
setClass("NetworkStructHRG", representation(tree_list = "list", expand = "list", counts = "numeric"), contains = "NetworkStruct")
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
  Nnodes = model_params$Nnodes
  
  ## Helper function for adjusting block model probabilities  
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
  
  # if block assignments are not pre-specified: 
  K = model_params$block_nclass
  if (is.null(model_params$block_assign)) {
    group_assign = sample(1:K, size = Nnodes, replace = TRUE)
  } else { # Use prespecified block assignments
    group_assign = model_params$block_assign
  }
  
  # if block probability matrix is not pre-specified: 
  if (is.null(model_params$block_probs)) {
    # If we want to control average density: 
    if (!is.null(model_params$block_avgdensity)) {
      prob_matrix = adjust_blockprobs(mod = res, avgden = model_params$block_avgdensity, plimit = c(model_params$pmin, model_params$pmax))
    } else {
      prob_matrix = matrix(0, nrow = K, ncol = K)
      for(j in 1:K) { for(i in 1:K) {
        if (i <= j) {
          prob_matrix[i,j] = runif(1, model_params$pmin, model_params$pmax)
          if (i != j) {prob_matrix[j,i] = prob_matrix[i,j]}
        }
      }}
    }
  } else { # Use prespecified block probabilities
    prob_matrix = model_params$block_probs
  }
  
  netm = new("NetworkModelSBM", Nnodes = Nnodes, assign = group_assign, probmat = prob_matrix)
  return(netm)
}



#' Constructor for HRG network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_params A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelHRG object
#' 
#' @export
#' 
NetworkModelHRG = function(Nnodes = NULL, model_params = set_model_param()) {
  if (is.null(Nnodes)) { Nnodes = model_params$Nnodes }
  
  ## TODO: - fill this in eventually 
  
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
  res_tree$prob = runif(Nnodes-1, min = random_plimit[1], max = random_plimit[2])
  
  netm = new("NetworkModelHRG", Nnodes = Nnodes, parents = res_tree$parents, 
             children = res_tree$children, prob = res_tree$prob)
  return(netm)
}




#' Constructor for LSM network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_params A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelLSM object
#' 
#' @export
#' 
NetworkModelLSM = function(Nnodes = NULL, model_params = set_model_param()) {
  if (is.null(Nnodes)) { Nnodes = model_params$Nnodes }
  
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
#' @param Nnodes Number of nodes in the network model
#' @param model_params A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkModelRND = function(Nnodes = NULL, model_params = set_model_param()) {
  if (is.null(Nnodes)) { Nnodes = model_params$Nnodes }
  
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
#' @param m1 temp
#' @param m2 temp
#' @param is_null temp
#' 
#' @return Object of class NetworkModelPair
#' 
#' @export
#' 
NetworkModelPair = function(m1, m2 = NULL, is_null = FALSE) {
  ## TODO: Make m1/m2 more extensible (=> can input param list)
  if (is_null) {
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m1, is_null = TRUE)
  } else {
    if (is.null(m2)) { stop("---You must provide the second model if the null hypothesis is FALSE.---") }
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m2, is_null = FALSE)
  }
  return(netmp)
}




# Constructors -- NetworkStruct & subclasses ------------------------------

#' Instantiates an object of class NetworkStruct
#' 
#' @param Nnodes Number of nodes in network
#' @param type Type of network model (accepts 'block', 'tree', 'random', 'latent') ('none' is also accepted, but the resulting object isn't really interesting; its mainly for testing the class)
#' @param model_params Model parameters specified by set_model_param()
#' 
#' @return Object of class NetworkStruct
#' 
#' @export
#' 
NetworkStruct = function(Nnodes = NULL, type = NULL, model_params = set_model_param()) {
  # creates a default NetworkStruct object -- this is the default constructor, but probably should never be used. the specific ones for a specific model should be used. 
  # maybe want to make this eventually call the network generation methods
  
  if (is.null(type)) { type = model_params$type }
  ## currently this does nothing but return a lame object, that satisfies the generic methods (although the generic methods would not do much?)
  if (type == "none") { return(new("NetworkStruct", Nnodes = Nnodes)) }
  if (type == "block") { return(NetworkStructSBM(Nnodes = Nnodes, model_params = model_params)) }
  if (type == "tree") { return(NetworkStructHRG(Nnodes = Nnodes, model_params = model_params)) }
  if (type == "latent") { stop("Not valid with latent space models") }
  if (type == "random") { return(NetworkStructRND(Nnodes = Nnodes, model_params = model_params)) }
  stop("Invalid 'type' specified")
}




#' Constructor for class NetworkStructList -- a list of network model structure information. 
#' 
#' These are the random edge partitions used in the network hypothesis testing. 
#' 
#' @param Nnodes Number of nodes for each network
#' @param Nmodels Number of random edge partitions to generate
#' @param type Type of random edge partition to generate
#' @param model_params List of model parameters; see set_model_param(). 
#' 
#' @return Object of class NetworkStructList
#' 
#' @export
#' 
NetworkStructList = function(Nnodes = NULL, Nmodels = 10, type = NULL, model_params = set_model_param()) {
  # TODO: [Improvement] 'type' can be a vector, and call it vectorized => resulting network list has multiple types
  if (is.null(Nnodes)) { Nnodes = model_params$Nnodes }
  
  res = replicate(n = Nmodels, expr = NetworkStruct(Nnodes = Nnodes, type = type, model_params = model_params))
  netsl = new("NetworkStructList", Nnodes = Nnodes, models = res)
  return(netsl)
}



#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_params A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkStructRND = function(Nnodes = NULL, model_params = set_model_param()) {
  if (is.null(Nnodes)) { Nnodes = model_params$Nnodes }
  
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelRND(Nnodes = Nnodes, model_params = model_params)
  return(extractStruct(NetM))
}


#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_params A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkStructHRG = function(Nnodes = NULL, model_params = set_model_param()) {
  if (is.null(Nnodes)) { Nnodes = model_params$Nnodes }
  
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelHRG(Nnodes = Nnodes, model_params = model_params)
  return(extractStruct(NetM))
}


#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_params A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelSBM object
#' 
#' @export
#' 
NetworkStructSBM = function(Nnodes = 10, model_params = set_model_param()) {
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelSBM(Nnodes = Nnodes, model_params = model_params)
  return(extractStruct(NetM))
}

