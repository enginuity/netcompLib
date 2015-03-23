## TODO: Reorganize this code, but figure out how to do it to not make it fail ?!
# Class Definitions -------------------------------------------------------

setClass("NetworkModel", representation(Nnodes = "numeric"))
setClass("NetworkModelSBM", representation(assign = "numeric", probmat = "matrix"), contains = "NetworkModel")
setClass("NetworkModelHRG", representation(parents = "numeric", children = "list", prob = "numeric"), contains = "NetworkModel")
setClass("NetworkModelLSM", representation(locs = "matrix", alpha = "numeric"), contains = "NetworkModel")
setClass("NetworkModelPair", representation(m1 = "NetworkModel", m2 = "NetworkModel", is_null = "logical"), contains = "NetworkModel")

setClass("NetworkStruct", representation(Nnodes = "numeric"))
setClass("NetworkStructSBM", representation(groups = "numeric", counts = "numeric", expand = "list", correct = "numeric"), contains = "NetworkStruct")
setClass("NetworkStructRND", representation(counts = "numeric", ids = "list"), contains = "NetworkStruct")
setClass("NetworkStructHRG", representation(tree_list = "list", expand = "list", counts = "numeric"), contains = "NetworkStruct")
setClass("NetworkStructList", representation(models = "list"), contains = "NetworkStruct")


# Generic Definitions -----------------------------------------------------

setGeneric("getNnodes", function(NetM) standardGeneric("getNnodes"))
setGeneric("sampleNetwork", function(NetM, Nobs = 1, ...) standardGeneric("sampleNetwork"))
setGeneric("computePval", function(NetS, adja1, adja2, Nobs, pl, mode) standardGeneric("computePval"))
setGeneric("extractStruct", function(NetM) standardGeneric("extractStruct"))
setGeneric("getEdgeProbMat", function(NetM) standardGeneric("getEdgeProbMat"))
setGeneric("getNetType", function(NetM) standardGeneric("getNetType"))


# Method Definitions ------------------------------------------------------

# This file defines the NetworkModel class.

## TODO: Add validity checks?



# Defines the NetworkModel Class and its subclasses -----------------------


#' Instantiates an object of class NetworkModel
#' 
#' @param Nnodes Number of nodes in network
#' @param type Type of network model (accepts 'block', 'tree', 'random', 'latent') ('none' is also accepted, but the resulting object isn't really interesting; its mainly for testing the class)
#' @param model_param Model parameters specified by set_model_param()
#' 
#' @return Object of class NetworkModel
#' 
#' @export
#' 
NetworkModel = function(Nnodes = 10, type = "none", model_param = set_model_param()) {
  # creates a default networkmodel object -- this is the default constructor, but probably should never be used. the specific ones for a specific model should be used. 
  # maybe want to make this eventually call the network generation methods
  
  ## currently this does nothing but return a lame object, that satisfies the generic methods (although the generic methods would not do much?)
  if (type == "none") { return(new("NetworkModel", Nnodes = Nnodes)) }
  if (type == "block") { return(NetworkModelSBM(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "tree") { return(NetworkModelHRG(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "latent") { return(NetworkModelLSM(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "random") { return(NetworkModelRND(Nnodes = Nnodes, model_param = model_param)) }
  stop("Invalid 'type' specified")
}


# Define Generic Methods --------------------------------------------------


### New method to extract the number of nodes from a network


#' Generic Function -- Extracts the number of nodes in the network model
#' 
#' @param NetM NetworkModel object
#' 
#' @return Numeric -- Number of nodes in network model
#' 
#' @export
#' 
getNnodes = function(NetM) { NULL }

#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkModel
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkModel = function(NetM) {
  # NetM should be object of type NetworkModel
  return(NetM@Nnodes)
}






#' Generic -- Samples an observed networks from a given network model(s)
#' 
#' Nobs is the number of network observations to simulate per simulation. This is built in here since callling sampleNetwork many times is bad (since it requires re-computing the edge probability matrix, which is the same each time given the same network model). There are several cases: 
#' If Nsim = 1, Nobs = 1 -> Result is a matrix
#' If Nsim = 1, Nobs > 1 -> Result is an array (with third dimension equal to Nobs)
#' If Nsim > 1 -> Result is a list of matrices or arrays (depends on Nobs)
#' 
#' @param NetM NetworkModel object
#' @param Nobs Number of network observations to simulate
#' @param Nsim Number of simulations
#' 
#' @return List or array of adjacency matrices
#' 
#' @export
#' 
sampleNetwork = function(NetM, Nobs = 1, Nsim = 1) { 
  return(NULL) 
}




#' Samples any number of observed networks from a given network model
#' 
#' Note -- this implementation requires the edge probability matrix. It's not really efficient (since it re-computes the edge probability matrix each time)
#' 
#' @param NetM Network Model object
#' @param Nobs Number of observed networks to sample
#' @param Nsim Number of simulations (number of times to do Nobs samples)
#' 
#' @return List or array of simulated networks
#' 
#' @export
#' 
sampleNetwork.NetworkModel = function(NetM, Nobs = 1, Nsim = 1) { 
  #Nsim -- if you want multiple sets of networks from the same model, generate it all now; will be given in a list if Nsim > 1. 
  
  # TODO: Use this style of code to improve on network simulation
  #   # Written by Andrew Thomas
  #   
  #   nn <- tm$nodes
  #   
  #   clo.anc <- closest_ancestor(tm)$anc.table
  #   out <- array(0, c(nn, nn, Nobs))
  #   series <- lower_diag(nn)
  #   for (kk in 1:Nobs) {
  #     out[series+nn^2*(kk-1)] <- rbinom(nn*(nn-1)/2, 1, tm$prob[clo.anc[series]-nn])
  #     out[,,kk] <- out[,,kk]+t(out[,,kk])
  #   }
  
  # uses variables out of scope (epmat, Nnodes)
  gen_one_network = function() {
    res = matrix(0, nrow = Nnodes, ncol = Nnodes) 
    for(k in 1:Nnodes) { for(j in 1:Nnodes) {
      if (k < j) {
        v = rbinom(n = 1, size = 1, prob = epmat[k,j])
        res[j,k] = v
        res[k,j] = v
      }
    }}
    return(res)
  }
  
  Nnodes = getNnodes(NetM)
  
  epmat = getEdgeProbMat(NetM)
  if (Nsim == 1) {
    
    res = array(0, dim = c(Nnodes, Nnodes, Nobs))
    for(d in 1:Nobs) {
      res[,,d] = gen_one_network()
    }
    return(res)
    
  } else {
    
    reslist = list()
    for(f in 1:Nsim) {
      res = array(0, dim = c(Nnodes, Nnodes, Nobs))
      for(d in 1:Nobs) {
        res[,,d] = gen_one_network()
      }
      reslist[[f]] = res
    }
    return(reslist)
    
  }
}




# Define Default print/summary methods ------------------------------------





## Methods for network-model specific information
## These are the functions that need to be defined for each specific network model. 


# Define generic functions that are model-specific ------------------------




#' Generic Function -- Extracts the type of network model
#' 
#' Allowable types are currently: "block", "tree", "latent", "random"
#' 
#' @param NetM Object of class NetworkModel (or inherits this type)
#' 
#' @return Character -- type of network model
#' 
#' @export
#' 
getNetType = function(NetM) { 
  NULL 
}






#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM Network model object
#' 
#' @return Matrix -- the edge probability matrix
#' 
#' @export
#' 
getEdgeProbMat = function(NetM) {
  return(NULL) 
}






#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct = function(NetM) {
  return(NULL) 
}



# Define generic null-model values so this doesn't crash ------------------

#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkModel = function(NetM) {
  "none" 
}


#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM Network Model object
#' 
#' @return NULL - no edge probability matrix for generic network model object
#' 
#' @export
#' 
getEdgeProbMat.NetworkModel = function(NetM) {
  NULL 
}


#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModel = function(NetM) {
  stop("This usage case is not defined. ")
}


# This file defines the NetworkStruct class.

## TODO: Add validity checks?



# Defines the NetworkStruct Class and its subclasses -----------------------

#' Instantiates an object of class NetworkStruct
#' 
#' @param Nnodes Number of nodes in network
#' @param type Type of network model (accepts 'block', 'tree', 'random', 'latent') ('none' is also accepted, but the resulting object isn't really interesting; its mainly for testing the class)
#' @param model_param Model parameters specified by set_model_param()
#' 
#' @return Object of class NetworkStruct
#' 
#' @export
#' 
NetworkStruct = function(Nnodes = 10, type = "none", model_param = set_model_param()) {
  # creates a default NetworkStruct object -- this is the default constructor, but probably should never be used. the specific ones for a specific model should be used. 
  # maybe want to make this eventually call the network generation methods
  
  ## currently this does nothing but return a lame object, that satisfies the generic methods (although the generic methods would not do much?)
  if (type == "none") { return(new("NetworkStruct", Nnodes = Nnodes)) }
  if (type == "block") { return(NetworkStructSBM(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "tree") { return(NetworkStructHRG(Nnodes = Nnodes, model_param = model_param)) }
  if (type == "latent") { stop("Not valid with latent space models") }
  if (type == "random") { return(NetworkStructRND(Nnodes = Nnodes, model_param = model_param)) }
  stop("Invalid 'type' specified")
}


# Define Generic Methods --------------------------------------------------

#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkStruct
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkStruct = function(NetM) {
  # NetM should be object of type NetworkStruct
  return(NetM@Nnodes)
}




# Define Default print/summary methods ------------------------------------





## Methods for network-model specific information
## These are the functions that need to be defined for each specific network model. 

#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkStruct = function(NetM) {
  "none" 
}





# Define new generics -----------------------------------------------------


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval)
#' <What does this function do>
#' 
#' @param NetS temp
#' @param adja1 temp
#' @param adja2 temp
#' @param Nobs temp
#' @param pl temp
#' @param mode temp
#' 
#' @return temp
#' 
#' @export
#' 
computePval = function(NetS, adja1, adja2, Nobs, pl, mode) {
  stop("Placeholder for documentation purposes")
}

#' Computes the p-values for each random partition in a NetworkStruct
#' 
#' @param NetS NetworkStruct to use (can be a NetworkStructList)
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of observations per input
#' @param pl parameter list, generated by set_sim_param()
#' @param mode Different modes -- 'default' just gives p-values; 'nodewise' gives chi-square contribs per node
#' 
#' @return List of arrays containing p-values
#' 
#' @export
#' 
computePval.NetworkStruct = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  stop("Not implemented for this case")
}

# This file defines the NetworkModelPair class.

## TODO: Add validity checks?

# Defines the NetworkModelPair Class and its subclasses -----------------------

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
  if (is_null) {
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m1, is_null = TRUE)
  } else {
    if (is.null(m2)) { stop("---You must provide the second model if the null hypothesis is FALSE.---") }
    netmp = new("NetworkModelPair", Nnodes = getNnodes(m1), m1 = m1, m2 = m2, is_null = FALSE)
  }
  return(netmp)
}


# Define Generic Methods --------------------------------------------------


#' Gets the number of nodes in the network model
#' 
#' @param NetM Object of class NetworkModelPair
#' 
#' @return Number of nodes
#' 
#' @export
#' 
getNnodes.NetworkModelPair = function(NetM) {
  # NetM should be object of type NetworkModelPair
  return(NetM@Nnodes)
}


#' Samples any number of observed networks from a given network model
#' 
#' Note -- this implementation requires the edge probability matrix. It's not really efficient (since it re-computes the edge probability matrix each time)
#' 
#' @param NetM Network Model object
#' @param Nobs Number of observed networks to sample
#' @param Nsim Number of simulations (number of times to do Nobs samples)
#' 
#' @return List of two lists/arrays of network models. 
#' 
#' @export
#' 
sampleNetwork.NetworkModelPair = function(NetM, Nobs = 1, Nsim = 1) { 
  return(list(sampleNetwork(NetM@m1, Nobs = Nobs, Nsim = Nsim), sampleNetwork(NetM@m2, Nobs = Nobs, Nsim = Nsim)))
}




# Define Default print/summary methods ------------------------------------





## Methods for network-model specific information
## These are the functions that need to be defined for each specific network model. 



# Define generic null-model values so this doesn't crash ------------------

#' Extracts the type of network model
#' 
#' @param NetM Network Model object
#' 
#' @return character 'none'
#' 
#' @export
#' 
getNetType.NetworkModelPair = function(NetM) {
  return(c(getNetType(NetM@m1), getNetType(NetM@m2)))
}


#' Extracts the edge probability matrix for a network model
#' 
#' @param NetM Network Model object
#' 
#' @return NULL - no edge probability matrix for generic network model object
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelPair = function(NetM) {
  return(list(getEdgeProbMat(NetM@m1),getEdgeProbMat(NetM@m2)))
}

#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelPair = function(NetM) {
  netsl = new("NetworkStructList", Nnodes = getNnodes(NetM), 
              models = list(extractStruct(NetM@m1), extractStruct(NetM@m2)))
  return(netsl)
}

# Defines the SBMNetworkModel class

#' Constructor for SBM network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelSBM object
#' 
#' @export
#' 
NetworkModelSBM = function(Nnodes = 10, model_param = set_model_param()) {
  
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
  K = model_param$block_nclass
  if (is.null(model_param$block_assign)) {
    group_assign = sample(1:K, size = Nnodes, replace = TRUE)
  } else { # Use prespecified block assignments
    group_assign = model_param$block_assign
  }
  
  # if block probability matrix is not pre-specified: 
  if (is.null(model_param$block_probs)) {
    # If we want to control average density: 
    if (!is.null(model_param$block_avgdensity)) {
      prob_matrix = adjust_blockprobs(mod = res, avgden = model_param$block_avgdensity, plimit = c(model_param$pmin, model_param$pmax))
    } else {
      prob_matrix = matrix(0, nrow = K, ncol = K)
      for(j in 1:K) { for(i in 1:K) {
        if (i <= j) {
          prob_matrix[i,j] = runif(1, model_param$pmin, model_param$pmax)
          if (i != j) {prob_matrix[j,i] = prob_matrix[i,j]}
        }
      }}
    }
  } else { # Use prespecified block probabilities
    prob_matrix = model_param$block_probs
  }
  
  netm = new("NetworkModelSBM", Nnodes = Nnodes, assign = group_assign, probmat = prob_matrix)
  return(netm)
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelSBM objects, this returns "block"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'block'
#' 
#' @export
#' 
getNetType.NetworkModelSBM = function(NetM) { "block" }


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model Object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelSBM = function(NetM) {
  res = matrix(0, nrow = NetM@Nnodes, ncol = NetM@Nnodes) 
  for(j in 1:(NetM@Nnodes-1)) { for(k in (j+1):NetM@Nnodes) { 
    res[j,k] = NetM@probmat[NetM@assign[j], NetM@assign[k]]
    res[k,j] = res[j,k]
  }}
  
  return(res)
}




#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelSBM = function(NetM) {
  # group assignments
  ga = NetM@assign
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

# Defines the NetworkModelLSM class

#' Constructor for LSM network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelLSM object
#' 
#' @export
#' 
NetworkModelLSM = function(Nnodes = 10, model_param = set_model_param()) {
  
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
  
  netm = new("NetworkModelLSM", Nnodes = Nnodes, locs = locs, alpha = runif(1,0,10))
  return(netm)
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelLSM objects, this returns "latent"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'latent'
#' 
#' @export
#' 
getNetType.NetworkModelLSM = function(NetM) { "latent" }



#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelLSM = function(NetM) {
  Nnodes = getNnodes(NetM)
  odds = exp(NetM@alpha - dist(NetM@locs))
  prob_mat = matrix(0, nrow = Nnodes, ncol = Nnodes)
  prob_mat[lower.tri(prob_mat)] = odds / (1 + odds)
  prob_mat = prob_mat + t(prob_mat)
  return(prob_mat) 
}




#' Constructor for class NetworkStructList -- a list of network model structure information. 
#' 
#' These are the random edge partitions used in the network hypothesis testing. 
#' 
#' @param Nnodes Number of nodes for each network
#' @param Nmodels Number of random edge partitions to generate
#' @param type Type of random edge partition to generate
#' @param model_param List of model parameters; see set_model_param(). 
#' 
#' @return Object of class NetworkStructList
#' 
#' @export
#' 
NetworkStructList = function(Nnodes = 10, Nmodels = 10, type = "none", model_param = set_model_param()) {
  # TODO: [Improvement] 'type' can be a vector, and call it vectorized => resulting network list has multiple types
  res = replicate(n = Nmodels, expr = NetworkStruct(Nnodes = Nnodes, type = type, model_param = model_param))
  netsl = new("NetworkStructList", Nnodes = Nnodes, models = res)
  return(netsl)
}


#' Extract the type of network model
#' 
#' @param NetM Object of class NetworkStructList
#' 
#' @return character vector of types
#' 
#' @export
#' 
getNetType.NetworkStructList = function(NetM) {
  return(sapply(NetM@models, getNetType))
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructList)
#' Computes the p-values for each random partition in a NetworkStructList
#' 
#' @param NetS NetworkStructList to use
#' @param adja1 Adjacency matrix/array
#' @param adja2 Adjacency matrix/array
#' @param Nobs Number of observations per input
#' @param pl parameter list, generated by set_sim_param()
#' @param mode temp
#' 
#' @return List of arrays containing p-values
#' 
#' @export
#' 
computePval.NetworkStructList = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  res = lapply(NetS@models, function(x) { computePval(x, adja1, adja2, Nobs, pl, mode = "default") } )
  return(res)
}


# Defines the NetworkModelHRG class

#' Constructor for HRG network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelHRG object
#' 
#' @export
#' 
NetworkModelHRG = function(Nnodes = 10, model_param = set_model_param()) {
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
  random_type = model_param$tree_type
  random_plimit = c(model_param$pmin, model_param$pmax)
  
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


#' Returns the type of network model
#' 
#' Specifically for NetworkModelHRG objects, this returns "tree"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'tree'
#' 
#' @export
#' 
getNetType.NetworkModelHRG = function(NetM) { "tree" }



#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelHRG = function(NetM) {
  nn = getNnodes(NetM)
  ## TODO: Rewrite closest_ancestor, so this section of code isn't this ugly. 
  
  tm = list(prob = NetM@prob, children = NetM@children, parents = NetM@parents, nodes = nn)
  clo.anc = closest_ancestor(tm)$anc.table
  
  out <- matrix(0, nn, nn)
  series <- lower_diag(nn)
  out[series] <- tm$prob[clo.anc[series]-nn]
  out = out + t(out)
  return(out)
}





#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelHRG = function(NetM) {
  tr = list(prob = NetM@prob, children = NetM@children, parents = NetM@parents, nodes = getNnodes(NetM))
  expc = expanded_children_from_tree(tr)
  
  nets = new("NetworkStructHRG", Nnodes = getNnodes(NetM), tree_list = tr, expand = expc, 
             counts = sapply(expc, function(x) {length(x[[1]]) * length(x[[2]]) }))
  return(nets)
}




# Defines the NetworkModelRND class
setClass("NetworkModelRND", representation(counts = "numeric", prob = "numeric", ids = "list"), contains = "NetworkModel")

#' Constructor for RND network model
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkModelRND = function(Nnodes = 10, model_param = set_model_param()) {
  rnd_Ngroups = model_param$random_ngroups
  
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
             prob = runif(n = rnd_Ngroups, min = model_param$pmin, max = model_param$pmax))
  return(netm)
}


#' Returns the type of network model
#' 
#' Specifically for NetworkModelRND objects, this returns "random"
#' 
#' @param NetM Network Model Object
#' 
#' @return character 'random'
#' 
#' @export
#' 
getNetType.NetworkModelRND = function(NetM) { "random" }
setMethod("getNetType", signature(NetM = "NetworkModelRND"), getNetType.NetworkModelRND)


#' Computes the edge probability matrix
#' 
#' @param NetM Network Model object
#' 
#' @return Edge probability matrix defined by model
#' 
#' @export
#' 
getEdgeProbMat.NetworkModelRND = function(NetM) {
  pmat = matrix(0, nrow = getNnodes(NetM), ncol = getNnodes(NetM))
  for(j in seq_along(NetM@counts)) {
    pmat[NetM@ids[[j]]] = NetM@prob[j]
  }
  return(pmat + t(pmat))
}
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelRND"), getEdgeProbMat.NetworkModelRND)


#' Extracts the Edge Structure from a Network Model
#' 
#' @param NetM network model obj
#' 
#' @return network struct obj
#' 
#' @export
#' 
extractStruct.NetworkModelRND = function(NetM) {
  nets = new("NetworkStructRND", Nnodes = getNnodes(NetM), counts = NetM@counts, ids = NetM@ids)
  return(nets)
}
setMethod("extractStruct", signature = (NetM = "NetworkModelRND"), extractStruct.NetworkModelRND)


# Defines the NetworkModelRND class

#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkStructRND = function(Nnodes = 10, model_param = set_model_param()) {
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelRND(Nnodes = Nnodes, model_param = model_param)
  return(extractStruct(NetM))
}



# Generic Function Definitions --------------------------------------------

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructRND)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructRND = function(NetM) {
  return("random")
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructRND)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructRND)
#' <What does this function do>
#' 
#' @param NetS temp
#' @param adja1 temp
#' @param adja2 temp
#' @param Nobs temp
#' @param pl temp
#' @param mode temp
#' 
#' @return temp
#' 
#' @export
#' 
computePval.NetworkStructRND = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  if (FALSE) {
    NetM = NetworkModel(Nnodes = 30, type = "block")
    adja1 = sampleNetwork(NetM)
    adja2 = sampleNetwork(NetM)
    Nobs = 1
    pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_models = c(1,20))
    NetS = NetworkStructRND(Nnodes = 30, model_param = set_model_param(block_nclass = 3))
    load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
  }
  
  ## Compute total count for each edge (and edgesumc is the sum of products, to be used in computing cell-wise correlations)
  if (Nobs > 1) {
    edgesum1 = apply(adja1, c(1,2), sum); edgesum2 = apply(adja2, c(1,2), sum); 
    edgesumc = apply(adja1*adja2, c(1,2), sum)
  } else if (Nobs == 1) {
    if (length(dim(adja1)) == 2) { edgesum1 = adja1 } else { edgesum1 = adja1[,,1] }
    if (length(dim(adja2)) == 2) { edgesum2 = adja2 } else { edgesum2 = adja2[,,1] }
    edgesumc = edgesum1 * edgesum2
  }
  
  ## Do counting
  vedgesum1 = as.vector(edgesum1); vedgesum2 = as.vector(edgesum2); vedgesumc = as.vector(edgesumc)
  obs1_count = sapply(NetS@ids, function(x) { sum(vedgesum1[x])} )
  obs2_count = sapply(NetS@ids, function(x) { sum(vedgesum2[x])} )
  obsc_count = obs1_count + obs2_count
  obsp_count = sapply(NetS@ids, function(x) { sum(vedgesumc[x])} )
  cell_sizes = NetS@counts * Nobs
  
  ## Compute MLE Estimates
  mle_p1 = obs1_count / cell_sizes
  mle_p2 = obs2_count / cell_sizes
  mle_pc = obsc_count / cell_sizes / 2
  
  ## Compute cell-correlation sample estimates
  mle_pxy = obsp_count / cell_sizes
  num = mle_pxy - (mle_p1 * mle_p2)
  den = sqrt(mle_p1 * (1 - mle_p1) * mle_p2 * (1 - mle_p2))
  cell_corrs = num/den
  
  ## Compute cellwise log-likelihoods
  LL_null = compute_loglik_fromPC(x = obsc_count, n = 2 * cell_sizes, p = mle_pc)
  LL_alt1 = compute_loglik_fromPC(x = obs1_count, n = cell_sizes, p = mle_p1)
  LL_alt2 = compute_loglik_fromPC(x = obs2_count, n = cell_sizes, p = mle_p2)
  
  cellwise_TS = -2 * (LL_null - LL_alt1 - LL_alt2)
  
  pval_matrix = matrix(-1, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
  df_adjs = sapply(pl$cc_adj, 
                   function(x) {compute_df_adjustment2(n = cell_sizes, cell_corr = cell_corrs, cc_adj = x)})
  for(i in seq_along(pl$thres_ignore)) {
    to_keep = which(cell_sizes >= pl$thres_ignore[i])
    csq = sum(cellwise_TS[to_keep])
    dfs = apply(df_adjs[to_keep,,drop = FALSE], 2, sum)
    pval_matrix[,i] = pchisq(csq, dfs, lower.tail = FALSE)
  }
  
  return(pval_matrix)  
}



# Defines the NetworkStructHRG class

#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelRND object
#' 
#' @export
#' 
NetworkStructHRG = function(Nnodes = 10, model_param = set_model_param()) {
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelHRG(Nnodes = Nnodes, model_param = model_param)
  return(extractStruct(NetM))
}







# Generic Function Definitions --------------------------------------------

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructHRG)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructHRG = function(NetM) {
  return("tree")
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructHRG)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructHRG)
#' <What does this function do>
#' 
#' @param NetS temp
#' @param adja1 temp
#' @param adja2 temp
#' @param Nobs temp
#' @param pl temp
#' @param mode temp
#' 
#' @return temp
#' 
#' @export
#' 
computePval.NetworkStructHRG = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  if (FALSE) {
    NetM = NetworkModel(Nnodes = 30, type = "block")
    adja1 = sampleNetwork(NetM)
    adja2 = sampleNetwork(NetM)
    Nobs = 1
    pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_models = c(1,20))
    NetS = NetworkStructHRG(Nnodes = 30, model_param = set_model_param(block_nclass = 3))
    load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
  }
  
  ## Compute total count for each edge (and edgesumc is the sum of products, to be used in computing cell-wise correlations)
  if (Nobs > 1) {
    edgesum1 = apply(adja1, c(1,2), sum); edgesum2 = apply(adja2, c(1,2), sum); 
    edgesumc = apply(adja1*adja2, c(1,2), sum)
  } else if (Nobs == 1) {
    if (length(dim(adja1)) == 2) { edgesum1 = adja1 } else { edgesum1 = adja1[,,1] }
    if (length(dim(adja2)) == 2) { edgesum2 = adja2 } else { edgesum2 = adja2[,,1] }
    edgesumc = edgesum1 * edgesum2
  }
  
  ## Do counting
  obs1_count = sapply(NetS@expand, function(x) { sum(edgesum1[x[[1]], x[[2]]]) })
  obs2_count = sapply(NetS@expand, function(x) { sum(edgesum2[x[[1]], x[[2]]]) })
  obsc_count = obs1_count + obs2_count
  obsp_count = sapply(NetS@expand, function(x) { sum(edgesumc[x[[1]], x[[2]]]) })
  cell_sizes = NetS@counts * Nobs
  
  
  ## Compute MLE Estimates
  mle_p1 = obs1_count / cell_sizes
  mle_p2 = obs2_count / cell_sizes
  mle_pc = obsc_count / cell_sizes / 2
  
  ## Compute cell-correlation sample estimates
  mle_pxy = obsp_count / cell_sizes
  num = mle_pxy - (mle_p1 * mle_p2)
  den = sqrt(mle_p1 * (1 - mle_p1) * mle_p2 * (1 - mle_p2))
  cell_corrs = num/den
  
  ## Compute cellwise log-likelihoods
  LL_null = compute_loglik_fromPC(x = obsc_count, n = 2 * cell_sizes, p = mle_pc)
  LL_alt1 = compute_loglik_fromPC(x = obs1_count, n = cell_sizes, p = mle_p1)
  LL_alt2 = compute_loglik_fromPC(x = obs2_count, n = cell_sizes, p = mle_p2)
  
  cellwise_TS = -2 * (LL_null - LL_alt1 - LL_alt2)
  
  pval_matrix = matrix(-1, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
  df_adjs = sapply(pl$cc_adj, 
                   function(x) {compute_df_adjustment2(n = cell_sizes, cell_corr = cell_corrs, cc_adj = x)})
  for(i in seq_along(pl$thres_ignore)) {
    to_keep = which(cell_sizes >= pl$thres_ignore[i])
    csq = sum(cellwise_TS[to_keep])
    dfs = apply(df_adjs[to_keep,,drop = FALSE], 2, sum)
    pval_matrix[,i] = pchisq(csq, dfs, lower.tail = FALSE)
  }
  
  return(pval_matrix)  
}


# Defines the NetworkStructSBM class

#' Constructor for RND network structure
#' 
#' This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
#' 
#' @param Nnodes Number of nodes in the network model
#' @param model_param A list of model parameters -- see set_model_param()
#' 
#' @return NetworkModelSBM object
#' 
#' @export
#' 
NetworkStructSBM = function(Nnodes = 10, model_param = set_model_param()) {
  # Just generate a model and then lose the probability information. 
  NetM = NetworkModelSBM(Nnodes = Nnodes, model_param = model_param)
  return(extractStruct(NetM))
}


# Generic Function Definitions --------------------------------------------

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (getNetType.NetworkStructSBM)
#' <What does this function do>
#' 
#' @param NetM temp
#' 
#' @return temp
#' 
#' @export
#' 
getNetType.NetworkStructSBM = function(NetM) {
  return("block")
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructSBM)
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (computePval.NetworkStructSBM)
#' <What does this function do>
#' 
#' @param NetS temp
#' @param adja1 temp
#' @param adja2 temp
#' @param Nobs temp
#' @param pl temp
#' @param mode temp
#' 
#' @return temp
#' 
#' @export
#' 
computePval.NetworkStructSBM = function(NetS, adja1, adja2, Nobs, pl, mode = "default") {
  if (FALSE) {
    NetM = NetworkModel(Nnodes = 30, type = "block")
    adja1 = sampleNetwork(NetM)
    adja2 = sampleNetwork(NetM)
    Nobs = 1
    pl = list(cc_adj = c(0,1,2), thres_ignore = c(2,5,10), alphas = 0.05, n_models = c(1,20))
    NetS = NetworkStructSBM(Nnodes = 30, model_param = set_model_param(block_nclass = 3))
    load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")
  }
  
  ## Compute total count for each edge (and edgesumc is the sum of products, to be used in computing cell-wise correlations)
  if (Nobs > 1) {
    edgesum1 = apply(adja1, c(1,2), sum); edgesum2 = apply(adja2, c(1,2), sum); 
    edgesumc = apply(adja1*adja2, c(1,2), sum)
  } else if (Nobs == 1) {
    if (length(dim(adja1)) == 2) { edgesum1 = adja1 } else { edgesum1 = adja1[,,1] }
    if (length(dim(adja2)) == 2) { edgesum2 = adja2 } else { edgesum2 = adja2[,,1] }
    edgesumc = edgesum1 * edgesum2
  }
  
  ## Do counting
  COR = NetS@correct # correction factor for block models
  obs1_count = sapply(NetS@expand, function(x) { sum(edgesum1[x[[1]], x[[2]]]) }) / COR
  obs2_count = sapply(NetS@expand, function(x) { sum(edgesum2[x[[1]], x[[2]]]) }) / COR
  obsc_count = obs1_count + obs2_count
  obsp_count = sapply(NetS@expand, function(x) { sum(edgesumc[x[[1]], x[[2]]]) }) / COR
  cell_sizes = NetS@counts * Nobs
  
  ## Compute MLE Estimates
  mle_p1 = obs1_count / cell_sizes
  mle_p2 = obs2_count / cell_sizes
  mle_pc = obsc_count / cell_sizes / 2
  
  ## Compute cell-correlation sample estimates
  mle_pxy = obsp_count / cell_sizes
  num = mle_pxy - (mle_p1 * mle_p2)
  den = sqrt(mle_p1 * (1 - mle_p1) * mle_p2 * (1 - mle_p2))
  cell_corrs = num/den
  
  ## Compute cellwise log-likelihoods
  LL_null = compute_loglik_fromPC(x = obsc_count, n = 2 * cell_sizes, p = mle_pc)
  LL_alt1 = compute_loglik_fromPC(x = obs1_count, n = cell_sizes, p = mle_p1)
  LL_alt2 = compute_loglik_fromPC(x = obs2_count, n = cell_sizes, p = mle_p2)
  
  cellwise_TS = -2 * (LL_null - LL_alt1 - LL_alt2)
  
  pval_matrix = matrix(-1, nrow = length(pl$cc_adj), ncol = length(pl$thres_ignore))
  df_adjs = sapply(pl$cc_adj, 
                   function(x) {compute_df_adjustment2(n = cell_sizes, cell_corr = cell_corrs, cc_adj = x)})
  for(i in seq_along(pl$thres_ignore)) {
    to_keep = which(cell_sizes >= pl$thres_ignore[i])
    csq = sum(cellwise_TS[to_keep])
    dfs = apply(df_adjs[to_keep,,drop = FALSE], 2, sum)
    pval_matrix[,i] = pchisq(csq, dfs, lower.tail = FALSE)
  }
  
  ## Compute edgewise LL's if desired
  if (mode == "nodewise") {
    ## TODO: Make this faster if necessary. 
    ncs = rep(0, times = nrow(edgesum1))
    mlenull = 0 * edgesum1
    mlealt1 = 0 * edgesum1
    mlealt2 = 0 * edgesum1
    for(pp in seq_along(NetS@expand)) {
      mlenull[NetS@expand[[pp]][[1]], NetS@expand[[pp]][[2]]] = log(mle_pc[pp])
      mlealt1[NetS@expand[[pp]][[1]], NetS@expand[[pp]][[2]]] = log(mle_p1[pp])
      mlealt2[NetS@expand[[pp]][[1]], NetS@expand[[pp]][[2]]] = log(mle_p2[pp])
    }
    llnull = edgesum1 * mlenull + (Nobs - edgesum1) * (1 - mlenull) + edgesum2 * mlenull + (Nobs - edgesum2) * (1 - mlenull)
    llalt1 = edgesum1 * mlealt1 + (Nobs - edgesum1) * (1 - mlealt1)
    llalt2 = edgesum2 * mlealt2 + (Nobs - edgesum2) * (1 - mlealt2)
    diag(llnull) = 0; diag(alt1) = 0; diag(llalt2) = 0
    
    res = -2 * (llnull - llalt1 - llalt2)
    ncs = apply(res, 1, sum)
    
    return(list(pvals = pval_matrix, nodecontrib = ncs))
  } else if (mode == "default") {
    return(pval_matrix)  
  }
  
  
}





# Set Methods to Generics -------------------------------------------------


setMethod("sampleNetwork", signature = (NetM = "NetworkModel"), sampleNetwork.NetworkModel)
setMethod("sampleNetwork", signature = (NetM = "NetworkModelPair"), sampleNetwork.NetworkModelPair)

setMethod("getNetType", signature(NetM = "NetworkModel"), getNetType.NetworkModel)
setMethod("getNetType", signature(NetM = "NetworkModelHRG"), getNetType.NetworkModelHRG)
setMethod("getNetType", signature(NetM = "NetworkModelLSM"), getNetType.NetworkModelLSM)
setMethod("getNetType", signature(NetM = "NetworkModelPair"), getNetType.NetworkModelPair)
setMethod("getNetType", signature(NetM = "NetworkModelSBM"), getNetType.NetworkModelSBM)
setMethod("getNetType", signature(NetM = "NetworkStruct"), getNetType.NetworkStruct)
setMethod("getNetType", signature(NetM = "NetworkStructList"), getNetType.NetworkStructList)
setMethod("getNetType", signature(NetM = "NetworkStructSBM"), getNetType.NetworkStructSBM)
setMethod("getNetType", signature(NetM = "NetworkStructRND"), getNetType.NetworkStructRND)
setMethod("getNetType", signature(NetM = "NetworkStructHRG"), getNetType.NetworkStructHRG)

setMethod("getEdgeProbMat", signature = (NetM = "NetworkModel"), getEdgeProbMat.NetworkModel)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelHRG"), getEdgeProbMat.NetworkModelHRG)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelLSM"), getEdgeProbMat.NetworkModelLSM)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelPair"), getEdgeProbMat.NetworkModelPair)
setMethod("getEdgeProbMat", signature = (NetM = "NetworkModelSBM"), getEdgeProbMat.NetworkModelSBM)

setMethod("extractStruct", signature = (NetM = "NetworkModel"), extractStruct.NetworkModel)
setMethod("extractStruct", signature = (NetM = "NetworkModelHRG"), extractStruct.NetworkModelHRG)
setMethod("extractStruct", signature = (NetM = "NetworkModelPair"), extractStruct.NetworkModelPair)
setMethod("extractStruct", signature = (NetM = "NetworkModelSBM"), extractStruct.NetworkModelSBM)

setMethod("computePval", signature(NetS = "NetworkStruct"), computePval.NetworkStruct)
setMethod("computePval", signature(NetS = "NetworkStructList"), computePval.NetworkStructList)
setMethod("computePval", signature(NetS = "NetworkStructSBM"), computePval.NetworkStructSBM)
setMethod("computePval", signature(NetS = "NetworkStructRND"), computePval.NetworkStructRND)
setMethod("computePval", signature(NetS = "NetworkStructHRG"), computePval.NetworkStructHRG)

setMethod("getNnodes", signature("NetworkModel"), getNnodes.NetworkModel)
setMethod("getNnodes", signature("NetworkModelPair"), getNnodes.NetworkModelPair)
setMethod("getNnodes", signature("NetworkStruct"), getNnodes.NetworkStruct)
