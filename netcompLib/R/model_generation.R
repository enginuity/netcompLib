##@S General model generation functions

#' Generates a pair of network observations from a pair of models (they can be the same model ... under null)
#' 
#' @param gen_model gen_model object containing both models
#' 
#' @return Array of two adjacency matrices
#' 
#' @export
#' 
gen_network_pair = function(gen_model) {
  require(abind)
  
  if (gen_model$mode == "tree") {
    adjmat = abind(cmn_network(gen_model$m1, 1), cmn_network(gen_model$m2, 1))
  } else if (gen_model$mode == "block" | gen_model$mode == "blockmodel") {
    adjmat = abind(block_network(gen_model$m1, 1), block_network(gen_model$m2, 1))
  } else if (gen_model$mode == "latent") {
    adjmat = abind(latent_network(gen_model$m1, 1), latent_network(gen_model$m2, 1))
  }
  return(adjmat)
}


#' Generate model parameters in several classes
#' 
#' @param NN Number of network nodes
#' @param mode Type of model
#' @param K Parameter for block model if applicable, otherwise, dimension of space for latent space model
#' @param is_null If TRUE: This assumes the null hypothesis is true (makes m2 = m1)
#' 
#' @return Model information in list format (gen_model)
#' 
#' @export
#' 
gen_model_fx = function(NN, mode, K = 3, is_null = FALSE) {
  if (mode == "tree") {
    m1 = gen_tree(NN, random_type = "random", random_plimit = c(0.1, 0.9))
    m2 = gen_tree(NN, random_type = "random", random_plimit = c(0.1, 0.9))
  } else if (mode == "block" | mode == "blockmodel") {
    m1 = gen_block(NN, K = K, pmin = 0.1, pmax = 0.9)
    m2 = gen_block(NN, K = K, pmin = 0.1, pmax = 0.9)
  } else if (mode == "latent") {
    m1 = gen_latent(NN, K = K, D = 3, sd_center = 5, gen_norm = FALSE)
    m2 = gen_latent(NN, K = K, D = 3, sd_center = 5, gen_norm = FALSE)
  }
  
  if (is_null) { m2 = m1 }
  return(list(mode = mode, m1 = m1, m2 = m2, is_null = is_null, Nnodes = NN))
}


#' Generate a number of fixed model structures
#' 
#' @param mode Either 'tree', 'blockmodel', or 'random'
#' @param Nnodes Number of nodes in the network
#' @param Nmodels Number of models to generate
#' @param bm_Nclasses Number of classes in block model
#' @param rnd_Ngroups Number of random edge groups in 'random' mode
#' 
#' @return A list of form fit_model of fixed models: If this is a 'tree', return a list of fixed trees and expanded children. If this is 'blockmodel', return groups, counts, expanded. If random, return group sizes & index vector. 
#' 
#' @export
#' 
generate_fitting_models = function(mode = c("tree", "block", "random"), Nnodes, Nmodels,
                                   bm_Nclasses = 3, rnd_Ngroups = 10) {
  res_list = list()
  for(j in 1:Nmodels) {
    if (mode == "tree") {
      tr = gen_tree(Nnodes = Nnodes, random_type = "random")
      expc = expanded_children_from_tree(tr)
      
      res_list[[j]] = list(
        fixed_tree = tr, ft_expand = expc,
        ft_counts = sapply(expc, function(x) {length(x[[1]]) * length(x[[2]]) }) )
      
    } else if (mode == "blockmodel" | mode == "block") {
      ga = sample(1:bm_Nclasses, size = Nnodes, replace = TRUE)
      
      counts = rep(0, times = bm_Nclasses^2)
      expanded = list()
      cur = 1
      for(k in 1:bm_Nclasses) {
        for(m in 1:bm_Nclasses) {
          if (k == m) {
            counts[cur] = sum(ga == k) * (sum(ga == k) - 1) / 2
          } else if (k != m) {
            counts[cur] = sum(ga == k) * sum(ga == m)
          }
          expanded[[cur]] = list(which(ga == k), which(ga == m))
          cur = cur + 1
        }}
      res_list[[j]] = list(bm_groups = ga, bm_counts = counts, bm_expand = expanded)
      
    } else if (mode == "random") {
      idm = matrix(1:(Nnodes^2), nrow = Nnodes)
      idm[lower.tri(x = idm, diag = TRUE)] = 0
      idm = as.vector(idm)
      good_ids = idm[idm != 0]
      
      rand_order = sample(good_ids, size = length(good_ids), replace = FALSE)
      sizes = rep(floor(length(good_ids) / rnd_Ngroups), times = rnd_Ngroups)
      m = length(good_ids) %% rnd_Ngroups
      if (m > 0) {
        sizes[1:m] = sizes[1:m] + 1
      }
      
      counted = 0
      temp = list()
      for(k in 1:rnd_Ngroups) {
        cur = 1 + counted
        temp[[k]] = sort(rand_order
                         [cur:(cur+sizes[k])])
        counted = counted + sizes[k]
      }
      
      res_list[[j]] = list(rs_counts = sizes, rs_ids = temp)
    }
  }
  return(list(Nnodes = Nnodes, Nmodels = Nmodels, mode = mode, model_list=res_list))
}


# Newer versions of code --------------------------------------------------


#' Create a list of model parameters
#' 
#' This function sets up the model parameters to be passed into model generation functions (eg. sets the max number of blocks in a SBM). There are default parameters that are used if this function is called with no arguments. 
#' 
#' @param pmin Minimal possible edge probability
#' @param pmax Maximal possible edge probability
#' @param block_nclass Number of blocks in block model
#' @param block_avgdensity Set average density in block model (ignored if NULL)
#' @param random_ngroups Number of groups for completely random edge partition
#' @param tree_type Randomization on structure of tree: can be "left", "random", or "original"
#' @param latent_dim Dimension of latent space in latent space models
#' @param latent_nclass Number of clusters in latent space model
#' @param latent_sdcenter SD on centers of latent space model
#' @param latent_isgennorm If TRUE: Uses normal distribution for latent locations. Otherwise, uses uniform distribution. 
#' 
#' @return Return list of parameters
#' 
#' @export
#' 
set_model_param = function(pmin = 0, pmax = 1, block_nclass = 3, block_avgdensity = NULL, random_ngroups = 10, tree_type = "random", latent_dim = 3, latent_nclass = 3, latent_sdcenter = 5, latent_isgennorm = TRUE) {
  return(list(pmin = pmin, pmax = pmax, block_nclass = block_nclass, block_avgdensity = block_avgdensity, random_ngroups = random_ngroups, tree_type = tree_type, latent_dim = latent_dim, latent_nclass = latent_nclass, latent_sdcenter = latent_sdcenter, latent_isgennorm = latent_isgennorm))
}


## TODO: Create a function that only samples a single network model?

#' Generates a pair of random network models
#' 
#' Generates a pair of network models of a specific type. If the null is true, then only one network model is generated. 
#' 
#' @param Nnodes Number of network nodes
#' @param mode Type of model
#' @param model_param List of model parameters (see set_model_param)
#' @param is_null Generate two models under null hypothesis? (TRUE means m2 = m1.)
#' 
#' @return Model information in list format (gen_model)
#' 
#' @export
#' 
sample_generating_models = function(Nnodes, mode, model_param = set_model_param(), is_null = FALSE) {
  ## TODO: [Obseletes] gen_model_fx. 
  
  if (mode == "tree") { samp_fx = sample_model_tree }
  if (mode == "block") { samp_fx = sample_model_block }
  if (mode == "latent") { samp_fx = sample_model_latent }
  
  m1 = samp_fx(Nnodes = Nnodes, model_param = model_param)
  if (is_null) { m2 = m1 } else { m2 = samp_fx(Nnodes = Nnodes, model_param = model_param) }
  
  return(list(mode = mode, m1 = m1, m2 = m2, is_null = is_null, Nnodes = Nnodes))
}


#' Samples a pair of network observations
#' 
#' From a pair of network models, this function samples network observations from each model. The number of networks sampled can vary. 
#' 
#' @param gen_model gen_model object containing both models
#' @param Nobs Number of network observation pairs to generate
#' 
#' @return Array of two adjacency matrices
#' 
#' @export
#' 
sample_network_pair = function(gen_model, Nobs = 1) {
  require(abind)
  ## TODO: [Obseletes] gen_network_pair
  
  if (gen_model$mode == "tree") {
    adjmat = list(cmn_network(gen_model$m1, Nobs), cmn_network(gen_model$m2, Nobs))
  } else if (gen_model$mode == "block" | gen_model$mode == "blockmodel") {
    adjmat = list(block_network(gen_model$m1, Nobs), block_network(gen_model$m2, Nobs))
  } else if (gen_model$mode == "latent") {
    adjmat = list(latent_network(gen_model$m1, Nobs), latent_network(gen_model$m2, Nobs))
  }
  
  return(adjmat)
}


#' Generates random network model structures
#' 
#' This function generates network model structures. Essentially, these store information only about the edge partitions, but with nothing given as the edge probabilities. Thus, it wouldn't be possible to sample networks from these models, since they lack probabilities. 
#' 
#' It's also possible to convert a more specific model (generating models) into this format, by passing it in. Note that latent network models cannot work in this way (since they usually do not correspond to a discrete edge partitioning).
#' 
#' @param mode Either 'tree', 'block', or 'random'
#' @param Nnodes Number of nodes in the network
#' @param Nmodels Number of models to generate
#' @param model_param List of model parameters
#' @param gen_model If NON-null: Generates fitting models for the format of input given gen_model (m1)
#' 
#' @return A list of form fit_model of fixed models: If this is a 'tree', return a list of fixed trees and expanded children. If this is 'blockmodel', return groups, counts, expanded. If random, return group sizes & index vector. 
#' 
#' @export
#' 
sample_fitting_models = function(mode, Nnodes, Nmodels, model_param = set_model_param(), gen_model = NULL) {
  ## if gen_model isn't null: IGNORE everything else, convert gen_model$m1 into fitting model. 
  if (!is.null(gen_model)) { 
    Nmodels = 1; mode = gen_model$mode
    print("Ignoring Nmodels: generating fitting form for given model") 
  }
  
  ## TODO: [Obselete old function] This function OBSELETES generate_fitting_models. So need to remove that function?
  res_list = list()
  for(j in 1:Nmodels) {
    if (mode == "tree") {
      tr = gen_tree(Nnodes = Nnodes, random_type = "random")
      
      ## Replace with m1 as necessary
      if (!is.null(gen_model)) { tr = gen_model$m1 }
      
      expc = expanded_children_from_tree(tr)
      
      res_list[[j]] = list(
        fixed_tree = tr, ft_expand = expc,
        ft_counts = sapply(expc, function(x) {length(x[[1]]) * length(x[[2]]) }) )
      
    } else if (mode == "blockmodel" | mode == "block") {
      NClass = model_param$block_nclass
      if (Nnodes <= NClass) {stop("Too few sbm classes for # of nodes. ")}
      ga = c(1:NClass, sample(1:NClass, size = (Nnodes - NClass), replace = TRUE))
      ga = sample(ga, size = length(ga), replace = FALSE)
      
      ## Replace with m1 as necessary
      if (!is.null(gen_model)) { ga = gen_model$m1$assign; NClass = length(unique(ga)) }
      
      counts = rep(0, times = NClass + NClass * (NClass - 1) / 2)
      correction = rep(0, times = NClass + NClass * (NClass - 1) / 2)
      expanded = list()
      cur = 1
      for(k in 1:NClass) {
        for(m in 1:NClass) {
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
      res_list[[j]] = list(bm_groups = ga, bm_counts = counts, bm_expand = expanded, bm_correct = correction)
      
    } else if (mode == "random") {
      NGroups = model_param$random_ngroups
      
      idm = matrix(1:(Nnodes^2), nrow = Nnodes)
      idm[lower.tri(x = idm, diag = TRUE)] = 0
      idm = as.vector(idm)
      good_ids = idm[idm != 0]
      if (length(good_ids) < NGroups) { stop("Too many groups for too few edges") }
      
      rand_order = sample(good_ids, size = length(good_ids), replace = FALSE)
      sizes = rep(floor(length(good_ids) / NGroups), times = NGroups)
      m = length(good_ids) %% NGroups
      if (m > 0) {
        sizes[1:m] = sizes[1:m] + 1
      }
      
      counted = 0
      temp = list()
      for(k in 1:NGroups) {
        temp[[k]] = sort(rand_order[seq(from = counted+1, length.out = sizes[k])])
        counted = counted + sizes[k]
      }
      
      res_list[[j]] = list(rs_counts = sizes, rs_ids = temp)
    }
  }
  return(list(Nnodes = Nnodes, Nmodels = Nmodels, mode = mode, model_list=res_list))
}


