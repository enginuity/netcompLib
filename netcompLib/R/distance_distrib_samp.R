##@S This file contains functions used to simulate distributions of distance functions

## General functions

## TODO: [Cleanup] Remove all functions in this file


## Simulation functions
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (sim_dist_true_model)
#' <<BasicInfo>> 
#' 
#' @param base_tree temp
#' @param B temp
#' @param obs temp
#' 
#' @return temp
#' 
#' @export
#' 
sim_dist_true_model = function(base_tree, B = 1000, obs = 1) {
  ## This doesnt need MCMC: assume we know the true tree, but dont know the probs. then see how good the ests are on average, using the measure of distance.

  probs = edge_probs(base_tree)
  probs_count = table(probs)[-1]

  probs = as.numeric(names(probs_count))

  res = rep(-1, times = B)
  for(k in 1:B) {
    accum_dist = 0
    for(j in 1:length(probs)) {
      samp1 = rbinom(n = probs_count[j]*obs, size = 1, prob = probs[j])
      samp2 = rbinom(n = probs_count[j]*obs, size = 1, prob = probs[j])
      accum_dist = accum_dist + probs_count[j] * abs(mean(samp1) - mean(samp2))
    }
    res[k] = accum_dist
  }

  return(res)
}














## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (compute_dist_distrib)
#' <<BasicInfo>> 
#' 
#' @param base_tree temp
#' @param num_sims temp
#' @param num_mcmc_iters temp
#' @param obs temp
#' @param convex_wt temp
#' @param m_tree_dist temp
#' @param dist_true temp
#' 
#' @return temp
#' 
#' @export
#' 
compute_dist_distrib = function(base_tree, num_sims = 100, num_mcmc_iters = 100000, obs = 1, convex_wt = 5, m_tree_dist, dist_true = FALSE) {
  ## For each i in 1:num_sims, this draws two graphs from base_tree and then estimates trees from them, then computes tree-distance.

  tenth = floor(num_sims/10)
  if(!(tenth > 2)) {tenth = 2}
  
  NN = base_tree$nodes
  res = array(dim = c(NN,NN,num_sims))
  
  ## REPEAT THIS:
  for(j in 1:num_sims) {
    if(j %% tenth == 0) {cat(".")}
    adj1 = cmn_network(base_tree, obs)
    adj2 = cmn_network(base_tree, obs)
    rt1 = gen_tree(base_tree$nodes, type = "random")
#|----##Check usage; this function updated to having 'random_plimit' as a parameter 9/15/2014 --Mon Sep 15 02:55:46 2014--
    rt2 = gen_tree(base_tree$nodes, type = "random")
#|----##Check usage; this function updated to having 'random_plimit' as a parameter 9/15/2014 --Mon Sep 15 02:55:46 2014--
    
    est1 = run_C_mcmc(
      par_tree = rt1$parents,
      adj_mat = adj1, num_iters = num_mcmc_iters,
      convex_wt = convex_wt)
    
    est2 = run_C_mcmc(
      par_tree = rt2$parents,
      adj_mat = adj2, num_iters = num_mcmc_iters,
      convex_wt = convex_wt)

    if (dist_true) {
      res[,,j] = m_tree_dist(est1$bestTree, est2$bestTree, base_tree)
    } else {
      res[,,j] = m_tree_dist(est1$bestTree, est2$bestTree)
    }
  }
  return(res)
  ##return(list(distrib = apply(res, 3, sum),
  ##           cellwise = apply(res,c(1,2),mean) ))
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (plot_tree)
#' <<BasicInfo>> 
#' 
#' @param tree temp
#' @param mean temp
#' @param cell_mean temp
#' @param cell_se temp
#' 
#' @return temp
#' 
#' @export
#' 
plot_tree = function(tree, mean, cell_mean, cell_se)  {
  NN = tree$nodes
  gNN = 2*NN - 1
  ftmat = matrix(0, ncol =2, nrow = (NN - 1) *2 )
  for(j in 1:(NN-1)) {
    ftmat[2*j - 1,] = c(NN+j, tree$children[[j]][1])
    ftmat[2*j,] =c(NN+j, tree$children[[j]][2])
  }
  g = ftM2graphNEL(ftmat)

  temp = expanded_children_from_tree(tree)
  child_nums = sapply(temp, function(x) {length(x[[1]]) * length(x[[2]])})

  contribs = rep(NA, times = NN-1)
  for(j in 1:(NN-1)) {
    cur = temp[[j]]
    ##zscores = (cell_mean[cur[[1]],cur[[2]]] - mean) / cell_se[cur[[1]],cur[[2]]]
    zscores = cell_mean[cur[[1]], cur[[2]]]
    contribs[j] = mean(zscores)
  }

  col_vec = c(rep("black", times = NN), rep("#FFDDAA", times = NN-1))
  label_vec = c(1:NN,
    paste((NN+1):gNN," : N=", child_nums,
          "\\\np=", round(tree$prob, 3),
          "\\\nc=", round(contribs, 2),
          sep = ""))
  height_vec = c(rep(0.4, times = NN), rep(1.2, times = NN-1))
  width_vec = c(rep(0.4, times = NN), rep(1.6, times = NN-1))
  shape_vec = c(rep("circle",times = NN), rep("rect", times = NN-1))
  ##fontsize_vec = c(rep(30, times = NN), rep(55, times = NN-1))
  ##names(fontsize_vec) = 1:gNN
  names(col_vec) = 1:gNN
  names(label_vec) = 1:gNN
  names(width_vec) = 1:gNN
  names(height_vec) = 1:gNN
  names(shape_vec) = 1:gNN
  nodeAttr = list(color = col_vec,
    label = label_vec,
    width = width_vec,
    height = height_vec,
    shape = shape_vec
    #fontsize = fontsize_vec
    )
  edge_width_vec = rep(8, times = 2*NN - 2)
  names(edge_width_vec) = paste(ftmat[,1], "~", ftmat[,2], sep = "")
  edgeAttr = list(penwidth=edge_width_vec)
  plot(g, nodeAttrs = nodeAttr, edgeAttrs = edgeAttr)
  return(NULL)
}

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (proc_res)
#' <<BasicInfo>> 
#' 
#' @param test_tree temp
#' @param res temp
#' @param plots temp
#' @param xmax temp
#' @param main temp
#' 
#' @return temp
#' 
#' @export
#' 
proc_res = function(test_tree,res, plots = 2, xmax = 40, main = "Plot") {
  n = dim(res)[1]
  N = dim(res)[3] 
  distrib = apply(res, 3, sum)
  cellwise_mean = apply(res, c(1,2), mean)
  overallmean = mean(distrib)/(n * (n-1))
  cellwise_se = apply(res, c(1,2), sd)/sqrt(N)

  if(plots == 2) {
    par(mfrow = c(1,2))
    plot_tree(test_tree, mean = overallmean, cell_mean = cellwise_mean,
              cell_se = cellwise_se)
    par(mar =  c(5, 4, 4, 2) + 0.1)
  }
  hist(distrib, xlim = c(0, xmax), freq = FALSE, main = main, xlab = "Distance")
  points(density(distrib), type = "l", col = "darkgreen", lwd = 2)
  
  
  return(NULL)
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (obtain_contribs)
#' <<BasicInfo>> 
#' 
#' @param NN temp
#' @param num_trees temp
#' @param distfx temp
#' @param avg temp
#' @param dist_true temp
#' 
#' @return temp
#' 
#' @export
#' 
obtain_contribs = function(NN, num_trees = 10, distfx, avg = TRUE, dist_true = FALSE) {
  ## dist_true => true if distance function takes as 3rd tree the true tree.
  ## this should never happen in practice, but is used to test theoretical...
  if (!avg) {
    xy_obtain = matrix(0, ncol = 2, nrow = (NN * (NN-1)/2))
    count = 1
    for(k in 1:(NN-1)) {
      for(l in (k+1):NN){
        xy_obtain[count,] = c(k,l)
        count = count + 1
      }}
    
    res_tab = matrix(-1, nrow = num_trees * 3 * (NN) * (NN-1)/2, ncol = 3)
    for(B in 1:num_trees) {
      cat("***",B,"***\n")
      tree = gen_tree(NN,"random")
#|----##Check usage; this function updated to having 'random_plimit' as a parameter 9/15/2014 --Mon Sep 15 02:55:46 2014--
      res = compute_dist_distrib(tree, num_sims = 3, num_mcmc_iters = 200, m_tree_dist = distfx, dist_true = dist_true)

      temp = expanded_children_from_tree(tree)
      child_nums = sapply(temp, function(x) {length(x[[1]]) * length(x[[2]])})
      ntree = tree
      ntree$prob = child_nums
      m_prob = edge_probs(tree)
      m_cnum = edge_probs(ntree)

      for(k in 1:3) {
        num_vals = NN * (NN-1)/2
        shift = ((B-1) * 3 + (k-1))* num_vals
        for(m in 1:num_vals) {
        
          res_tab[shift + m,1] = m_cnum[xy_obtain[m,1], xy_obtain[m,2]]
          res_tab[shift + m,2] = m_prob[xy_obtain[m,1], xy_obtain[m,2]]
          res_tab[shift + m,3] = res[xy_obtain[m,1], xy_obtain[m,2],k]
        }
      }
    }
  } else {
    res_tab = matrix(-1, nrow = num_trees * (NN-1), ncol = 3)
    
    for(B in 1:num_trees) {
      cat("***",B,"***\n")
      tree = gen_tree(NN,"random")
#|----##Check usage; this function updated to having 'random_plimit' as a parameter 9/15/2014 --Mon Sep 15 02:55:46 2014--
      res = compute_dist_distrib(tree, num_sims = 100, num_mcmc_iters = 200, m_tree_dist = distfx, dist_true = dist_true)
      
      n = dim(res)[1]
      N = dim(res)[3] 
      distrib = apply(res, 3, sum)
      cellwise_mean = apply(res, c(1,2), mean)
      overallmean = mean(distrib)/(n * (n-1))
      cellwise_se = apply(res, c(1,2), sd)/sqrt(N)
      
      NN = tree$nodes
      gNN = 2*NN - 1
      
      temp = expanded_children_from_tree(tree)
      child_nums = sapply(temp, function(x) {length(x[[1]]) * length(x[[2]])})
      contribs = rep(NA, times = NN-1)
      for(j in 1:(NN-1)) {
        cur = temp[[j]]
        zscores = cellwise_mean[cur[[1]], cur[[2]]]
        contribs[j] = mean(zscores)
      }
      res_tab[1:(NN-1) + (NN-1)*(B-1),] = cbind(child_nums, tree$prob, contribs)
    }
  }
  
  colnames(res_tab) = c("Num_edges", "Prob", "Contrib")
  return(res_tab)
}


