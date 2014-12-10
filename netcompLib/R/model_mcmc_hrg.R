##@S Functions for doing the model MCMC (for HRG's)


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (cmn_mcmc)
#' Estimate a HRG using MCMC
#'  
#' prob_est values:
#'   default = Use the MLE estimator for each node
#'   convex = uses convex combinatino of the MLE estimate & the parent node estimate 
#'     estimate = (1-prob_est_param) * MLE + prob_est_param * Parent node
#'   convex_scaled = Weights come from number of values used to estimate the MLE
#'     estimate = wt*MLE + (1-wt)*parent node
#'     wt = num_MLE / (prob_est_param + num_MLE)
#'     where num_MLE is the number of potential edges used to estimate MLE
#' 
#' @param obsadjs Input observed network models
#' @param iterations Number of MCMC iterations
#' @param tree Starting HRG tree
#' @param verbose T/F: verbose output?
#' @param only.the.max T/F: Only return best tree?
#' @param prob_est Probability estimation method: "default", "convex", or "convex_scaled". 
#' @param prob_smooth_param Parameter for probability smoothing (if appropriate) [0 means no smoothing]
#' @param smooth_root T/F: Smooth the root? (depends on probability estimation method)
#' 
#' @return Best-fitting HRG model (or a number of these)
#' 
#' @export
#' 
cmn_mcmc = function(obsadjs, iterations=1000, tree=NULL, verbose=TRUE, only.the.max=TRUE, prob_est = c("default", "convex", "convex_scaled"), prob_smooth_param = 0.5, smooth_root = TRUE) {
  # Originally written by Andrew Thomas, but modified
  #Presently, the only moves we allow are exchanges -- no mergers or splits.
  
  WT = prob_smooth_param
  
  ## smooth_root : whether to apply smoothing to root of tree (smoothed
  ##   towards overall probability of edge
  
  #nodes=100; obsadjs=array(rbinom(nodes*nodes,1,0.3), rep(nodes,2)); diag(obsadjs)=0; obsadjs=1*(obsadjs+t(obsadjs)>0); tree=NULL; verbose=TRUE; iterations=10000
  
  #obsadjs=netdraw;  iterations=100; tree=NULL; verbose=TRUE; only.the.max=TRUE
  
  #prep data.
  if (is.na(dim(obsadjs)[3])) obsadjs <- array(c(obsadjs), c(dim(obsadjs)[1:2], 1))
  
  nn <- dim(obsadjs)[1]
  nets <- dim(obsadjs)[3]
  s.one <- apply(obsadjs, c(1,2), sum)
  
  if (is.null(tree)) tree <- starter_tree(nn)
  ancestor.info <- closest_ancestor(tree)
  
  tree <- update_probabilities (tree, ancestor.info$anc.table,
                                s.one, nets)
  int.nodes <- length(tree$prob)
  
  #log-likelihood pieces.
  log.entropy <- function(pp) ifelse(pp==0, 0, ifelse(pp==1, 0, pp*log(pp)+(1-pp)*log(1-pp)))  #because we're fixing pp = yy/nn, it's implicit here.
  
  loglik <- function(which.internal.nodes, this.anc.table, this.tree)
    #which.internal.nodes=c(73, 61); this.anc.table=ancestor.info$anc.table; this.tree=tree
    sapply(which.internal.nodes, function(rr) sum(this.anc.table==rr))*nets*  #this is the binomial "n".
    log.entropy(this.tree$prob[which.internal.nodes-nn])
  
  pp.part <- function(node.id, anc.table) {
    return(mean(s.one[anc.table==node.id])/nets)
  }  
  
  entry.ge.k <- function(vec, k) {
    return(vec[vec > k])
  }
  
  tree.ret <- as.list(rep(numeric(0), iterations))
  loglik.ret <- rep(NA, iterations)
  
  tree.prop <- tree
  
  for (kk in 1:iterations) {

    #pick an internal node that isn't on top.
    cur <- nn + sample(2:int.nodes, 1)
    cur.parent <- nn + which(sapply(tree$children, function(chch) sum(chch == cur))>0)
    
    #pick a sibling of "cur".
    sibling <- tree$children[[cur.parent-nn]]; sibling <- sibling[sibling != cur]
    if (length(sibling)>1) {
      sibling.pos <- sample(length(sibling),1);
      sibling <- sibling[sibling.pos]
    } else {
      sibling.pos <- which(tree$children[[cur.parent-nn]]==sibling)
    }
    
    #pick a child of "cur".
    child.pos <- sample(length(tree$children[[cur-nn]]), 1)
    child <- tree$children[[cur-nn]][child.pos]
    
    #c(cur, cur.parent, sibling, child)
    
    #make the switch.
    tree.prop$children[[cur.parent-nn]][sibling.pos] <- child
    tree.prop$children[[cur-nn]][child.pos] <- sibling
    
    tree.prop$parents[sibling] <- cur
    tree.prop$parents[child] <- cur.parent
        
    ancestor.info.prop <- closest_ancestor(tree.prop)
    
    if (prob_est == "default") {
      tree.prop$prob[cur-nn] <- pp.part (cur, ancestor.info.prop$anc.table)
      tree.prop$prob[cur.parent-nn] <- pp.part (cur.parent, ancestor.info.prop$anc.table)
    } else {
      overallprob = (sum(s.one)/(nn*(nn-1)))/nets
      if (prob_est == "convex") {
        MLE = pp.part(nn+1, ancestor.info.prop$anc.table)
        if (smooth_root) {
          tree.prop$prob[1] <- (1-WT)*MLE + (WT)*overallprob
        } else {
          tree.prop$prob[1] <- MLE
        }
        to_estimate = entry.ge.k(tree.prop$children[[1]], nn)
        
        while(length(to_estimate) > 0) {
          CUR = to_estimate[1]
          to_estimate = to_estimate[-1]
          
          MLE = pp.part(CUR, ancestor.info.prop$anc.table)
          tree.prop$prob[CUR-nn] =
            (1-WT)*MLE + (WT)*tree.prop$prob[tree.prop$parents[CUR]-nn]
          
          to_estimate = c(to_estimate,
                          entry.ge.k(tree.prop$children[[CUR-nn]], nn) )
        }    
      } else if (prob_est == "convex_scaled") {
        MLE = pp.part(nn+1, ancestor.info.prop$anc.table)
        NUM_MLE = sum(ancestor.info.prop$anc.table == (nn+1))*nets/2
        ADJ_WT = NUM_MLE / (WT + NUM_MLE)
        
        if (smooth_root) {
          tree.prop$prob[1] <- ADJ_WT*MLE + (1-ADJ_WT)*overallprob
        } else {
          tree.prop$prob[1] <- MLE
        }
        to_estimate = entry.ge.k(tree.prop$children[[1]], nn)
        
        while(length(to_estimate) > 0) {
          CUR = to_estimate[1]
          to_estimate = to_estimate[-1]
          
          MLE = pp.part(CUR, ancestor.info.prop$anc.table)
          NUM_MLE = sum(ancestor.info.prop$anc.table == CUR)*nets/2
          ADJ_WT = NUM_MLE / (WT + NUM_MLE)
          tree.prop$prob[CUR-nn] <-
            ADJ_WT*MLE + (1-ADJ_WT)*tree.prop$prob[tree.prop$parents[CUR]-nn]
          
          to_estimate = c(to_estimate,
                          entry.ge.k(tree.prop$children[[CUR-nn]], nn) )
        }
        
        
      }
    }
    
    #compare likelihoods.
    loglike.orig <- sum(loglik(c(cur, cur.parent), ancestor.info$anc.table, tree))
    loglike.prop <- sum(loglik(c(cur, cur.parent), ancestor.info.prop$anc.table, tree.prop))
    
    pick.proposal = loglike.prop-loglike.orig > -rexp(1)
    if(is.na(pick.proposal)) {
      print(tree)
      break
    }
    if (pick.proposal) {
      tree <- tree.prop; ancestor.info <- ancestor.info.prop
    } else {
      tree.prop <- tree; ancestor.info.prop <- ancestor.info
    }
    
    tree.ret[[kk]] <- tree
    loglik.ret[kk] <- sum(loglik(nn+1:int.nodes, ancestor.info$anc.table, tree))
    
    if (verbose & (kk %% (iterations/10) == 0)) print(kk)
  }
  
  ml.achieved <- min(which(loglik.ret==max(loglik.ret)))
  
  out.obj <- list(tree.ret=tree.ret,
                  loglik.ret=loglik.ret,
                  ancestor.info.final=ancestor.info,
                  ml.achieved=ml.achieved,
                  choice.tree.ret=tree.ret[[ml.achieved]],
                  iterations=iterations)
  
  if (only.the.max) {
    out.obj$tree.ret <- NULL
  }
  return(out.obj)
}




## TODO: [Figure out] what this function does... 
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (update_probabilities)
#' <<BasicInfo>> 
#' 
#' @param tree temp
#' @param anc.table temp
#' @param net.compressed temp
#' @param max.tie.single temp
#' @param subset temp
#' 
#' @return temp
#' 
#' @export
#' 
update_probabilities = function(tree, anc.table,
                                net.compressed, max.tie.single,
                                subset=1:length(tree$children)) {
  
  # Written by Andrew Thomas
  
  #anc.table=anc; net.compressed=s.one; max.tie.single=nets; subset=1:length(tree$children)
  
  
  ## This is the original code. I've commented it out... 
  #  tree$prob[subset] <- sapply(subset, function(kk) mean(net.compressed[anc.table==kk+tree$nodes])/max.tie.single)
  N = dim(net.compressed)[1]
  tree$prob[subset] <- 0.5
  ## TODO: Input appropriate estimate here? 
  return(tree)
  
}

## Function to call the C version of MCMC
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (run_C_mcmc)
#' <<BasicInfo>> 
#' 
#' @param par_tree temp
#' @param adj_mat temp
#' @param num_iters temp
#' @param convex_wt temp
#' @param verbose temp
#' @param rngseed temp
#' @param simplify temp
#' 
#' @return temp
#' 
#' @export
#' 
run_C_mcmc = function(
  par_tree, adj_mat,
  num_iters = 1000, convex_wt = 0,
  verbose = FALSE, rngseed = NULL,
  simplify = TRUE) {
  
  num_obs = 1
  if (length(dim(adj_mat)) == 3) {
    num_obs = dim(adj_mat)[3]
    adj_mat = apply(adj_mat, 1:2, sum)
  }
  
  plen = length(par_tree)
  
  if (is.null(rngseed)) {
    rngseed = sample(1:300000, size = 1)
  }
  
  res = .C("cmn_mcmc",
           #|----##Updated this function (see model_mcmc.R), so need to update usages as necessary --Tue Sep 16 16:35:10 2014--
           par = as.integer(par_tree), pvlen = as.integer(plen),
           edgemat = as.integer(as.vector(adj_mat)), num_obs = as.integer(num_obs),
           mcmc_iter = as.integer(num_iters),
           wt = as.double(convex_wt),
           o_lik = double(num_iters), o_sib = double(num_iters), o_chi = double(num_iters),
           o_bestTree = integer(plen), o_bestProbs = double(plen), o_bestIter = integer(1),
           verbose = as.integer(verbose),
           rngseed = rngseed)
  
  if (simplify) { ## eventually, dont return everything. for now , ignore
  } else {
  }
  
  res$bestTree = list()
  res$bestTree$nodes = (plen+1)/2
  res$bestTree$prob = res$o_bestProbs[(res$bestTree$nodes+1):plen]
  res$bestTree$parents = res$o_bestTree +1
  
  res$bestTree$children = try( tree_from_parents(res$bestTree$parents) )
  if(class(res$bestTree$children) == "try-error") {
    print(res$bestTree$parents)
    print(paste("Best iteration: ",res$o_bestIter))
    stop("Caught error in vbestTree$parents ^^")
  }
  
  return(res)
}


