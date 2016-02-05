##@S This file contains the main network comparison function. 


#' Compares two networks
#' 
#' This function compares a pair of networks [add stuff on multiple pairs of networks if implemented properly]
#' 
#' @return p-value of some sort? FIX THIS
#' 
#' @export
#' 
compareNetworks = function() {
  ## Also -- allow for random partitions? 
  
  ## Step 0 -- Clean / check input data
  
  ## Step 1 -- Hide part of the input adjacency arrays
  
  ## Step 2 -- Fit appropriate dyad parition to the non-hidden part of the data
  
#   ## method = 'spectral' or 'mf'
#   ## method 'mf' for mean field EM approach, 'spectral' for just doing spectral clustering
#   if (cl$SBM_method == "spectral") {
#     ## For spectral clustering, simply apply spectral clustering 
#     return(specClust(adjm, cl$SBM_Nclass, cl$Ntries))
#     
#   } else 
  
  ## Step 3 -- Compute test statistic on the fitted model
  
  ## Step 4 -- Apply multiple comparisons methods if needed for doing multiple data splits
  
  
}
