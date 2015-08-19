

#' Plots a heatmap of a probability matrix
#' 
#' Generates a heatmap of a probability matrix. This allows the probabilities to also be a difference between two probability matrices, and allows custom input of color scheme. 
#' 
#' @param pm [matrix-double] :: Probability matrix
#' @param add_legend [logical; DEFAULT = T] :: If TRUE, adds a legend plot
#' @param is_dm [logical; DEFAULT = F] :: If TRUE, input matrix is on a [-1,1] scale (as a difference of probability matrices)
#' @param only_legend [logical; DEFAULT = F] :: If TRUE: ONLY plot the legend
#' @param colors [vector-char] :: Colors to use in the heat map (as a vector of hex values)
#' @param ... [] :: Other variables to be passed to 'plot'
#' 
#' @return [nothing] :: Produces a plot (standard R plot functions)
#' 
#' @export
#' 
plot_probmatrix = function(pm, add_legend = TRUE, is_dm = FALSE, only_legend = FALSE, 
                           colors = rainbow(n = 200, s = 1, v = 1, start = 0, end = 1/3, alpha = 1), ...) {
  
  ncol = length(colors)
  diag(pm) = NA
  if (add_legend & !only_legend) par(mfrow = c(1,2))
  
  zlim = c(0,1)
  if (is_dm) zlim = c(-1,1) 
  image(x = 1:dim(pm)[1], y = 1:dim(pm)[1], z = pm, zlim = zlim, xlab = "Node Number", ylab = "Node Number", col = colors, ...)
  
  if (add_legend | only_legend) {
    ys = seq(ifelse(is_dm, -1, 0), 1, length.out = ncol)
    plot(x = rep(1, times = length(ys)), y = ys, pch = 15, col = colors, xlab = "", ylab = "Edge Probability", main = "Legend")
    for(j in 1:50) {
      points(x = rep(1 + j/100, times = length(ys)), y = ys, pch = 15, col = colors)
      points(x = rep(1 - j/100, times = length(ys)), y = ys, pch = 15, col = colors)
    }
  }
  invisible(0)
}


#' Plot the Edge Group structure of a network model
#' 
#' @param NetM [\code{\link{NetworkModel}}] :: Model to plot edge group of
#' @param ... [] :: Further parameters to pass on to the plot function
#' 
#' @return [int] :: Ignore output; a plot is produced
#' 
#' @export
#' 
plot_edgegroup = function(NetM, ...) {
  z = getEdgeProbMat(NetM)
  u = sort(unique(as.vector(z)), decreasing = FALSE)
  for(j in 1:nrow(z)) {
    z[j,] = match(z[j,], u)
  }
  
  mat = z - 1
  
  N = nrow(mat)
  plot(-5, -5, xlim = c(0,N) + 0.5, ylim = c(0,N) + 0.5, ...)
  for(j in 1:N) { for(k in 1:N) { 
    rect(j - 0.5, k - 0.5, j + 0.5, k + 0.5, col =mat[j,k])
  }}
  invisible(0)
}
