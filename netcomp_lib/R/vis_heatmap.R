
## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (plot_probmatrix)
#' Plots a heatmap of a probability matrix
#' 
#' @param pm Probability matrix
#' @param add_legend T/F: If TRUE: add a legend plot
#' @param is_dm T/F: If TRUE, input matrix is on a [-1,1] scale (as a difference of probability matrices)
#' @param only_legend T/F: If TRUE: ONLY plot the legend
#' @param colors Colors to use in the heat map (as a vector of hex values)
#' @param ... temp
#' 
#' @return Nothing; produces a plot (standard R plot functions)
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

