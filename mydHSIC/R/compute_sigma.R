#' Compute the sigma parameter for the RBF kernel.
#' 
#' This function computes the sigma parameter for the RBF kernel.
#' 
#' @name compute_sigma
#' @param data A matrix of data points.
#' @return The sigma parameter.
#' @export

compute_sigma <-
function(data) {
  n <- nrow(data)
  dists <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      dists[i, j] <- sum((data[i, ] - data[j, ])^2)
    }
  }
  return(sqrt(median(dists)))
}
