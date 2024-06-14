#' Compute the Gaussian kernel density estimate using the median rule
#' 
#' This function computes the Gaussian kernel density estimate using the median rule.
#' 
#' @name gauss_kern_median
#' @param data A matrix of data points.
#' @return The Gaussian kernel density estimate.
#' @export

gauss_kern_median <-
function(data) {
  sigma <- compute_sigma(data)
  return(gauss_kern(data, sigma))
}
