#' Linear kernel function
#' 
#' This function computes the linear kernel.
#' 
#' @name linear_kern
#' @param data A matrix of data points.
#' @return The linear kernel.
#' @export

linear_kern <-
function(data) {
  return(data %*% t(data))
}
