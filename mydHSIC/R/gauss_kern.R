gauss_kern <-
function(data, sigma) {
  n <- nrow(data)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- exp(-sum((data[i, ] - data[j, ])^2) / (2 * sigma^2))
    }
  }
  return(K)
}
