
# Linear kernel
linear_kernel <- function(x_1, x_2) {
  return(sum(x_1 * x_2))
}

# Gaussian kernels with range of sigma values
gauss_kern_1 <- function(x_1, x_2) {
  return(exp(-sum((x_1 - x_2)^2) / (2 * 0.01^2)))
}

gauss_kern_2 <- function(x_1, x_2) {
  return(exp(-sum((x_1 - x_2)^2) / (2 * 0.1^2)))
}

gauss_kern_3 <- function(x_1, x_2) {
  return(exp(-sum((x_1 - x_2)^2) / (2 * 0.3^2)))
}

gauss_kern_4 <- function(x_1, x_2) {
  return(exp(-sum((x_1 - x_2)^2) / (2 * 1^2)))
}

gauss_kern_5 <- function(x_1, x_2) {
  return(exp(-sum((x_1 - x_2)^2) / (2 * 10^2)))
}

# Gaussian kernel with optimised sigma - REPLACE SIGMA WITH THE OPTIMISED VALUE
gaussian_kernel <- function(x_1, x_2) {
  return(exp(-sum((x_1 - x_2)^2) / (2 * sigma^2)))
}
