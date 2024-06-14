generate_gaussian_samples <- function(p = 2, q = 2, rho = 0.5, I = NULL, n_obs = 10) {
  # Ensure I is specified, otherwise default to 1:p
  if (is.null(I)) {
    I <- 1:p
  }
  
  # Create the covariance matrix
  d <- p + q
  Gamma_rho <- matrix(0, nrow = d, ncol = d)

  # Fill in the covariance matrix
  for (i in 1:d) {
    for (j in 1:d) {
      if (i == j) {
        Gamma_rho[i, j] <- 1
      } else if ((i <= p) && (j > p) && ((i %in% I) && ((j - p) %in% 1:q))) {
        Gamma_rho[i, j] <- rho
        Gamma_rho[j, i] <- rho
      }
    }
  }
  
  # Generate the mean vector
  e_d <- rep(0, d)
  
  # Generate the samples from the multivariate normal distribution
  samples <- mvrnorm(n = n_obs, mu = e_d, Sigma = Gamma_rho)
  
  # Split the samples into X and Y
  X <- samples[, 1:p]
  Y <- samples[, (p + 1):(p + q)]
  
  return(list(X = X, Y = Y))
}

# Example usage
# set.seed(42) # For reproducibility
# result <- generate_gaussian_samples(p = 50, q = 50, rho = 0.02, n_obs = 5)
# 
# For p=q=2, max rho is 0.5
# For p=q=5, max rho is 0.2
# For p=q=15, max rho is 0.066
# For p=q=50, max rho is 0.02


