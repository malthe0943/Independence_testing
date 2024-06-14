generate_non_monotonic_samples_M3 <- function(n, a) {
  # Define the joint density function
  joint_density <- function(x, y, a) {
    return((1 / (4 * pi^2)) * (1 + sin(a * x) * sin(a * y)))
  }
  
  # Initialize storage for accepted samples
  accepted_X <- numeric(n)
  accepted_Y <- numeric(n)
  num_accepted <- 0
  
  # Generate more samples than needed to ensure we get enough accepted ones
  oversample_factor <- 10
  while (num_accepted < n) {
    num_to_generate <- (n - num_accepted) * oversample_factor
    x_candidates <- runif(num_to_generate, -pi, pi)
    y_candidates <- runif(num_to_generate, -pi, pi)
    u_candidates <- runif(num_to_generate)
    
    densities <- joint_density(x_candidates, y_candidates, a)
    accepted_indices <- which(u_candidates < densities / 2)
    
    num_new_accepted <- min(length(accepted_indices), n - num_accepted)
    if (num_new_accepted > 0) {
      accepted_X[(num_accepted + 1):(num_accepted + num_new_accepted)] <- x_candidates[accepted_indices[1:num_new_accepted]]
      accepted_Y[(num_accepted + 1):(num_accepted + num_new_accepted)] <- y_candidates[accepted_indices[1:num_new_accepted]]
      num_accepted <- num_accepted + num_new_accepted
    }
  }
  
  # Generate X2 and Y2 independently from uniform distribution
  X2 <- runif(n)
  Y2 <- runif(n)
  
  # Combine into X and Y
  X_full <- cbind(accepted_X, X2)
  Y_full <- cbind(accepted_Y, Y2)
  
  return(list(X = X_full, Y = Y_full))
}

# Example usage
# set.seed(123) # For reproducibility
# result_M3 <- generate_non_monotonic_samples_M3(n = 100, a = 5)
# print(result_M3)
