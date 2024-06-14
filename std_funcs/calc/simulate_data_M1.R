generate_non_monotonic_samples_M1 <- function(n, A) {
  # Generate random values for A, Theta, and epsilon
  Theta <- runif(n, 0, 2 * pi)
  epsilon1 <- rnorm(n)
  epsilon2 <- rnorm(n)
  
  # Generate X1 and Y1
  X1 <- A * cos(Theta) + epsilon1 / 4
  Y1 <- A * sin(Theta) + epsilon2 / 4
  
  # Generate X2 and Y2 independently from uniform distribution
  X2 <- runif(n)
  Y2 <- runif(n)
  
  # Combine into X and Y
  X <- cbind(X1, X2)
  Y <- cbind(Y1, Y2)
  
  return(list(X = X, Y = Y))
}

# # Example usage
# set.seed(42) # For reproducibility
# result_M1 <- generate_non_monotonic_samples_M1(n = 100)
# print(result_M1)
# 
# x <- result_M1$X
# y <- result_M1$Y




