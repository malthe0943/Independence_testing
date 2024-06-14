generate_non_monotonic_samples_M2 <- function(n, rho) {
  # Generate X1 from U([-1, 1])
  X1 <- runif(n, -1, 1)
  
  # Generate epsilon from N(0, 1)
  epsilon <- rnorm(n)
  
  # Generate Y1 using the model Y1 = |X1| / rho + epsilon
  Y1 <- (abs(X1)^rho)*epsilon
  
  # Generate X2 and Y2 independently from uniform distribution
  X2 <- runif(n)
  Y2 <- runif(n)
  
  # Combine into X and Y
  X <- cbind(X1, X2)
  Y <- cbind(Y1, Y2)
  
  return(list(X = X, Y = Y))
}

# Example usage
set.seed(42) # For reproducibility
result_M2 <- generate_non_monotonic_samples_M2(n = 100, rho = 0.5)

x <- result_M2$X
y <- result_M2$Y

