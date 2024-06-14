library(GiniDistance)
# Function to compute dCOR
compute_dCOR <- function(x, y, alpha = 1) {
  if (!is.matrix(x)) stop("x should be a matrix")
  if (!is.vector(y)) stop("y should be a vector")
  if (!is.numeric(alpha) || alpha <= 0 || alpha > 2) stop("alpha should be in (0, 2]")
  
  n <- nrow(x)
  
  # Compute pairwise Euclidean distances raised to the power alpha for x and y
  a <- as.matrix(dist(x)^alpha)
  b <- as.matrix(dist(y)^alpha)
  
  # Double-centering
  A <- a - rowMeans(a) - matrix(rep(colMeans(a), each = n), n, n) + mean(a)
  B <- b - rowMeans(b) - matrix(rep(colMeans(b), each = n), n, n) + mean(b)
  
  # Compute dCov and dCor
  dCov <- sum(A * B) / (n^2)
  dVarX <- sum(A * A) / (n^2)
  dVarY <- sum(B * B) / (n^2)
  
  dCor <- dCov / sqrt(dVarX * dVarY)
  
  return(dCor)
}

# Function to compute dCOR
compute_dCOR <- function(x, y, alpha = 1) {
  if (!is.matrix(x)) stop("x should be a matrix")
  if (!is.vector(y)) stop("y should be a vector")
  if (!is.numeric(alpha) || alpha <= 0 || alpha > 2) stop("alpha should be in (0, 2]")
  
  n <- nrow(x)
  
  # Compute pairwise Euclidean distances raised to the power alpha for x and y
  a <- as.matrix(dist(x)^alpha)
  b <- as.matrix(dist(y)^alpha)
  
  # Double-centering
  H <- diag(n) - (1/n) * matrix(1, n, n)
  A <- H %*% a %*% H
  B <- H %*% b %*% H
  
  # Compute dCov and dCor
  dCov <- sum(A * B) / (n^2)
  dVarX <- sum(A * A) / (n^2)
  dVarY <- sum(B * B) / (n^2)
  
  dCor <- dCov / sqrt(dVarX * dVarY)
  
  return(dCor)
}

# # Independent
# print("Independent")
# set.seed(42)
# x <- matrix(rnorm(500),ncol=1)
# y <- rbinom(500,30,0.1)
# cat("dCor library:", dCor(x = x, y = y, alpha = 1), "\n")
# cat("my dCOR:", compute_dCOR(x = x, y = y, alpha = 1), "\n")
# 
# # Dependent
# print("Dependent")
# set.seed(42)
# x <- matrix(rnorm(500),ncol=1)
# y <- as.vector(x)+0.02*rnorm(500)
# cat("dCor library:", dCor(x = x, y = y, alpha = 1), "\n")
# cat("my dCOR:", compute_dCOR(x = x, y = y, alpha = 1), "\n")



