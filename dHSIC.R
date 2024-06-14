# library(kernlab)
library(dHSIC)

# Function to compute median of data for sigma
compute_sigma <- function(data) {
  n <- nrow(data)
  dists <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      dists[i, j] <- sum((data[i, ] - data[j, ])^2)
    }
  }
  return(sqrt(median(dists)))
}

# Function to compute the Gaussian kernel matrix with defined sigma
gauss_kern <- function(data, sigma) {
  n <- nrow(data)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- exp(-sum((data[i, ] - data[j, ])^2) / (2 * sigma^2))
    }
  }
  return(K)
}

# Function to compute the Gaussian kernel matrix with median heuristic sigma
gauss_kern_median <- function(data) {
  sigma <- compute_sigma(data)
  return(gauss_kern(data, sigma))
}

# Function to compute the linear kernel matrix
linear_kern <- function(data) {
  return(data %*% t(data))
}

# Function to compute dHSIC, with cpp kernels
compute_dHSIC <- function(data, hyperparameters) {
  # Assuming data is a list of matrices, each matrix representing a variable set
  if (!is.list(data)) stop("Data should be a list of matrices")
  
  n <- nrow(data[[1]])
  m <- length(data)
  
  # Initialize kernel matrices list
  K_list <- vector("list", m)
  
  # Compute kernel matrices for each variable set
  for (i in 1:m) {
    if (hyperparameters$kernel == "linear_cpp") {
      K_list[[i]] <- cpp_linearKern(data[[i]])
    } 
    else if (hyperparameters$kernel == "linear") {
      K_list[[i]] <- linear_kern(data[[i]])
    } 
    else if (hyperparameters$kernel == "gaussian_cpp") {
      K_list[[i]] <- cpp_gaussKern(data[[i]], hyperparameters$sigma)
    } 
    else if (hyperparameters$kernel == "gaussian") {
      K_list[[i]] <- gauss_kern(data[[i]], hyperparameters$sigma)
    } 
    else if (hyperparameters$kernel == "gaussian_median") {
      K_list[[i]] <- gauss_kern_median(data[[i]])
    } 
    else if (hyperparameters$kernel == "gaussian_median_cpp") {
      sigma <- compute_sigma(data[[i]])
      K_list[[i]] <- cpp_gaussKern(data[[i]], sigma)
    } 
    else {
      stop("Unknown kernel type")
    }
  }
  
  # Centering matrix
  H <- diag(n) - (1/n) * matrix(1, n, n)
  
  # Compute the centered kernel matrices
  K_centered <- lapply(K_list, function(K) H %*% K %*% H)

  # Compute the dHSIC value
  HSIC_value <- 0
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      HSIC_value <- HSIC_value + sum(K_centered[[i]] * K_centered[[j]])
    }
  }
  HSIC_value <- HSIC_value / (n - 1)^2
  
  return(HSIC_value)
}

