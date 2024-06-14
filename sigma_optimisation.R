library(dHSIC)
library(ggplot2)
library(profvis)
library(foreach)
library(doParallel)
library(doSNOW)
library(microbenchmark)
library(progress)
library(inline)
library(Rcpp)
library(mydHSIC)

optimise_sigma_par <- function(n_obs = 200, n_dim = 1, num_sim = 100, num_perm = 100, num_reps = 5) {
  set.seed(42)
  
  alpha_values <- seq(0, 1, by = 0.05)
  results <- data.frame()
  
  total_iterations <- length(gauss_kerns) * length(alpha_values) * num_reps
  
  pb <- txtProgressBar(max = total_iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)

  progress_env <- new.env()
  progress_env$progress <- 0
  
  for (k in 1:length(gauss_kerns)) {
    kernel_name <- gauss_kern_names[k]
    rep_rejection_rates <- matrix(0, nrow = length(alpha_values), ncol = num_reps)
    
    for (rep in 1:num_reps) {
      rejection_rates <- foreach(i = seq_along(alpha_values), .combine = c 
                                 #.packages = c("dHSIC", "Rcpp", "mydHSIC"), 
                                 ) %dopar% {
        alpha <- alpha_values[i]
        rejection_count <- 0
        for (j in 1:num_sim) {
          # Simulate data, independent
          x <- matrix(rnorm(n_obs * n_dim), ncol = n_dim)
          y <- matrix(rnorm(n_obs), ncol = 1)
          
          if (gauss_kerns[k] == "gaussian") {
          observed_statistic <- dhsic(list(x, y), kernel = "gaussian")$dHSIC
          } else {
          observed_statistic <- mydHSIC::compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = gauss_kerns[[k]]))
          }
          # Permutation test
          permuted_statistics <- numeric(num_perm)
          for (p in 1:num_perm) {
            y_permuted <- as.matrix(sample(y))
            if (gauss_kerns[k] == "gaussian") {
              permuted_statistics[p] <- dhsic(list(x, y_permuted), kernel = "gaussian")$dHSIC
            } else {
              permuted_statistics[p] <- mydHSIC::compute_dHSIC(list(x, y_permuted), list(kernel = "gaussian_cpp", sigma = gauss_kerns[[k]]))
            }
          }
          
          # Calculate the critical value
          critical_value <- quantile(permuted_statistics, 1 - alpha)
          
          # Check if the observed test statistic exceeds the critical value
          if (observed_statistic > critical_value) {
            rejection_count <- rejection_count + 1
          }
        }
        
        # Calculate the rejection rate for this alpha
        rejection_rate <- rejection_count / num_sim
        rejection_rate
      }
      
      rep_rejection_rates[, rep] <- rejection_rates
      
      progress_env$progress <- progress_env$progress + length(alpha_values)
      setTxtProgressBar(pb, progress_env$progress)
    }
    
    mean_rejection_rates <- rowMeans(rep_rejection_rates)
    se_rejection_rates <- apply(rep_rejection_rates, 1, sd) / sqrt(num_reps)
    
    temp_results <- data.frame(alpha = alpha_values, mean_rejection_rate = mean_rejection_rates, se_rejection_rate = se_rejection_rates, kernel = kernel_name)
    results <- rbind(results, temp_results)
  }
  return(results)
}

optimise_sigma <- function(n_obs = 200, n_dim = 1, num_sim = 100, num_perm = 100, num_reps = 5) {
  set.seed(42)
  
  alpha_values <- seq(0, 1, by = 0.05)
  results <- data.frame()
  
  for (k in 1:length(gauss_kerns)) {
    kernel_name <- gauss_kern_names[k]
    rep_rejection_rates <- matrix(0, nrow = length(alpha_values), ncol = num_reps)
    
    for (rep in 1:num_reps) {
    rejection_rates <- numeric(length(alpha_values))

    for (i in seq_along(alpha_values)) {
      alpha <- alpha_values[i]
      if ((i - 1) %% 5 == 0) {
        cat("Checking alpha consistency for", kernel_name, "and rep =", rep,
            "and alpha =", alpha, "\n")
      }
      rejection_count <- 0

      for (j in 1:num_sim) {
        # Simulate data, independent
        x <- matrix(rnorm(n_obs * n_dim), ncol = n_dim)
        y <- matrix(rnorm(n_obs), ncol = 1)

        if (gauss_kerns[k] == "gaussian") {
          observed_statistic <- dhsic(list(x, y), kernel = "gaussian")$dHSIC
        } else {
          observed_statistic <- mydHSIC::compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = gauss_kerns[[k]]))
        }
        # Permutation test
        permuted_statistics <- numeric(num_perm)
        for (p in 1:num_perm) {
          y_permuted <- as.matrix(sample(y))
          if (gauss_kerns[k] == "gaussian") {
            permuted_statistics[p] <- dhsic(list(x, y_permuted), kernel = "gaussian")$dHSIC
          } else {
            permuted_statistics[p] <- mydHSIC::compute_dHSIC(list(x, y_permuted), list(kernel = "gaussian_cpp", sigma = gauss_kerns[[k]]))
          }
        }

        # Calculate the critical value
        critical_value <- quantile(permuted_statistics, 1 - alpha)

        # Check if the observed test statistic exceeds the critical value
        if (observed_statistic > critical_value) {
          rejection_count <- rejection_count + 1
        }
      }

      # Calculate the rejection rate for this alpha
      rejection_rate <- rejection_count / num_sim
      rejection_rates[i] <- rejection_rate
    }

    rep_rejection_rates[, rep] <- rejection_rates
  }
    mean_rejection_rates <- rowMeans(rep_rejection_rates)
    se_rejection_rates <- apply(rep_rejection_rates, 1, sd) / sqrt(num_reps)
    
    temp_results <- data.frame(alpha = alpha_values, mean_rejection_rate = mean_rejection_rates, se_rejection_rate = se_rejection_rates, kernel = kernel_name)
    results <- rbind(results, temp_results)
    #cat("Checking alpha consistency for kernel =", kernel_name, ", repetition =", rep, " and alpha =", alpha, "\n")
    
  }
  return(results)
}

##################################################################
#
#     Under H0
#
##################################################################

# Register parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the kernels and their names
gauss_kerns <- list("gaussian", 0.003, 0.01, 0.1, 0.3, 1, 10)
gauss_kern_names <- c("default", "sigma_0.003", "sigma_0.01", "sigma_0.1", "sigma_0.3", "sigma_1", "sigma_10")

# Export necessary variables and functions to each worker node
clusterExport(cl, c("gauss_kerns", "dhsic"))

# Running test function with parallelisation
results <- optimise_sigma_par(n_obs = 10, n_dim = 1, num_sim = 100, num_perm = 100, num_reps = 5)

# Stop the parallel backend
stopCluster(cl)

# Function to plot rejection rates
plot_rejection_rates <- function(results, n) {
  ggplot(results, aes(x = alpha, y = mean_rejection_rate, color = kernel, fill = kernel)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_rejection_rate - se_rejection_rate, ymax = mean_rejection_rate + se_rejection_rate), alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste("Rejection Rates by Alpha with SE Bands, n =", n), x = "Significance Level (alpha)", y = "Rejection Rate") +
    theme_minimal()
}

plot_rejection_rates(results, 200)

##################################################################
#                                                                #
#     Under H1                                                   #
#                                                                #
##################################################################

optimise_sigma_eta <- function(n_obs = 200, dims, eta, num_sim = 100, num_perm = 100, num_reps = 5) {
  set.seed(42)
  
  alpha <- 0.05
  results <- data.frame()
  
  total_iterations <- length(gauss_kerns) * length(dims) * num_reps
  
  pb <- txtProgressBar(max = total_iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  
  progress_env <- new.env()
  progress_env$progress <- 0
  
  for (k in 1:length(gauss_kerns)) {
    kernel_name <- gauss_kern_names[k]
    rep_rejection_rates <- matrix(0, nrow = length(dims), ncol = num_reps)
    
    for (rep in 1:num_reps) {
      rejection_rates <- foreach(i = seq_along(dims), .combine = c 
                                 #.packages = c("dHSIC", "Rcpp", "mydHSIC"), 
      ) %dopar% {
        dim <- dims[i]
        rejection_count <- 0
        for (j in 1:num_sim) {
          # Simulate data, linear dependence
          x <- matrix(rnorm(n_obs * dim), ncol = dim)
          y <- matrix(rnorm(n_obs), ncol = 1)
          y <- y + 0.5 * x[, 1]
          
          if (gauss_kerns[k] == "gaussian") {
            observed_statistic <- dhsic(list(x, y), kernel = "gaussian")$dHSIC
          } else {
            observed_statistic <- mydHSIC::compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = gauss_kerns[[k]]*(dim^eta)))
          }
          # Permutation test
          permuted_statistics <- numeric(num_perm)
          for (p in 1:num_perm) {
            y_permuted <- as.matrix(sample(y))
            if (gauss_kerns[k] == "gaussian") {
              permuted_statistics[p] <- dhsic(list(x, y_permuted), kernel = "gaussian")$dHSIC
            } else {
              permuted_statistics[p] <- mydHSIC::compute_dHSIC(list(x, y_permuted), list(kernel = "gaussian_cpp", sigma = gauss_kerns[[k]]*(dim^eta)))
            }
          }
          
          # Calculate the critical value
          critical_value <- quantile(permuted_statistics, 1 - alpha)
          
          # Check if the observed test statistic exceeds the critical value
          if (observed_statistic > critical_value) {
            rejection_count <- rejection_count + 1
          }
        }
        
        # Calculate the rejection rate for this alpha
        rejection_rate <- rejection_count / num_sim
        rejection_rate
      }
      
      rep_rejection_rates[, rep] <- rejection_rates
      
      progress_env$progress <- progress_env$progress + length(dims)
      setTxtProgressBar(pb, progress_env$progress)
    }
    
    mean_rejection_rates <- rowMeans(rep_rejection_rates)
    se_rejection_rates <- apply(rep_rejection_rates, 1, sd) / sqrt(num_reps)
    
    temp_results <- data.frame(dim = dims, mean_rejection_rate = mean_rejection_rates, se_rejection_rate = se_rejection_rates, kernel = kernel_name)
    results <- rbind(results, temp_results)
  }
  return(results)
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the kernels and their names
gauss_kerns <- list("gaussian", 0.01, 0.1, 1, 10, 100)
gauss_kern_names <- c("default", "sigma_0.01", "sigma_0.1", "sigma_1", "sigma_10", "sigma_100")

# Defining etas
eta1 <- 0.25
eta2 <- 0.5
eta3 <- 0.75
eta4 <- 1

# Export necessary variables and functions to each worker node
clusterExport(cl, c("gauss_kerns", "dhsic"))

# Running test function with parallelisation
results1 <- optimise_sigma_eta(n_obs = 100, dims = c(4, 10, 30, 100), eta1, 
                              num_sim = 30, num_perm = 30, num_reps = 5)
results2 <- optimise_sigma_eta(n_obs = 100, dims = c(4, 10, 30, 100), eta2, 
                              num_sim = 30, num_perm = 30, num_reps = 5)
results3 <- optimise_sigma_eta(n_obs = 100, dims = c(4, 10, 30, 100), eta3, 
                              num_sim = 30, num_perm = 30, num_reps = 5)
results4 <- optimise_sigma_eta(n_obs = 100, dims = c(4, 10, 30, 100), eta4, 
                              num_sim = 30, num_perm = 30, num_reps = 5)

# Stop the parallel backend
stopCluster(cl)

# Function to plot rejection rates
plot_rejection_rates_eta <- function(results, eta) {
  ggplot(results, aes(x = dim, y = mean_rejection_rate, color = kernel, fill = kernel)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_rejection_rate - se_rejection_rate, ymax = mean_rejection_rate + se_rejection_rate), alpha = 0.2) +
    scale_x_continuous(trans='log10', breaks=c(4, 10, 30, 100)) +
    labs(title = paste("Rejection Rates by d with SE Bands, eta =", eta), 
         x = "Dimension, d (log scale)", 
         y = "Rejection Rate") +
    theme_minimal()
}

plot_rejection_rates_eta(results1, eta1)
plot_rejection_rates_eta(results2, eta2)
plot_rejection_rates_eta(results3, eta3)
plot_rejection_rates_eta(results4, eta4)

# Conclude that eta doesn't matter much for sigma large enough.
# So for simplicity, we choose eta = 1.
# It appears that sigma = 10 is slightly better than the alternatives.

##################################################################
#
#     Testing time difference w and w/o parallelisation
#
##################################################################

gauss_kerns <- list("gaussian", 0.003, 0.01, 0.1, 0.3, 1, 10)
gauss_kern_names <- c("default", "sigma_0.003", "sigma_0.01", "sigma_0.1", "sigma_0.3", "sigma_1", "sigma_10")

par_start <- Sys.time()

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("gauss_kerns", "dhsic"))
results <- optimise_sigma_par(n_obs = 200, n_dim = 1, num_sim = 20, num_perm = 20, num_reps = 3)
stopCluster(cl)

par_end <- Sys.time()

single_start <- Sys.time()

results <- optimise_sigma(n_obs = 200, n_dim = 1, num_sim = 20, num_perm = 20, num_reps = 3)

single_end <- Sys.time()

par_end - par_start
single_end - single_start

# > par_end - par_start
# Time difference of 8.078227 mins
# > single_end - single_start
# Time difference of 53.79159 mins

# Parallelisation faster by around a factor 6

profvis({
  results <- optimise_sigma(n_obs = 50, n_dim = 1, num_sim = 20, num_perm = 20, num_reps = 5)
})










