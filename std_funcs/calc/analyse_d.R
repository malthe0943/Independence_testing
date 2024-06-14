# Define the statistics and their names
statistics <- list("dCor", "linear", "gaussian", "gaussian_opt")

analyse_d <- function(data_gen_method,
                      n_obs = 200, 
                      dims, 
                      n_sim = 100, 
                      n_perm = 100, 
                      n_reps = 5) {
  
  alpha <- 0.05
  results <- data.frame()
  
  rho_range <- c(0.5, 0.2, 0.066, 0.02)*2/5

  total_iterations <- length(gauss_kerns) * length(dims) * n_reps
  
  pb <- txtProgressBar(max = total_iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  
  progress_env <- new.env()
  progress_env$progress <- 0
  
  for (k in 1:length(gauss_kerns)) {
    kernel_name <- gauss_kern_names[k]
    rep_rejection_rates <- matrix(0, nrow = length(dims), ncol = n_reps)
    
    for (rep in 1:n_reps) {
      rejection_rates <- foreach(i = seq_along(dims), .combine = c 
                                 #.packages = c("dHSIC", "Rcpp", "mydHSIC"), 
      ) %dopar% {
        dim <- dims[i]
        sigma <- 10*dim
        rejection_count <- 0
        for (j in 1:n_sim) {
          # Simulate data, linear dependence
          data <- simulate_data(data_gen_method = data_gen_method, 
                                n_obs = n_obs,
                                n_dim = dim,
                                rho = rho_range[i])
          x <- data$x
          y <- data$y
          
          observed_statistic <- calculate_statistic(x, y, statistic = statistics[k], sigma = sigma)
          
          # Permutation test
          permuted_statistics <- numeric(n_perm)
          for (p in 1:n_perm) {
            y_permuted <- apply(y, 2, sample)
            permuted_statistics[p] <- calculate_statistic(x, y_permuted, statistic = statistics[k], sigma = sigma)
          }
          
          # Calculate the critical value
          critical_value <- quantile(permuted_statistics, 1 - alpha)
          
          # Check if the observed test statistic exceeds the critical value
          if (observed_statistic > critical_value) {
            rejection_count <- rejection_count + 1
          }
        }
        
        # Calculate the rejection rate for this alpha
        rejection_rate <- rejection_count / n_sim
        rejection_rate
      }
      
      rep_rejection_rates[, rep] <- rejection_rates
      
      progress_env$progress <- progress_env$progress + length(dims)
      setTxtProgressBar(pb, progress_env$progress)
    }
    
    mean_rejection_rates <- rowMeans(rep_rejection_rates)
    se_rejection_rates <- apply(rep_rejection_rates, 1, sd) / sqrt(n_reps)
    
    temp_results <- data.frame(dim = dims, mean_rejection_rate = mean_rejection_rates, se_rejection_rate = se_rejection_rates, statistic = kernel_name)
    results <- rbind(results, temp_results)
  }
  return(results)
}



