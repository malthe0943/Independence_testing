# Define the statistics and their names
statistics <- list("dCor", "linear", "gaussian", "gaussian_opt")

calibrate <- function( data_gen_method,
                       n_obs = 200, 
                       n_dim = 1, 
                       n_sim = 100, 
                       n_perm = 100, 
                       n_reps = 5,
                       sigma = 0.01) {
  
  alpha_values <- seq(0, 1, by = 0.05)
  results <- data.frame()
  
  total_iterations <- length(statistics) * length(alpha_values) * n_reps
  
  pb <- txtProgressBar(max = total_iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  
  progress_env <- new.env()
  progress_env$progress <- 0
  
  for (k in 1:length(statistics)) {
    stat_name <- statistics[[k]]
    rep_rejection_rates <- matrix(0, nrow = length(alpha_values), ncol = n_reps)
    
    for (rep in 1:n_reps) {
      rejection_rates <- foreach(i = seq_along(alpha_values), .combine = c 
                                 #.packages = c("dHSIC", "Rcpp", "mydHSIC"), 
      ) %dopar% {
        alpha <- alpha_values[i]
        rejection_count <- 0
        for (j in 1:n_sim) {
          data <- simulate_data(data_gen_method = data_gen_method, 
                                n_obs = n_obs,
                                n_dim = n_dim)
          x <- data$x
          y <- data$y
          
          observed_statistic <- calculate_statistic(x, y, statistic = statistics[k], sigma = sigma)
          
          # Permutation test
          permuted_statistics <- numeric(n_perm)
          for (p in 1:n_perm) {
            y_permuted <- as.matrix(sample(y))
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
      
      progress_env$progress <- progress_env$progress + length(alpha_values)
      setTxtProgressBar(pb, progress_env$progress)
    }
    
    mean_rejection_rates <- rowMeans(rep_rejection_rates)
    se_rejection_rates <- apply(rep_rejection_rates, 1, sd) / sqrt(n_reps)
    
    temp_results <- data.frame(alpha = alpha_values, mean_rejection_rate = mean_rejection_rates, se_rejection_rate = se_rejection_rates, statistic = stat_name)
    results <- rbind(results, temp_results)
  }
  return(results)
}

