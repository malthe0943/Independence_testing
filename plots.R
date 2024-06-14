source("rr_by_p.R")

# Hyperparameters
statistics <- c(
                "dCOR",
                "dHSIC with gaussian kernel",
                "dHSIC with gaussian kernel and optimized sigma",
                "dHSIC with linear kernel"
                )
# hyperparameters <- NULL
hyperparameters <- list(
  dHSIC = list(sigma = 0.01))  # Example value for sigma

num_observations <- 200
dimensions <- c(1, 5, 20, 50, 100)
data_gen_method <- "independent" #c("independent", "linear","normal", "M1", "M2", "M3")

# Plots for calibration of the rejection rate vs alpha
experiment <- Experiment$new(data_gen_method, statistics, hyperparameters, num_observations, dimensions)
alpha_check_results <- experiment$run_alpha_check()
experiment$plot_alpha_check(alpha_check_results)







# Plots for Rejection Rate by increasing p

experiment$run_all_experiments()
experiment$plot_results()

experiment$rejection_rates


