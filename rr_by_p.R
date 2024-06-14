library(R6)
library(dHSIC)
library(mydHSIC)
library(energy)
library(doParallel)

# source("dHSIC.R")
# source("dCOR.R")
source("./misc/custom_kernel_functions.R")

IndividualExperiment <- R6Class("IndividualExperiment",
  public = list(
    data_gen_method = NULL,
    statistic = NULL,
    hyperparameters = NULL,
    num_observations = NULL,
    dimension = NULL,
    rejection_rate = NULL,
    
    initialize = function(data_gen_method, statistic, hyperparameters, num_observations, dimension) {
      self$data_gen_method <- data_gen_method
      self$statistic <- statistic
      self$hyperparameters <- hyperparameters
      self$num_observations <- num_observations
      self$dimension <- dimension
    },
    
    run = function() {
      self$rejection_rate <- self$rr_by_p(self$dimension, self$statistic)
    },

    simulate_data = function() {
        # Generate independent data
      if (self$data_gen_method == "independent") {
        x <- matrix(rnorm(self$num_observations * self$dimension), ncol = self$dimension)
        y <- matrix(rnorm(self$num_observations), ncol = 1)
        # Generate linearly dependent data
      } else if (self$data_gen_method == "linear") {
        x <- matrix(rnorm(self$num_observations * self$dimension), ncol = self$dimension)
        y <- matrix(rnorm(self$num_observations), ncol = 1)
        y <- y + 0.5 * x[, 1]
        # Generate normally distributed data
      } else if (self$data_gen_method == "normal") {
        # Generate data according to the M1 procedure
      } else if (self$data_gen_method == "M1") {
        # Generate data according to the M2 procedure
      } else if (self$data_gen_method == "M2") {
        # Generate data according to the M3 procedure
      } else if (self$data_gen_method == "M3") {
      } else {
        stop("Unknown data generation method")
      }
      list(x = x, y = y)
    },
    
    # Method to calculate statistics based on c("linear", "gauss", "gauss_opt_sigma", "dCOR")
    calculate_statistics = function(x, y) {
      if (self$statistic == "dCOR") {
        return(dcor(x, y))
      } 
      else if (self$statistic == "dHSIC with linear kernel") {
          return(mydHSIC::compute_dHSIC(list(x, y), list(kernel = "linear_cpp")))
      } 
      else if (self$statistic == "dHSIC with gaussian kernel") {
          return(dhsic(list(x, y), kernel = c("gaussian"))$dHSIC)
      } 
      else if (self$statistic == "dHSIC with gaussian kernel and optimized sigma") {
        return(mydHSIC::compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = self$hyperparameters$dHSIC$sigma)))
      } 
      else {
          stop("Unknown statistic")
      }
    },
    
    # Register parallel backend
    start_parallel = function() {
      num_cores <- detectCores() - 1
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      # Export necessary variables and functions to each worker node
      #clusterExport(cl, c())
    },
    
    stop_parallel = function() {
      stopCluster(cl)
    },
    
    rr_by_p = function(dimension, statistic) {
      set.seed(123)
      num_simulations <- 100  # Number of simulations to run
      num_permutations <- 100  # Number of permutations for the test
      alpha <- 0.05  # Significance level
      rejection_count <- 0
      
      for (i in 1:num_simulations) {
        if (i %% 10 == 0) cat("Dimension ", dimension,
                              ", statistic ", statistic, 
                              ", simulation ", i, "\n", sep = "")
        
        data <- self$simulate_data()
        x <- data$x
        y <- data$y

        observed_statistic <- self$calculate_statistics(x, y)
        
        # Permutation test (for consistency)
        permuted_statistics <- numeric(num_permutations)
        for (p in 1:num_permutations) {
          y_permuted <- as.matrix(sample(y))
          permuted_statistics[p] <- self$calculate_statistics(x, y_permuted)
        }
        
        # Calculate the critical value
        critical_value <- quantile(permuted_statistics, 1 - alpha)
        
        # Check if the observed test statistic exceeds the critical value
        if (observed_statistic > critical_value) {
          rejection_count <- rejection_count + 1
        }
      }
      
      # Calculate the rejection rate
      rejection_rate <- rejection_count / num_simulations
      return(rejection_rate)
    },
    
    # Rejection rate calculator using permutations
    check_alpha_consistency = function() {
      start_parallel()
      
      set.seed(42)
      
      alpha_values <- seq(0.05, 0.95, by = 0.05)
      rejection_rates <- numeric(length(alpha_values))

        for (i in seq_along(alpha_values)) {
        alpha <- alpha_values[i]
        cat("Checking alpha consistency for alpha =", alpha, "\n")
        
        rejection_count <- 0
        num_simulations <- 10 # 100  # Number of simulations to run
        num_permutations <- 10 # 100  # Number of permutations for the test
        
        for (j in 1:num_simulations) {
          data <- self$simulate_data()
          x <- data$x
          y <- data$y
          
          observed_statistic <- self$calculate_statistics(x, y)
          
          # Permutation test
          permuted_statistics <- numeric(num_permutations)
          for (p in 1:num_permutations) {
            y_permuted <- as.matrix(sample(y))
            permuted_statistics[p] <- self$calculate_statistics(x, y_permuted)
          }
          
          # Calculate the critical value
          critical_value <- quantile(permuted_statistics, 1 - alpha)
          
          # Check if the observed test statistic exceeds the critical value
          if (observed_statistic > critical_value) {
            rejection_count <- rejection_count + 1
          }
        }
        
        # Calculate the rejection rate for this alpha
        rejection_rate <- rejection_count / num_simulations
        rejection_rates[i] <- rejection_rate
      }
      
      return(list(alpha_values = alpha_values, rejection_rates = rejection_rates))
    }
  )
)

Experiment <- R6Class("Experiment",
  public = list(
    data_gen_method = NULL,
    statistics = NULL,
    hyperparameters = NULL,
    num_observations = NULL,
    dimension = NULL,
    rejection_rates = NULL,
    
    initialize = function(data_gen_method, statistics, hyperparameters, num_observations, dimension) {
      self$data_gen_method <- data_gen_method
      self$statistics <- statistics
      self$hyperparameters <- hyperparameters
      self$num_observations <- num_observations
      self$dimension <- dimension
      self$rejection_rates <- matrix(0, nrow = length(dimension), ncol = length(statistics))
      colnames(self$rejection_rates) <- statistics
      rownames(self$rejection_rates) <- dimension
    },
    
    run_all_experiments = function() {
      for (d in self$dimension) {
        for (s in self$statistics) {
          experiment <- IndividualExperiment$new(self$data_gen_method, s, self$hyperparameters, self$num_observations, d)
          experiment$run()
          self$rejection_rates[as.character(d), s] <- experiment$rejection_rate
        }
      }
    },
    
    plot_results = function() {
      matplot(as.numeric(rownames(self$rejection_rates)), self$rejection_rates, type = "b", pch = 1, col = 1:length(self$statistics), lty = 1:length(self$statistics),
              xlab = "Dimension (p)", ylab = "Rejection Rate", main = "Rejection Rate as a Function of Dimension")
      legend("topright", legend = self$statistics, col = 1:length(self$statistics), lty = 1:length(self$statistics), pch = 1)
    },
    
    run_alpha_check = function() {
      results <- list()
      for (s in self$statistics) {
          experiment <- IndividualExperiment$new(self$data_gen_method, s, self$hyperparameters, self$num_observations, 1)
          results[[paste0("stat_", s)]] <- experiment$check_alpha_consistency()# experiment$alpha_check()
      }
      
      return(results)
    },
    
    plot_alpha_check = function(alpha_check_results) {
      for (result_name in names(alpha_check_results)) {
        result <- alpha_check_results[[result_name]]
        main_title <- paste("Alpha Consistency Check for", result_name)
        plot(result$alpha_values, result$rejection_rates, type = "b", pch = 19, col = "blue",
             xlab = "Significance Level (alpha)", ylab = "Rejection Rate", main = main_title)
        abline(a = 0, b = 1, col = "red", lty = 2)  # Add a y = x reference line
        legend("topleft", legend = c("Observed Rejection Rates", "Ideal: y = x"), col = c("blue", "red"), lty = c(1, 2), pch = c(19, NA))
      }
    }
  )
)


