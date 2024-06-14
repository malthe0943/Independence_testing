simulate_data <- function(data_gen_method, 
                          n_obs = 200, 
                          n_dim = 1,
                          rho = 0.01) {
  # Generate independent data
  if (data_gen_method == "independent") {
    x <- matrix(rnorm(n_obs * n_dim), ncol = n_dim)
    y <- matrix(rnorm(n_obs), ncol = 1)
    # Generate linearly dependent data
  } else if (data_gen_method == "linear") {
    x <- matrix(rnorm(n_obs * n_dim), ncol = n_dim)
    y <- matrix(rnorm(n_obs), ncol = 1)
    y <- y + 0.15 * x[, 1]
    # Generate normally distributed data
  } else if (data_gen_method == "gaussian") {
    data <- generate_gaussian_samples(p = n_dim/2, q = n_dim/2, 
                                      rho = rho, n_obs = n_obs)
    x <- as.matrix(data$X)
    y <- as.matrix(data$Y)
    # Generate data according to the M1 procedure
  } else if (data_gen_method == "M1") {
    data <- generate_non_monotonic_samples_M1(n_obs, rho)
    x <- data$X
    y <- data$Y
    # Generate data according to the M2 procedure
  } else if (data_gen_method == "M2") {
    data <- generate_non_monotonic_samples_M2(n_obs, rho)
    x <- data$X
    y <- data$Y
    # Generate data according to the M3 procedure
  } else if (data_gen_method == "M3") {
    data <- generate_non_monotonic_samples_M3(n_obs, rho)
    x <- data$X
    y <- data$Y
  } else {
    stop("Unknown data generation method")
  }
  list(x = x, y = y)
}


