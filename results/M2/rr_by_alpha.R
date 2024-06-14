library(miceadds)
library(mydHSIC)
library(Rcpp)
library(MASS)
library(doParallel)
library(energy)
library(dHSIC)
library(ggplot2)

start_time <- Sys.time()

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

set.seed(42)

data_gen_method <- "M2"
n_obs <- 100 # 200
n_dim <- 1 # 1
n_sim <- 40 # 100
n_perm <- 40 # 100
n_reps <- 5 # 5
rho <- 0.8

# Register parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("generate_gaussian_samples",
                    "mvrnorm",
                    "generate_non_monotonic_samples_M2",
                    "simulate_data", 
                    "calculate_statistic",
                    "statistics",
                    "dcor",
                    "dhsic"))

res_alt1 <- calibrate_M123(data_gen_method, n_obs, n_dim, n_sim, n_perm, n_reps, rho = rho)

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

plot_rr_M2(res_alt1, rho)


