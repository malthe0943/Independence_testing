library(miceadds)
library(mydHSIC)
library(Rcpp)
library(MASS)
library(doParallel)
library(energy)
library(dHSIC)
library(ggplot2)

Sys.time()

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

set.seed(42)

data_gen_method <- "gaussian"
n_obs <- 5 # 200
n_dim <- 1 # 1
n_sim <- 5 # 100
n_perm <- 5 # 100
n_reps <- 5 # 5

# Register parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("generate_gaussian_samples",
                    "mvrnorm",
                    "generate_non_monotonic_samples_M1",
                    "generate_non_monotonic_samples_M2",
                    "generate_non_monotonic_samples_M3",
                    "simulate_data", 
                    "calculate_statistic",
                    "statistics",
                    "dcor",
                    "dhsic"))

res_alt1 <- calibrate_gaussian(data_gen_method, n_obs, n_dim, n_sim, n_perm, n_reps)

stopCluster(cl)

Sys.time()

plot_rr(res_alt1)


