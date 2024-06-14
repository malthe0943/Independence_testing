library(miceadds)
library(mydHSIC)
library(Rcpp)
library(MASS)
library(doParallel)
library(energy)
library(dHSIC)
library(ggplot2)
library(stringr)

start_time <- Sys.time()

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

set.seed(42)

# List of kernels
gauss_kerns <- c("dCor", "linear", "gaussian", "gaussian_opt")
gauss_kern_names <- c("dCor", "linear", "gaussian", "gaussian_opt")

data_gen_method <- "gaussian"
n_obs <- 10 # 200
dims <- c(4, 10, 30, 100)
n_sim <- 5 # 100
n_perm <- 5 # 100
n_reps <- 5 # 100

# Register parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("generate_gaussian_samples",
                    "mvrnorm",
                    "simulate_data", 
                    "calculate_statistic",
                    "statistics",
                    "dcor",
                    "dhsic"))

res_gauss_d <- analyse_d(data_gen_method, n_obs, dims, n_sim, n_perm, n_reps)

stopCluster(cl)

end_time <- Sys.time()

plot_analysis_d(res_gauss_d, data_gen_method)

end_time - start_time

