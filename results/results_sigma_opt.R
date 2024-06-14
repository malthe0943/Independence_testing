library(miceadds)
library(mydHSIC)
library(Rcpp)
library(MASS)
library(doParallel)
library(energy)
library(dHSIC)
library(ggplot2)
library(stringr)

eta_range <- c(0.25, 0.5, 0.75, 1)

Sys.time()

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

set.seed(42)

data_gen_method <- "gaussian"
n_obs <- 20 # 200
n_sim <- 3 # 100
n_perm <- 3 # 100
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

res_gauss1 <- analyse_d(data_gen_method, d_range1, n_obs, n_dim1, n_sim, n_perm, n_reps)
res_gauss2 <- analyse_d(data_gen_method, d_range2, n_obs, n_dim2, n_sim, n_perm, n_reps)
res_gauss3 <- analyse_d(data_gen_method, d_range3, n_obs, n_dim3, n_sim, n_perm, n_reps)
res_gauss4 <- analyse_d(data_gen_method, d_range4, n_obs, n_dim4, n_sim, n_perm, n_reps)

stopCluster(cl)

Sys.time()

plot_analysis(res_gauss1, data_gen_method)
plot_analysis(res_gauss2, data_gen_method)
plot_analysis(res_gauss3, data_gen_method)
plot_analysis(res_gauss4, data_gen_method)
