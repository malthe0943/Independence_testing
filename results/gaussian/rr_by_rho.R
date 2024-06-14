library(miceadds)
library(mydHSIC)
library(Rcpp)
library(MASS)
library(doParallel)
library(energy)
library(dHSIC)
library(ggplot2)
library(stringr)

# For p=q=2, max rho is 0.5
# For p=q=5, max rho is 0.2
# For p=q=15, max rho is 0.066
# For p=q=50, max rho is 0.02

n_dim1 <- 4
rho_range1 <- seq(0, 0.5, length.out = 21)

n_dim2 <- 10
rho_range2 <- seq(0, 0.2, length.out = 21)

n_dim3 <- 30
rho_range3 <- seq(0, 0.066, length.out = 21)

n_dim4 <- 100
rho_range4 <- seq(0, 0.02, length.out = 21)

start_time <- Sys.time()

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

set.seed(42)

data_gen_method <- "gaussian"
n_obs <- 100 # 200
n_sim <- 35 # 100
n_perm <- 35 # 100
n_reps <- 5 # 5

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

res_gauss1 <- analyse_rho(data_gen_method, rho_range1, n_obs, n_dim1, n_sim, n_perm, n_reps)
res_gauss2 <- analyse_rho(data_gen_method, rho_range2, n_obs, n_dim2, n_sim, n_perm, n_reps)
res_gauss3 <- analyse_rho(data_gen_method, rho_range3, n_obs, n_dim3, n_sim, n_perm, n_reps)
res_gauss4 <- analyse_rho(data_gen_method, rho_range4, n_obs, n_dim4, n_sim, n_perm, n_reps)

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

plot_analysis_rho(res_gauss1, data_gen_method, n_dim1)
plot_analysis_rho(res_gauss2, data_gen_method, n_dim2)
plot_analysis_rho(res_gauss3, data_gen_method, n_dim3)
plot_analysis_rho(res_gauss4, data_gen_method, n_dim4)
