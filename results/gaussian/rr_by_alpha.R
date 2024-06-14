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

data_gen_method <- "gaussian"
n_obs <- 100 # 200
# n_dim <- 4 # 1
n_sim <- 40 # 100
n_perm <- 40 # 100
n_reps <- 5 # 5

rho1 <- 0.5*2/5
n_dim1 <- 4 
rho2 <- 0.2*2/5
n_dim2 <- 10
rho3 <- 0.066*2/5
n_dim3 <- 30
rho4 <- 0.02*2/5
n_dim4 <- 100

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

res_alt1 <- calibrate_gaussian(data_gen_method, n_obs, n_dim1, n_sim, n_perm, n_reps, rho = rho1)
res_alt2 <- calibrate_gaussian(data_gen_method, n_obs, n_dim2, n_sim, n_perm, n_reps, rho = rho2)
res_alt3 <- calibrate_gaussian(data_gen_method, n_obs, n_dim3, n_sim, n_perm, n_reps, rho = rho3)
res_alt4 <- calibrate_gaussian(data_gen_method, n_obs, n_dim4, n_sim, n_perm, n_reps, rho = rho4)

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

plot_rr_gauss(res_alt1, n_dim1, rho1)
plot_rr_gauss(res_alt2, n_dim2, rho2)
plot_rr_gauss(res_alt3, n_dim3, rho3)
plot_rr_gauss(res_alt4, n_dim4, rho4)


