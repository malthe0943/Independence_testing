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

rho_range <- seq(0.1, 1, 0.1)
data_gen_method <- "M2"
n_obs <- 200 # 200
n_dim <- 4
n_sim <- 40 # 100
n_perm <- 40 # 100
n_reps <- 5 # 5

# Register parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("generate_non_monotonic_samples_M2",
                    "simulate_data", 
                    "calculate_statistic",
                    "statistics",
                    "dcor",
                    "dhsic"))

res_M2 <- analyse_rho(data_gen_method, rho_range, n_obs, n_dim, n_sim, n_perm, n_reps)

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

plot_analysis_rho_M2(res_M2, data_gen_method)
