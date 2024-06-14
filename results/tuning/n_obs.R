library(miceadds)
library(mydHSIC)
library(Rcpp)

set.seed(42)

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

n_obs1 <- 10
n_obs2 <- 30
n_obs3 <- 70
n_obs4 <- 200

data_gen_method <- "independent"
n_dim <- 1
n_sim <- 40 # 100
n_perm <- 40 # 100
n_reps <- 5 # 5

start_time <- Sys.time()

# Register parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("simulate_data", 
                    "calculate_statistic",
                    "statistics",
                    "dcor",
                    "dhsic")) # c("gauss_kerns", "dhsic"))

res_calibrate1 <- calibrate(data_gen_method, n_obs1, n_dim, n_sim, n_perm, n_reps)
res_calibrate2 <- calibrate(data_gen_method, n_obs2, n_dim, n_sim, n_perm, n_reps)
res_calibrate3 <- calibrate(data_gen_method, n_obs3, n_dim, n_sim, n_perm, n_reps)
res_calibrate4 <- calibrate(data_gen_method, n_obs4, n_dim, n_sim, n_perm, n_reps)

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

plot_rr_n_obs(res_calibrate1, n_obs1)
plot_rr_n_obs(res_calibrate2, n_obs2)
plot_rr_n_obs(res_calibrate3, n_obs3)
plot_rr_n_obs(res_calibrate4, n_obs4)




