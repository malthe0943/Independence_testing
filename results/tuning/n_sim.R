library(miceadds)
library(mydHSIC)
library(Rcpp)

set.seed(42)

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

n_sim1 <- 5
n_sim2 <- 15
n_sim3 <- 40
n_sim4 <- 100

data_gen_method <- "independent"
n_obs <- 100 # 200
n_dim <- 1
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

res_calibrate1 <- calibrate(data_gen_method, n_obs, n_dim, n_sim1, n_perm, n_reps)
res_calibrate2 <- calibrate(data_gen_method, n_obs, n_dim, n_sim2, n_perm, n_reps)
res_calibrate3 <- calibrate(data_gen_method, n_obs, n_dim, n_sim3, n_perm, n_reps)
res_calibrate4 <- calibrate(data_gen_method, n_obs, n_dim, n_sim4, n_perm, n_reps)

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

plot_rr_n_sim(res_calibrate1, n_sim1)
plot_rr_n_sim(res_calibrate2, n_sim2)
plot_rr_n_sim(res_calibrate3, n_sim3)
plot_rr_n_sim(res_calibrate4, n_sim4)




