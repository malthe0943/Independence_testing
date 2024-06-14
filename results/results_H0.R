library(miceadds)
library(mydHSIC)
library(Rcpp)

set.seed(42)

source.all("./std_funcs/plot", ".R")
source.all("./std_funcs/calc", ".R")

data_gen_method <- "independent"
n_obs <- 200 # 200
n_dim <- 1
n_sim <- 100 # 100
n_perm <- 100 # 100
n_reps <- 5 # 5

# Register parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("simulate_data", 
                    "calculate_statistic",
                    "statistics",
                    "dcor",
                    "dhsic")) # c("gauss_kerns", "dhsic"))

res_calibrate <- calibrate(data_gen_method, n_obs, n_dim, n_sim, n_perm, n_reps)

stopCluster(cl)

Sys.time()

plot_rr(res_calibrate)




# data_gen_method <- "linear"
# data_gen_method <- "normal"
# data_gen_method <- "M1"
# data_gen_method <- "M2"
# data_gen_method <- "M3"
