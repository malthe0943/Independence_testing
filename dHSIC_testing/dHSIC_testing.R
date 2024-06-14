library(microbenchmark)
library(Rcpp)
library(ggplot2)

source("./dHSIC.R")
source("./misc/custom_kernel_functions.R")
sourceCpp("./dHSIC_testing/cpp_functions.cpp")

x <- matrix(rnorm(200 * 1), ncol = 1)
y <- matrix(rnorm(200), ncol = 1)

# Test result of Gaussian kernel with sigma = 0.3 with own implementation,
# the dHSIC package and the C++ implementation
dhsic(list(x, y), kernel = c("gauss_kern_3"))$dHSIC

compute_dHSIC(list(x, y), list(kernel = "gaussian", sigma = 0.3))

compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = 0.3))

cpp_compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = 0.3, n = nrow(x)))

# Test result of Gaussian kernel with sigma = median with own implementation,
# the dHSIC package and the C++ implementation
# NOT IMPORTANT! dHSIC is fine
dhsic(list(x, y), kernel = c("gaussian"))$dHSIC

# compute_dHSIC(list(x, y), list(kernel = "gaussian_median", sigma = NULL))

# compute_dHSIC(list(x, y), list(kernel = "gaussian_median_cpp", sigma = NULL))

# cpp_compute_dHSIC(list(x, y), list(kernel = "gaussian_median_cpp", sigma = NULL, n = nrow(x)))

# Test result of linear kernel with own implementation, the dHSIC package
# and the C++ implementation
dhsic(list(x, y), kernel = c("linear_kernel"))$dHSIC

compute_dHSIC(list(x, y), list(kernel = "linear"))

compute_dHSIC(list(x, y), list(kernel = "linear_cpp"))

cpp_compute_dHSIC(list(x, y), list(kernel = "linear_cpp", sigma = NULL, n = nrow(x)))

# Profiling our own dHSIC implementations to see where to optimise
profvis::profvis({
  compute_dHSIC(list(x, y), list(kernel = "linear", sigma = 0.3))
})

benchmark_results <- microbenchmark(
  "dHSIC Custom Gaussian Kernel" = dhsic(list(x, y), kernel = c("gauss_kern_3"))$dHSIC,
  "My dHSIC R Gaussian" = mydHSIC::compute_dHSIC(list(x, y), list(kernel = "gaussian", sigma = 0.3)),
  "My dHSIC R C++ Gaussian" = mydHSIC::compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = 0.3)),
  "My dHSIC Only C++ Gaussian" = cpp_compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = 0.3, n = nrow(x))),
  "dHSIC Default Gaussian Kernel" = dhsic(list(x, y), kernel = c("gaussian"))$dHSIC,
  "dHSIC Custom Linear Kernel" = dhsic(list(x, y), kernel = c("linear_kernel"))$dHSIC,
  "My dHSIC R Linear" = mydHSIC::compute_dHSIC(list(x, y), list(kernel = "linear")),
  "My dHSIC R C++ Linear" = mydHSIC::compute_dHSIC(list(x, y), list(kernel = "linear_cpp")),
  "My dHSIC C++ Linear" = cpp_compute_dHSIC(list(x, y), list(kernel = "linear_cpp", sigma = NULL, n = nrow(x))),
  times = 100
)

# Define custom colors
color_mapping <- c(
  "dHSIC Custom Gaussian Kernel" = "green",
  "My dHSIC R Gaussian" = "green",
  "My dHSIC R C++ Gaussian" = "orange",
  "My dHSIC Only C++ Gaussian" = "red",
  "dHSIC Default Gaussian Kernel" = "green",
  "dHSIC Custom Linear Kernel" = "green",
  "My dHSIC R Linear" = "green",
  "My dHSIC R C++ Linear" = "orange",
  "My dHSIC C++ Linear" = "red"
)

# Convert benchmark results to data frame for manual plotting
benchmark_df <- as.data.frame(benchmark_results)

# Manually plot using ggplot2
ggplot(benchmark_df, aes(x = expr, y = time / 1e6, color = expr, fill = expr)) +  # Convert nanoseconds to milliseconds
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  ggtitle("Performance Comparison of Different dHSIC Computations") +
  theme_bw() +
  xlab("Function") +
  ylab("Time (milliseconds)") +
  scale_y_log10() +
  scale_colour_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1.01), legend.position = "none")

# From this use for
# dHSIC linear kernel: My dHSIC R Linear
# dHSIC gaussian kernel: dHSIC R Gaussian default
# dHSIC optimised sigma kernel: My dHSIC R with C++ kernel

# PUFH: Continue with the optimisation of sigma
# (small scale and prepare for overnight run)
# Then do the calibration
# (small scale and prepare for overnight run)
# Then do the H1 case
