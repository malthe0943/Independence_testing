# Load necessary libraries
library(ggplot2)

# Function to simulate data and plot with correlation coefficients
simulate_and_plot <- function(x, y, title_prefix) {
  # Calculate Pearson's correlation
  pearson_r <- cor(x, y, method = "pearson")

  # Calculate Spearman's correlation
  spearman_r <- cor(x, y, method = "spearman")

  # Create data frame
  data <- data.frame(x = x, y = y)

  # Plot the data with correlation coefficients in the title
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    ggtitle(paste(title_prefix, "\nPearson's ρ =", round(pearson_r, 2),
                  ", Spearman's ρ =", round(spearman_r, 2))) +
    theme_bw()

  print(p)
}

# Set seed for reproducibility
set.seed(42)

# Generate data for different frameworks

# Independent data
x_independent <- rnorm(100)
y_independent <- rnorm(100)
simulate_and_plot(x_independent, y_independent, "Independent")

# Linear dependence
x_linear <- rnorm(100)
y_linear <- 2 * x_linear + rnorm(100)
simulate_and_plot(x_linear, y_linear, "Linear Dependence")

# Exponential dependence
x_exponential <- rnorm(100)
y_exponential <- exp(x_exponential) + rnorm(100, mean = 0, sd = 0.5)
simulate_and_plot(x_exponential, y_exponential, "Exponential Dependence")

# U-shaped dependence
x_u_shaped <- rnorm(100)
y_u_shaped <- (x_u_shaped)^2 + rnorm(100)
simulate_and_plot(x_u_shaped, y_u_shaped, "U-shaped Dependence")
