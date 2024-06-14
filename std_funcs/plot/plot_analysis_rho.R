plot_analysis_rho <- function(results, data_gen_method, dim) {
  ggplot(results, aes(x = rho, y = mean_rejection_rate, color = statistic, fill = statistic)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_rejection_rate - se_rejection_rate, ymax = mean_rejection_rate + se_rejection_rate), alpha = 0.2) +
    labs(title = paste("Rejection Rates by ρ with SE Bands, Data Generated\nby", 
                       str_to_title(data_gen_method), "Method, Dimension =", dim), 
         x = "Level of Dependence (ρ)", y = "Rejection Rate") +
    theme_minimal()
}
