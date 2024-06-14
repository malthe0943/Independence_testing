plot_analysis_a_M3 <- function(results, data_gen_method, dim) {
  ggplot(results, aes(x = rho, y = mean_rejection_rate, color = statistic, fill = statistic)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_rejection_rate - se_rejection_rate, ymax = mean_rejection_rate + se_rejection_rate), alpha = 0.2) +
    scale_x_continuous(breaks = seq(1,10,1)) +
    labs(title = paste("Rejection Rates by a with SE Bands, Data Generated\nby", 
                       str_to_title(data_gen_method), "Method"), 
         x = "Level of Dependence (a)", y = "Rejection Rate") +
    theme_minimal()
}
