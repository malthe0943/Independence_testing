plot_analysis_d <- function(results, data_gen_method) {
  ggplot(results, aes(x = dim, y = mean_rejection_rate, color = statistic, fill = statistic)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_rejection_rate - se_rejection_rate, ymax = mean_rejection_rate + se_rejection_rate), alpha = 0.2) +
    scale_x_continuous(trans='log10', breaks=c(4, 10, 30, 100)) +
    labs(title = paste("Rejection Rates by dimension with SE Bands, Data Generated\nby", 
                       str_to_title(data_gen_method), "Method, ρ is 2/5 of maximum possible value"), 
         x = "Dimension, d", y = "Rejection Rate") +
    theme_minimal()
}
