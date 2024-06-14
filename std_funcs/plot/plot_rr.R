plot_rr <- function(results) {
  ggplot(results, aes(x = alpha, y = mean_rejection_rate, color = statistic, fill = statistic)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_rejection_rate - se_rejection_rate, ymax = mean_rejection_rate + se_rejection_rate), alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = "Rejection Rates by Alpha with SE Bands", x = "Significance Level (alpha)", y = "Rejection Rate") +
    theme_minimal()
}
