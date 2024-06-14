plot_rr_n_sim <- function(results, n_sim) {
  ggplot(results, aes(x = alpha, y = mean_rejection_rate, color = statistic, fill = statistic)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_rejection_rate - se_rejection_rate, ymax = mean_rejection_rate + se_rejection_rate), alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste("Rejection Rates by Alpha with SE Bands,", n_sim, "simulations" ), x = "Significance Level (alpha)", y = "Rejection Rate") +
    theme_minimal()
}
