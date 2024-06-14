calculate_statistic <- function( x, 
                                 y, 
                                 statistic,
                                 sigma) {
  if (statistic == "dCor") {
    return(dcor(x, y))
  } 
  else if (statistic == "linear") {
    return(mydHSIC::compute_dHSIC(list(x, y), list(kernel = "linear_cpp")))
  } 
  else if (statistic == "gaussian") {
    return(dhsic(list(x, y), kernel = c("gaussian"))$dHSIC)
  } 
  else if (statistic == "gaussian_opt") {
    return(mydHSIC::compute_dHSIC(list(x, y), list(kernel = "gaussian_cpp", sigma = sigma)))
  } 
  else {
    stop("Unknown statistic")
  }
}

