#include <Rcpp.h>
using namespace Rcpp;

//' Function to compute the Gaussian kernel matrix
//' 
//' @param data A numeric matrix of data points
//' @param sigma A numeric value for the bandwidth of the kernel
//' @return A numeric matrix of the kernel matrix
// [[Rcpp::export]]
NumericMatrix cpp_gaussKern(NumericMatrix data, double sigma) {
  int n = data.nrow();
  NumericMatrix K(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum_sq = 0;
      for (int k = 0; k < data.ncol(); ++k) {
        sum_sq += pow(data(i, k) - data(j, k), 2);
      }
      K(i, j) = exp(-sum_sq / (2 * sigma * sigma));
    }
  }
  return K;
}

//' Function to compute the linear kernel matrix
//' 
//' @param data A numeric matrix of data points
//' @return A numeric matrix of the kernel matrix
// [[Rcpp::export]]
NumericMatrix cpp_linearKern(NumericMatrix data) {
  int n = data.nrow();
  NumericMatrix K(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double dot_product = 0;
      for (int k = 0; k < data.ncol(); ++k) {
        dot_product += data(i, k) * data(j, k);
      }
      K(i, j) = dot_product;
    }
  }
  return K;
}
