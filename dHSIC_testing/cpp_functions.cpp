#include <Rcpp.h>
using namespace Rcpp;

// Function to compute the Gaussian kernel matrix
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

// Function to compute the linear kernel matrix
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

// Function to compute dHSIC
// [[Rcpp::export]]
double cpp_compute_dHSIC(List data, List hyperparameters) {
  int n = as<int>(hyperparameters["n"]);
  int m = data.size();
  
  std::vector<NumericMatrix> K_list(m);
  
  for (int i = 0; i < m; ++i) {
    NumericMatrix data_i = data[i];
    std::string kernel = as<std::string>(hyperparameters["kernel"]);
    if (kernel == "linear_cpp") {
      K_list[i] = cpp_linearKern(data_i);
    } else if (kernel == "gaussian_cpp") {
      double sigma = as<double>(hyperparameters["sigma"]);
      K_list[i] = cpp_gaussKern(data_i, sigma);
    } else {
      stop("Unknown kernel type");
    }
  }
  
  // Centering matrix
  NumericMatrix H(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      H(i, j) = -1.0 / n;
    }
    H(i, i) += 1;
  }
  
  // Center the kernel matrices
  std::vector<NumericMatrix> K_centered(m);
  for (int i = 0; i < m; ++i) {
    NumericMatrix K = K_list[i];
    NumericMatrix HK = no_init_matrix(n, n);
    NumericMatrix HKH = no_init_matrix(n, n);
    
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        double sum_hkj = 0;
        for (int l = 0; l < n; ++l) {
          sum_hkj += H(k, l) * K(l, j);
        }
        HK(k, j) = sum_hkj;
      }
    }
    
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        double sum_hik = 0;
        for (int l = 0; l < n; ++l) {
          sum_hik += H(j, l) * HK(l, k);
        }
        HKH(j, k) = sum_hik;
      }
    }
    
    K_centered[i] = HKH;
  }
  
  // Compute the dHSIC value
  double HSIC_value = 0;
  for (int i = 0; i < m - 1; ++i) {
    for (int j = i + 1; j < m; ++j) {
      for (int k = 0; k < n; ++k) {
        for (int l = 0; l < n; ++l) {
          HSIC_value += K_centered[i](k, l) * K_centered[j](k, l);
        }
      }
    }
  }
  HSIC_value /= pow(n - 1, 2);
  
  return HSIC_value;
}