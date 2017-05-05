// Compute the distance matrix in neighborhood smoothing method

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat compute_dist_nbs(arma::mat A) {
  int n = A.n_rows;
  rowvec temp(n);
  mat d = zeros<mat>(n, n);
  for (int i = 0; i < n - 1; i ++) {
    for (int j = i + 1; j < n; j ++) {
      temp = A.row(i) - A.row(j);
      int k = 0;
      double cur_max = -1;
      while (k < n) {
        if (k != i && k != j) {
          cur_max = std::max(cur_max, std::abs(accu(temp % A.row(k))) / n);
        }
        k ++;
      }
      d(i, j) = cur_max;
    }
  }
  d = sqrt(d + d.t());
  return d;
}