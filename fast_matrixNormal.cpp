#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sum_of_squares(NumericMatrix mat) {
  double total = 0;
  int n = mat.size();
  for (int i = 0; i < n; i++) {
    total += mat[i] * mat[i];
  }
  return total;
}