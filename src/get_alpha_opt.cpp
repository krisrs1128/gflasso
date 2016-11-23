// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat get_alpha_opt(arma::mat C, arma::mat W, double mu) {
  arma::mat alpha = C * W.t();
  for (int i = 0; i < alpha.n_rows; i++) {
    for (int j = 0; j < alpha.n_cols; j++) {
      if (alpha(i, j) > 1) {
	alpha(i, j) = 1;
      }
      if (alpha(i, j) < -1) {
	alpha(i, j) = -1;
      }
    }
  }
  return alpha;
}
