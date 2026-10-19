//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// Cressie-Hawkins robust semivariogram estimator and the covariance derived
// from it, C_CH(h) = Chat(0,0) - gamma_CH(h).
//
// Convention: h1 = vertical lag (rows), h2 = horizontal lag (columns).
//-----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "naive.h"

using namespace Rcpp;

//-----------------------------------------------------------------------------
// CRESSIE HAWKINS SEMIVARIOGRAM ESTIMATOR
//-----------------------------------------------------------------------------

// [[Rcpp::export]]
double cressiehawkins_semivariogram(NumericMatrix X, NumericVector hvec,
                                    bool Nh_normalize = true) {
  // Model: 
  // sv_CH(h) = [mean(|X(t+h) - X(t)|^{1/2})]^4 / (2 * (0.457 + 0.494/N_h)).                                  

  if (hvec.size() != 2)
    stop("hvec must have length 2");

  int h1 = (int)hvec[0];
  int h2 = (int)hvec[1];
  int n1 = X.nrow();
  int n2 = X.ncol();
  double N = static_cast<double>(n1 * n2);

  int T11 = std::max(1, 1 - h1);
  int T12 = std::min(n1, n1 - h1);
  int T21 = std::max(1, 1 - h2);
  int T22 = std::min(n2, n2 - h2);

  if (T12 < T11 || T22 < T21) return NA_REAL;

  double Nh = static_cast<double>((T12 - T11 + 1) * (T22 - T21 + 1));

  double sum_fourth_root = 0.0;
  for (int i = T11; i <= T12; ++i) {
    for (int j = T21; j <= T22; ++j) {
      double abs_diff = std::abs(X(i + h1 - 1, j + h2 - 1) - X(i - 1, j - 1));
      sum_fourth_root += std::sqrt(abs_diff);
    }
  }

  double denom = Nh_normalize ? Nh : N;
  double mean_fr = sum_fourth_root / denom;
  double fourth_power = mean_fr * mean_fr * mean_fr * mean_fr;

  double correction = 0.457 + 0.494 / denom;

  return 0.5 * fourth_power / correction;
}

//-----------------------------------------------------------------------------
// Cressie-Hawkins robust covariance via C(h) = C(0) - gamma(h)
//-----------------------------------------------------------------------------

// Compute C_CH(h) = Chat(0,0) - gamma_CH(h).
// [[Rcpp::export]]
double cressiehawkins_cov(NumericMatrix X, NumericVector hvec,
                          bool Nh_normalize = true) {

  if (hvec.size() != 2)
    stop("hvec must have length 2");

  double sv_h = cressiehawkins_semivariogram(X, hvec, Nh_normalize);
  if (ISNAN(sv_h)) return NA_REAL;

  NumericVector zero_lag = NumericVector::create(0, 0);
  double C_0 = naive_sample_cov(X, zero_lag, Nh_normalize);

  return C_0 - sv_h;
}
