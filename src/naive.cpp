//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// Naive (sample) autocovariance estimators computed from an observed n1 x n2
// data matrix: the scalar Chat(h), the Gamma_hat matrix of sample covariances
// over all lags, and the N x N field covariance built from it.
//
// Convention: h1 = vertical lag (rows), h2 = horizontal lag (columns).
//-----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "naive.h"

using namespace Rcpp;

//-----------------------------------------------------------------------------
// NAIVE SAMPLE COVARIANCE ESTIMATOR
//-----------------------------------------------------------------------------
// [[Rcpp::export]]
double naive_sample_cov(NumericMatrix X, NumericVector hvec,
                        bool Nh_normalize = true) {
// Naive sample covariance estimator: 
// Chat(h) = (1/N_h) * sum_{t} (X(t+h) - Xbar) * (X(t) - Xbar)
// Parameters:
//   X: n1 x n2 data matrix.
//   hvec: lag vector (h1, h2).
//   Nh_normalize: if true, normalise by N_h; if false, normalise by N = n1 * n2.

  if (hvec.size() != 2)
    stop("hvec must have length 2");
 
  int h1 = (int)hvec[0];
  int h2 = (int)hvec[1];
  int n1 = X.nrow();
  int n2 = X.ncol();
  double N = static_cast<double>(n1 * n2);

  // Computes the sample mean for centering
  double X_mean = 0.0;
  for (int i = 0; i < n1; ++i)
    for (int j = 0; j < n2; ++j)
      X_mean += X(i, j);
  X_mean /= N;


  // Summation bounds
  int T11 = std::max(1, 1 - h1);
  int T12 = std::min(n1, n1 - h1);
  int T21 = std::max(1, 1 - h2);
  int T22 = std::min(n2, n2 - h2);

  double Nh = static_cast<double>((T12 - T11 + 1) * (T22 - T21 + 1));

  double cov_value = 0.0;
  for (int i = T11; i <= T12; ++i) {
    for (int j = T21; j <= T22; ++j) {
      double x1 = X(i - 1, j - 1) - X_mean;
      double x2 = X(i + h1 - 1, j + h2 - 1) - X_mean;
      cov_value += x1 * x2;
    }
  }

  return cov_value / (Nh_normalize ? Nh : N);
}
//-----------------------------------------------------------------------------
// FULL EMPIRICAL GAMMA (n1-1) x (n2-1) MATRIX OF SAMPLE COVARIANCES
//-----------------------------------------------------------------------------

// Build the L1 x L2 matrix Gamma_hat(h1, h2) for all lags from data X.
// [[Rcpp::export]]
NumericMatrix build_Gamma_naive(NumericMatrix X,
                                bool Nh_normalize = true) {

  int n1 = X.nrow();
  int n2 = X.ncol();
  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;

  NumericMatrix Gamma_hat(L1, L2);

  for (int i = 0; i < L1; ++i) {
    int lag1 = i - (n1 - 1);
    for (int j = 0; j < L2; ++j) {
      int lag2 = j - (n2 - 1);
      NumericVector hvec = NumericVector::create(lag1, lag2);
      Gamma_hat(i, j) = naive_sample_cov(X, hvec, Nh_normalize);
    }
  }

  CharacterVector rownames(L1), colnames(L2);
  for (int i = 0; i < L1; ++i)
    rownames[i] = "h1=" + std::to_string(i - (n1 - 1));
  for (int j = 0; j < L2; ++j)
    colnames[j] = "h2=" + std::to_string(j - (n2 - 1));
  Gamma_hat.attr("dimnames") = List::create(rownames, colnames);

  return Gamma_hat;
}
//-----------------------------------------------------------------------------
// NAIVE SIGMA_HAT: N x N covariance matrix from the empirical lag matrix.
// Build the N x N variance covariance matrix Sigma_hat

// [[Rcpp::export]]
NumericMatrix build_Sigma_naive(NumericMatrix X, bool Nh_normalize = true) {

  int n1 = X.nrow();
  int n2 = X.ncol();
  int N  = n1 * n2;

  NumericMatrix Gamma_hat = build_Gamma_naive(X, Nh_normalize);
  NumericMatrix Sigma_hat(N, N);

  for (int s2 = 0; s2 < n2; ++s2) {
    for (int s1 = 0; s1 < n1; ++s1) {
      int s_idx = s2 * n1 + s1;
      for (int t2 = 0; t2 < n2; ++t2) {
        for (int t1 = 0; t1 < n1; ++t1) {
          int t_idx = t2 * n1 + t1;
          int row_g = (s1 - t1) + (n1 - 1);   // 0-based row in Gamma_hat
          int col_g = (s2 - t2) + (n2 - 1);   // 0-based col in Gamma_hat
          Sigma_hat(s_idx, t_idx) = Gamma_hat(row_g, col_g);
        }
      }
    }
  }

  return Sigma_hat;
}
