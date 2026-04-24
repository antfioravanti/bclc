// ============================================================================
// Author: Antonio Fioravanti
// Empirical (sample) estimators computed from an observed n1 x n2 data matrix.
// Includes the standard autocovariance,
// the empirical Gamma matrix of sample covariances for all lags,
// Cressie-Hawkins robust covariance estimator,
// and the weight/selection matrices for the bias-corrected estimator.
//
// Convention: h1 = vertical lag (rows), h2 = horizontal lag (columns).
// ============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "counting.h"

using namespace Rcpp;


// ----------------------------------------------------------------------------
// NAIVE SAMPLE COVARIANCE ESTIMATOR
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

// Build the L1 x L2 matrix Chat(h1, h2) for all lags from data X.
// [[Rcpp::export]]
NumericMatrix build_Gamma_naive(NumericMatrix X,
                                bool Nh_normalize = true) {

  int n1 = X.nrow();
  int n2 = X.ncol();
  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;

  NumericMatrix Chat(L1, L2);

  for (int i = 0; i < L1; ++i) {
    int lag1 = i - (n1 - 1);
    for (int j = 0; j < L2; ++j) {
      int lag2 = j - (n2 - 1);
      NumericVector hvec = NumericVector::create(lag1, lag2);
      Chat(i, j) = naive_sample_cov(X, hvec, Nh_normalize);
    }
  }

  CharacterVector rownames(L1), colnames(L2);
  for (int i = 0; i < L1; ++i)
    rownames[i] = "h1=" + std::to_string(i - (n1 - 1));
  for (int j = 0; j < L2; ++j)
    colnames[j] = "h2=" + std::to_string(j - (n2 - 1));
  Chat.attr("dimnames") = List::create(rownames, colnames);

  return Chat;
}
// -----------------------------------------------------------------------------
// NAIVE SIGMA_HAT: N x N covariance matrix from the empirical lag matrix.
//
// Computes Gamma_hat internally from X, then fills
// Sigma_hat[s,t] = Gamma_hat(s1-t1, s2-t2) for all pairs of spatial
// locations indexed in column-major order.
// No truncation is applied — all lags are used.
//
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

// -----------------------------------------------------------------------------
// RESHAPE CORRECTED ESTIMATOR VECTOR INTO LAG MATRIX AND COVARIANCE MATRIX
//
// Takes Chat_star_m (vec(Gamma_hat_corr_m), length (2m+1)^2, column-major in
// lag space) and produces:
//   Gamma_hat_corr_m : (2m+1) x (2m+1) lag-indexed matrix (rows=h1, cols=h2)
//   Sigma_hat_m      : n1*n2 x n1*n2 plug-in field covariance matrix where
//                        Sigma_hat_m[s,t] = Chat_star_m(s1-t1, s2-t2)
//                      with spatial locations in column-major order.
// [[Rcpp::export]]
List reshape_Sigma_hat_m(NumericVector Chat_star_m, int m, int n1, int n2) {

  int L_m = 2 * m + 1;

  if (Chat_star_m.size() != L_m * L_m)
    stop("Chat_star_m must have length (2m+1)^2");

  // --- (1) Lag matrix: reshape Chat_star_m column-major into (2m+1) x (2m+1) ---
  NumericMatrix Gamma_hat_corr_m(L_m, L_m);
  for (int j = 0; j < L_m; ++j)
    for (int i = 0; i < L_m; ++i)
      Gamma_hat_corr_m(i, j) = Chat_star_m[j * L_m + i];

  CharacterVector rown(L_m), coln(L_m);
  for (int k = 0; k < L_m; ++k) {
    rown[k] = "h1=" + std::to_string(k - m);
    coln[k] = "h2=" + std::to_string(k - m);
  }
  Gamma_hat_corr_m.attr("dimnames") = List::create(rown, coln);

  // --- (2) Sigma_hat_m: n1*n2 x n1*n2, column-major spatial indexing ---
  int N = n1 * n2;
  NumericMatrix Sigma_hat_m(N, N);

  for (int s2 = 0; s2 < n2; ++s2) {
    for (int s1 = 0; s1 < n1; ++s1) {
      int s_idx = s2 * n1 + s1;            // column-major, 0-based
      for (int t2 = 0; t2 < n2; ++t2) {
        for (int t1 = 0; t1 < n1; ++t1) {
          int h1 = s1 - t1;
          int h2 = s2 - t2;
          if (std::abs(h1) <= m && std::abs(h2) <= m) {
            int t_idx = t2 * n1 + t1;      // column-major, 0-based
            int row_g = h1 + m;            // row in Gamma_hat_corr_m
            int col_g = h2 + m;            // col in Gamma_hat_corr_m
            Sigma_hat_m(s_idx, t_idx) = Gamma_hat_corr_m(row_g, col_g);
          }
        }
      }
    }
  }

  return List::create(
    Named("Gamma_hat_corr_m") = Gamma_hat_corr_m,
    Named("Sigma_hat_m")      = Sigma_hat_m
  );
}

//-----------------------------------------------------------------------------
// CRESSIE HAWKINS SEMIVARIOGRAM ESTIMATOR

// [[Rcpp::export]]
double cressiehawkins_semivariogram(NumericMatrix X, NumericVector hvec,
                                    bool Nh_normalize = true) {
  // Model: 
  // sv_CH(h) = [mean(|X(t+h) - X(t)|^{1/2})]^4 / (2 * (0.457 + 0.494/N_h)).
  //                                  

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

// --- Cressie-Hawkins robust covariance via C(h) = C(0) - gamma(h) -----------

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


// =============================================================================
// WEIGHT / SELECTION MATRICES FOR THE BIAS-CORRECTED ESTIMATOR
//
// Replaces build_Pm(), build_Sm(), build_W_m() from bias_correction.R.
//
// Key optimisation: since P_m is a pure selection matrix (one 1 per row),
//   P %*% A %*% t(P)  ==  A[central_rows, central_rows]
// so we never form P explicitly and avoid all matrix multiplications.
// =============================================================================

// Internal helper: 0-based central indices for lag dimension of size L, half-width m.
static arma::uvec central_idx(int L, int m) {
  int center = (int)std::ceil(L / 2.0) - 1;   // 0-based; equals n_k - 1
  return arma::regspace<arma::uvec>(center - m, center + m);
}

// -----------------------------------------------------------------------------
// build_Pm_cpp
//   Returns the (2m+1) x L selection matrix P_m.
//   Exported mainly for inspection; not needed internally.
// [[Rcpp::export]]
arma::mat build_Pm(int L, int m) {
  int p_size = 2 * m + 1;
  arma::mat P = arma::zeros<arma::mat>(p_size, L);
  arma::uvec idx = central_idx(L, m);
  for (int p = 0; p < p_size; ++p)
    P(p, idx(p)) = 1.0;
  return P;
}

// -----------------------------------------------------------------------------
// build_Sm_cpp
//   Returns the (2m+1)^2 x (L1*L2) selection matrix S_m = kron(P2, P1).
//   In practice, prefer slice_indices() to avoid materialising this matrix.
// [[Rcpp::export]]
arma::mat build_Sm(int L1, int L2, int m) {
  arma::mat P1 = build_Pm(L1, m);
  arma::mat P2 = build_Pm(L2, m);
  return arma::kron(P2, P1);
}

// -----------------------------------------------------------------------------
// build_W_m_cpp
//   Builds W_m and W*_m directly from the count-matrix lists returned by
//   build_count_matrices().  Uses index slicing instead of P %*% A %*% t(P).
// [[Rcpp::export]]
List build_W_m(List mat1, List mat2, int m) {

  // Pull matrices out of the R lists
  arma::mat A1  = as<arma::mat>(mat1["A"]);
  arma::mat B1  = as<arma::mat>(mat1["B"]);
  arma::mat D1  = as<arma::mat>(mat1["D"]);
  arma::mat A1n = as<arma::mat>(mat1["A_norm"]);
  arma::mat B1n = as<arma::mat>(mat1["B_norm"]);
  arma::mat D1n = as<arma::mat>(mat1["D_norm"]);
  int L1        = as<int>(mat1["L"]);

  arma::mat A2  = as<arma::mat>(mat2["A"]);
  arma::mat B2  = as<arma::mat>(mat2["B"]);
  arma::mat D2  = as<arma::mat>(mat2["D"]);
  arma::mat A2n = as<arma::mat>(mat2["A_norm"]);
  arma::mat B2n = as<arma::mat>(mat2["B_norm"]);
  arma::mat D2n = as<arma::mat>(mat2["D_norm"]);
  int L2        = as<int>(mat2["L"]);

  // Central index ranges (replaces P %*% A %*% t(P) with a submatrix slice)
  arma::uvec idx1 = central_idx(L1, m);
  arma::uvec idx2 = central_idx(L2, m);

  int krondim = (2 * m + 1) * (2 * m + 1);
  arma::mat I  = arma::eye<arma::mat>(krondim, krondim);

  // Unnormalised W_m
  arma::mat W_m = I
    - arma::kron(A2.submat(idx2, idx2), A1.submat(idx1, idx1))
    - arma::kron(B2.submat(idx2, idx2), B1.submat(idx1, idx1))
    + arma::kron(D2.submat(idx2, idx2), D1.submat(idx1, idx1));

  // Normalised W*_m
  arma::mat W_star_m = I
    - arma::kron(A2n.submat(idx2, idx2), A1n.submat(idx1, idx1))
    - arma::kron(B2n.submat(idx2, idx2), B1n.submat(idx1, idx1))
    + arma::kron(D2n.submat(idx2, idx2), D1n.submat(idx1, idx1));

  return List::create(
    Named("W_m")      = wrap(W_m),
    Named("W_star_m") = wrap(W_star_m)
  );
}

// -----------------------------------------------------------------------------
// bias_corrected_estimator_cpp
//   Full bias-corrected estimator in C++.  Optimisations over the R version:
//     (1) W_m built via index slicing (no P matrix, no matrix multiplications)
//     (2) Chat_m extracted via submatrix slice of Gamma_hat (no S_m matrix)
//     (3) arma::solve(W_star_m, Chat_m) — no explicit matrix inverse
//
//   mat1, mat2: lists returned by build_count_matrices() called from R.
//   (build_count_matrices lives in counting.cpp — a separate DLL — so it
//    cannot be called directly here; pass the lists in from R instead.)
// [[Rcpp::export]]
List bias_corrected_estimator(NumericMatrix X, int m,
                                  bool Nh_normalize = true) {
  int n1 = X.nrow();
  int n2 = X.ncol();

  // Count matrices — now available via counting.h
  List mat1 = build_count_matrices(n1, Nh_normalize);
  List mat2 = build_count_matrices(n2, Nh_normalize);

  // W*_m via optimised C++ builder
  List W_mats = build_W_m(mat1, mat2, m);
  arma::mat W_star_m = as<arma::mat>(W_mats["W_star_m"]);

  // Gamma_hat (C++) — reuse existing function
  NumericMatrix Gamma_r = build_Gamma_naive(X, Nh_normalize);
  arma::mat Gamma_hat   = as<arma::mat>(Gamma_r);

  // Extract Chat_m: central (2m+1)x(2m+1) submatrix, vectorised column-major
  // rows of Gamma_hat = h1 lags, cols = h2 lags; center is (n1-1, n2-1) 0-based
  arma::uvec ridx = arma::regspace<arma::uvec>(n1 - 1 - m, n1 - 1 + m);
  arma::uvec cidx = arma::regspace<arma::uvec>(n2 - 1 - m, n2 - 1 + m);
  arma::vec Chat_m = arma::vectorise(Gamma_hat.submat(ridx, cidx));

  // Solve W*_m x = Chat_m  (no inverse formed)
  arma::vec Chat_star_m = arma::solve(W_star_m, Chat_m);

  // Lag names: same convention as R version — "C(h1,h2)", column-major order
  int L_m = 2 * m + 1;
  CharacterVector vec_names(L_m * L_m);
  int k = 0;
  for (int j = -m; j <= m; ++j)
    for (int i = -m; i <= m; ++i)
      vec_names[k++] = "C(" + std::to_string(i) + "," + std::to_string(j) + ")";

  NumericVector Chat_m_r    = wrap(Chat_m);
  NumericVector Chat_star_r = wrap(Chat_star_m);
  Chat_m_r.attr("names")    = vec_names;
  Chat_star_r.attr("names") = vec_names;

  return List::create(
    Named("Chat_star_m") = Chat_star_r,
    Named("Chat_m")      = Chat_m_r,
    Named("W_star_m")    = wrap(W_star_m),
    Named("m")           = m,
    Named("n1")          = n1,
    Named("n2")          = n2
  );
}
