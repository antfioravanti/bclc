//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// Empirical (sample) estimators computed from an observed n1 x n2 data matrix.
// Includes the standard autocovariance,
// the empirical Gamma matrix of sample covariances for all lags,
// Cressie-Hawkins robust covariance estimator,
// and the weight/selection matrices for the bias-corrected estimator.
//
// Convention: h1 = vertical lag (rows), h2 = horizontal lag (columns).
//-----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "counting.h"

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

//-----------------------------------------------------------------------------
// LAG VECTOR TO LAG MATRIX CONVERSION
//-----------------------------------------------------------------------------

static NumericMatrix lagvec_to_matrix(NumericVector Chat_star_m,
                                       int m1, int m2) {

  int M1 = 2 * m1 + 1;   // number of h1 (dim-1) lags
  int M2 = 2 * m2 + 1;   // number of h2 (dim-2) lags

  if (Chat_star_m.size() != M1 * M2)
    stop("Chat_star_m must have length (2*m1+1)*(2*m2+1)");

  NumericMatrix Gamma_hat_corr_m(M1, M2);
  for (int j = 0; j < M2; ++j)
    for (int i = 0; i < M1; ++i)
      Gamma_hat_corr_m(i, j) = Chat_star_m[j * M1 + i];

  CharacterVector rown(M1), coln(M2);
  for (int i = 0; i < M1; ++i)
    rown[i] = "h1=" + std::to_string(i - m1);
  for (int j = 0; j < M2; ++j)
    coln[j] = "h2=" + std::to_string(j - m2);
  Gamma_hat_corr_m.attr("dimnames") = List::create(rown, coln);

  return Gamma_hat_corr_m;
}

// Pick the truncated lag vector to reshape out of a bias_corrected_estimator()
// fit: which = "corrected" -> Chat_star_m (bias-corrected); which = "naive" ->
// Chat_m (tapered/truncated sample vector, no correction). Both have the same
// (2m1+1)*(2m2+1) column-major layout (h1 fastest).
static NumericVector pick_lagvec(List fit, std::string which) {
  if (which == "corrected") return fit["Chat_star_m"];
  if (which == "naive")     return fit["Chat_m"];
  stop("which must be \"corrected\" or \"naive\"");
}

//-----------------------------------------------------------------------------
// NEAREST PSD PROJECTION
//-----------------------------------------------------------------------------

// Project a symmetric matrix to the nearest positive semi-definite matrix
// (Frobenius-optimal). Symmetrises M, floors all negative eigenvalues to zero
// and reconstructs. Parameters:
//   M: numeric symmetric matrix.
// Returns the nearest PSD matrix.
// [[Rcpp::export]]
arma::mat psd_project(arma::mat M) {

  arma::mat M_sym = 0.5 * (M + M.t());     // symmetrise numerically

  arma::vec lam;
  arma::mat evec;
  arma::eig_sym(lam, evec, M_sym);         // ascending eigenvalues, orthonormal

  lam.transform([](double v) { return v > 0.0 ? v : 0.0; });  // floor at 0

  return evec * arma::diagmat(lam) * evec.t();
}


//-----------------------------------------------------------------------------
// BIAS CORRECTED ESTIMATOR
//-----------------------------------------------------------------------------
// Compute the bias-corrected estimator from a data matrix X, with truncation
// cut-offs m1, m2 for the two dimensions (m2 defaults to m1). Returns a list
// with the corrected lag vector (Chat_star_m), the uncorrected truncated lag
// vector (Chat_m), the full lag vector (Chat), the truncated weight matrix 
// (W_star_m) and the scalars
// m1, m2, n1, n2. 
// [[Rcpp::export]]
List bias_corrected_estimator(NumericMatrix X, int m1, int m2 = -1,
                                  bool Nh_normalize = true) {

  if (m2 < 0) m2 = m1; // default: same cut-off in both dims
  int n1 = X.nrow();
  int n2 = X.ncol();
  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;
  int M1 = 2 * m1 + 1;
  int M2 = 2 * m2 + 1;
  // Count matrices
  List mat1 = build_count_matrices(n1);
  List mat2 = build_count_matrices(n2);

  // W^*_m
  List W_mats = build_weights(mat1, mat2, m1, m2, Nh_normalize);
  arma::mat W_star_m = as<arma::mat>(W_mats["W_star_m"]);

  // Gamma_hat (full lag grid) and its column-major vectorisation
  NumericMatrix Gamma_hat = build_Gamma_naive(X, Nh_normalize);
  arma::vec Chat = arma::vectorise(as<arma::mat>(Gamma_hat));

  // Sm selects the central (M1*M2) block; h1 (dim-1) varies fastest
  arma::mat Sm = build_Sm(L1, L2, m1, m2);   // (M1*M2) x (L1*L2)
  arma::vec Chat_m = Sm * Chat;

  // Solve W*_m x = Chat_m  (bias-corrected lag vector)
  arma::vec Chat_star_m = arma::solve(W_star_m, Chat_m);

  // Lag names, column-major over {-m1..m1} x {-m2..m2} with h1 varying fastest.
  CharacterVector vec_names(M1 * M2);
  int k = 0;
  for (int h2 = -m2; h2 <= m2; ++h2)
    for (int h1 = -m1; h1 <= m1; ++h1)
      vec_names[k++] = "C(" + std::to_string(h1) + "," + std::to_string(h2) + ")";

  NumericVector Chat_m_named    = Rcpp::wrap(Chat_m);
  NumericVector Chat_star_m_named = Rcpp::wrap(Chat_star_m);
  Chat_m_named.attr("names")      = vec_names;
  Chat_star_m_named.attr("names") = vec_names;

  return List::create(
    Named("Chat_star_m")        = Chat_star_m_named,
    Named("Chat_m")             = Chat_m_named,
    Named("Chat")               = Rcpp::wrap(Chat),
    Named("W_star_m")           = Rcpp::wrap(W_star_m),
    Named("m1")                 = m1,
    Named("m2")                 = m2,
    Named("n1")                 = n1,
    Named("n2")                 = n2
  );
}

//-----------------------------------------------------------------------------
// LAG MATRIX from a bias_corrected_estimator() FIT
//-----------------------------------------------------------------------------
//
// Reshape a truncated lag vector of a bias_corrected_estimator() fit into the
// M1 x M2 lag matrix (rows = h1, cols = h2).
// Parameters:
//   fit:   output list of bias_corrected_estimator().
//   which: "corrected" (default) -> Chat_star_m (bias-corrected);
//          "naive"               -> Chat_m (tapered sample vector, no correction).
// [[Rcpp::export]]
NumericMatrix fit_to_Gamma(List fit, std::string which = "corrected") {

  NumericVector Chat = pick_lagvec(fit, which);
  int m1 = as<int>(fit["m1"]);
  int m2 = as<int>(fit["m2"]);

  return lagvec_to_matrix(Chat, m1, m2);
}

//-----------------------------------------------------------------------------
// SIGMA MATRIX from a bias_corrected_estimator() FIT
//-----------------------------------------------------------------------------
//
// Builds the N x N plug-in Sigma hat covariance from a bias_corrected_estimator()
// fit result. Reads m1, m2, n1, n2 and the chosen lag vector from the fit list;
// reshapes it (via lagvec_to_matrix) into the (2m1+1) x (2m2+1) lag matrix and
// fills the N x N covariance from it.
// Returns:
//   Sigma_hat_m_O     : N x N field covariance where
//                         Sigma_hat_m_O[s,t] = C(s1-t1, s2-t2),
//                       spatial locations in column-major order, 0 outside the
//                       {-m1..m1} x {-m2..m2} window.
//   Sigma_hat_m_O_psd : nearest-PSD projection of Sigma_hat_m_O (a valid
//                         covariance matrix).
// Parameters:
//   fit:   output list of bias_corrected_estimator().
//   which: "corrected" (default) -> Chat_star_m (bias-corrected);
//          "naive"               -> Chat_m (tapered sample vector, no correction).
// [[Rcpp::export]]
List fit_to_Sigma(List fit, std::string which = "corrected") {

  NumericVector Chat = pick_lagvec(fit, which);
  int m1 = as<int>(fit["m1"]);
  int m2 = as<int>(fit["m2"]);
  int n1 = as<int>(fit["n1"]);
  int n2 = as<int>(fit["n2"]);

  // Lag matrix (rows = h1, cols = h2), shared reshape helper.
  NumericMatrix Gamma_hat_corr_m = lagvec_to_matrix(Chat, m1, m2);

  // N x N field covariance, column-major spatial indexing over the n1 x n2
  // lattice; 0 outside the {-m1..m1} x {-m2..m2} lag window.
  int N = n1 * n2;
  NumericMatrix Sigma_hat_m_O(N, N);

  for (int s2 = 0; s2 < n2; ++s2) {
    for (int s1 = 0; s1 < n1; ++s1) {
      int s_idx = s2 * n1 + s1;            // column-major, 0-based
      for (int t2 = 0; t2 < n2; ++t2) {
        for (int t1 = 0; t1 < n1; ++t1) {
          int h1 = s1 - t1;
          int h2 = s2 - t2;
          if (std::abs(h1) <= m1 && std::abs(h2) <= m2) {
            int t_idx = t2 * n1 + t1;      // column-major, 0-based
            int row_g = h1 + m1;           // row in Gamma_hat_corr_m
            int col_g = h2 + m2;           // col in Gamma_hat_corr_m
            Sigma_hat_m_O(s_idx, t_idx) = Gamma_hat_corr_m(row_g, col_g);
          }
        }
      }
    }
  }

  NumericMatrix Sigma_hat_m_O_psd =
    as<NumericMatrix>(Rcpp::wrap(psd_project(as<arma::mat>(Sigma_hat_m_O))));

  return List::create(
    Named("Sigma_hat_m_O")     = Sigma_hat_m_O,
    Named("Sigma_hat_m_O_psd") = Sigma_hat_m_O_psd
  );
}

//-----------------------------------------------------------------------------
// CRESSIE HAWKINS SEMIVARIOGRAM ESTIMATOR
//-----------------------------------------------------------------------------

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


