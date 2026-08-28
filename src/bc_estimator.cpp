//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// The bias-corrected lattice covariance estimator: C*_m = (W*_m)^{-1} Chat_m,
// and helpers to turn a fit into the lag matrix Gamma or Sigma
//
// Convention: h1 = vertical lag (rows), h2 = horizontal lag (columns).
//-----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "counting.h"
#include "naive.h"
#include "select_m.h"
#include "bc_estimator.h"
using namespace Rcpp;

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

// Pick which estimator's truncated lag vector to reshape out of a
// bias_corrected_estimator() fit:
//  "corrected"  Chat_star_m (bias-corrected); 
//  "truncated" Chat_m (tapered/truncated sample vector, no correction). 
//Both have the same (2m1+1)*(2m2+1) column-major layout.
static NumericVector pick_estimator(List fit, std::string which) {
  if (which == "corrected") return fit["Chat_star_m"];
  if (which == "truncated") return fit["Chat_m"];
  stop("which must be \"corrected\" or \"truncated\"");
}

//-----------------------------------------------------------------------------
// NEAREST PSD PROJECTION
//-----------------------------------------------------------------------------

// Project a symmetric matrix to the nearest positive semi-definite matrix
// by flooring all negative eigenvalues to zero
// Parameters:
//   M: numeric symmetric matrix.
// Returns the nearest PSD matrix.
// [[Rcpp::export]]
arma::mat psd_project(arma::mat M) {
  arma::mat M_sym = 0.5 * (M + M.t());
  arma::vec lam;
  arma::mat evec;
  arma::eig_sym(lam, evec, M_sym);
  // floor at 0
  lam.transform([](double eigenv) { return eigenv > 0.0 ? eigenv : 0.0; });  
  return evec * arma::diagmat(lam) * evec.t();
}

//-----------------------------------------------------------------------------
// BIAS CORRECTED ESTIMATOR
//-----------------------------------------------------------------------------
// Compute the bias-corrected estimator C*_m = (W*_m)^{-1} Chat_m.
// Parameters:
//   X:            n1 x n2 data matrix.
//   m1, m2:       truncation cut-offs per dimension (m2 defaults to m1).
//   Nh_normalize: N_h (TRUE) or constant 1/N (FALSE) normalisation.
// Returns:
//   Chat_star_m    : corrected truncated lag vector.
//   Chat_m         : uncorrected truncated lag vector.
//   Chat           : full lag vector.
//   W_star_m       : truncated weight matrix.
//   m1, m2, n1, n2 : cut-offs and grid dimensions.
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
// fit_lattice_cov() = select_m() + bias_corrected_estimator()
//-----------------------------------------------------------------------------
// One-call automatic fit: pick the cut-off with select_m() and run
// bias_corrected_estimator() at it. The cut-off is ALWAYS data-driven -- there
// is deliberately no manual m argument here; to set m1/m2 by hand, call
// bias_corrected_estimator() directly.
// Parameters:
//   X:            n1 x n2 data matrix.
//   method:       forwarded to select_m() ("rule" or "cv").
//   pooled:       forwarded to select_m() (single pooled cut-off for both axes;
//                 with method = "cv" this is also what makes the search cheap).
//   kappa, eta:   forwarded to select_m() (tuning constants of the rule).
//   K:            forwarded to select_m() (correlogram lags used in the WLS fit).
//   cv_grid, cv_block, cv_shift, cv_nugget: forwarded to select_m(), used only
//                 when method = "cv".
//   Nh_normalize: forwarded to both select_m() and bias_corrected_estimator().
// Returns the bias_corrected_estimator() fit with the select_m() list appended
// as m_selection, recording how the cut-off was reached.
// [[Rcpp::export]]
List fit_lattice_cov(NumericMatrix X,
                     std::string method = "rule",
                     bool pooled = false,
                     double kappa = 0.445,
                     double eta = 0.5,
                     int K = -1,
                     bool Nh_normalize = true,
                     IntegerVector cv_grid = IntegerVector::create(),
                     int cv_block = -1,
                     int cv_shift = 2,
                     double cv_nugget = 1e-6) {

  // Cut-off is always chosen by the selection rule.
  List sel = select_m(X, method, pooled, kappa, eta, K, -1, Nh_normalize,
                      cv_grid, cv_block, cv_shift, cv_nugget);
  int  m1  = as<int>(sel["m1"]);
  int  m2  = as<int>(sel["m2"]);

  // Copy the fit, then append the selection record.
  List fit = bias_corrected_estimator(X, m1, m2, Nh_normalize);
  fit["m_selection"] = sel;
  return fit;
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
//          "truncated"           -> Chat_m (tapered sample vector, no correction).
// Returns:
//   Gamma_hat_corr_m : (2m1+1) x (2m2+1) lag matrix, dimnames "h1=..." / "h2=...".
// [[Rcpp::export]]
NumericMatrix fit_to_Gamma(List fit, std::string which = "corrected") {

  NumericVector Chat = pick_estimator(fit, which);
  int m1 = as<int>(fit["m1"]);
  int m2 = as<int>(fit["m2"]);

  return lagvec_to_matrix(Chat, m1, m2);
}

//-----------------------------------------------------------------------------
// SIGMA MATRIX from a bias_corrected_estimator() FIT
//-----------------------------------------------------------------------------
//
// Builds the plug-in Sigma hat covariance from a bias_corrected_estimator()
// fit.
// Parameters:
//   fit:   output list of bias_corrected_estimator().
//   which: "corrected" (default) -> Chat_star_m (bias-corrected);
//          "truncated"           -> Chat_m (tapered sample vector).
//   psd:   TRUE (default) -> the PSD projection, a valid covariance matrix;
//          FALSE -> the raw plug-in, which may be indefinite.
// Returns:
//   Sigma_hat_m : field covariance where Sigma_hat_m[s,t] = C(s1-t1, s2-t2),
//                 spatial locations in column-major order, 0 outside the
//                 {-m1..m1} x {-m2..m2} window;
//   psd_applied : the psd argument reported back, so a caller holding only the
//                 result can tell which of the two it has.

// [[Rcpp::export]]
List fit_to_Sigma(List fit, std::string which = "corrected", bool psd = true) {

  NumericVector Chat = pick_estimator(fit, which);
  int m1 = as<int>(fit["m1"]);
  int m2 = as<int>(fit["m2"]);
  int n1 = as<int>(fit["n1"]);
  int n2 = as<int>(fit["n2"]);

  // Lag matrix (rows = h1, cols = h2), shared reshape helper.
  NumericMatrix Gamma_hat_corr_m = lagvec_to_matrix(Chat, m1, m2);

  // N x N field covariance, column-major spatial indexing over the n1 x n2
  // lattice; 0 outside the {-m1..m1} x {-m2..m2} lag window.
  int N = n1 * n2;
  NumericMatrix Sigma_hat_m(N, N);

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
            Sigma_hat_m(s_idx, t_idx) = Gamma_hat_corr_m(row_g, col_g);
          }
        }
      }
    }
  }
  if (psd)
    Sigma_hat_m =
      as<NumericMatrix>(Rcpp::wrap(psd_project(as<arma::mat>(Sigma_hat_m))));

  return List::create(
    Named("Sigma_hat_m") = Sigma_hat_m,
    Named("psd_applied") = psd
  );
}
