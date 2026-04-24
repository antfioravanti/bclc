// ============================================================================
// Author: Antonio Fioravanti
// Analytical (true-parameter) computation of the covariance structure of
// sample autocovariance estimators on a 2D rectangular lattice.
//
// Core idea: precompute C(d1, d2) and the partial sums
//   S1(a, b) = sum_{t1,t2} C(a - t1, b - t2)
// in lookup tables, then evaluate the 16-term Isserlis expansion of
//   Cov(Chat(h), Chat(g))
// using fast O(1) table lookups inside the quadruple loop.
//
// Provides builders for:
//   - Single entry:   cov_sample_autocov_opt, bias_Chat
//   - Full matrices:  build_Omega_m / build_Omega_full
//   - Bias vectors:   build_bias_m / build_bias_full
//   - MSE = Omega + b b':  build_MSE_m / build_MSE_full, expand_MSE variants
//
// Depends on covariance.h for the true_cov() inline function.
// Uses RcppArmadillo for the precomputed arma::mat lookup tables.
// ============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "covariance.h"

using namespace Rcpp;
using namespace arma;


// ============================================================================
// SECTION 1: Precomputed lookup tables (file-scoped, persist across R calls)
// ============================================================================

// Lookup table for C(d1, d2) over all needed differences.
static mat C_matrix;
static int C_offset1, C_offset2;

// Cached parameters to detect when reinitialisation is needed.
static int cached_n1 = -1, cached_n2 = -1;
static double cached_beta = -999, cached_l1 = -999, cached_l2 = -999;
static double cached_a1 = -999, cached_a2 = -999, cached_sigma = -999;

// Precompute C(d1, d2) for all differences that can arise.
void init_C_matrix(int n1, int n2, double beta,
                   double lambda1, double lambda2,
                   double alpha1, double alpha2,
                   double sigma) {

  int max_d1 = 2 * n1;
  int max_d2 = 2 * n2;
  C_offset1 = max_d1;
  C_offset2 = max_d2;

  C_matrix.set_size(2 * max_d1 + 1, 2 * max_d2 + 1);

  for (int d1 = -max_d1; d1 <= max_d1; d1++) {
    for (int d2 = -max_d2; d2 <= max_d2; d2++) {
      C_matrix(d1 + max_d1, d2 + max_d2) =
        true_cov(d1, d2, beta, lambda1, lambda2, alpha1, alpha2, sigma);
    }
  }

  // Update cache
  cached_n1 = n1; cached_n2 = n2;
  cached_beta = beta; cached_l1 = lambda1; cached_l2 = lambda2;
  cached_a1 = alpha1; cached_a2 = alpha2; cached_sigma = sigma;
}

// O(1) covariance lookup.
inline double C_lookup(int d1, int d2) {
  return C_matrix(d1 + C_offset1, d2 + C_offset2);
}


// Lookup table for S1(a, b) = sum_{t1=1}^{n1} sum_{t2=1}^{n2} C(a-t1, b-t2).
static mat S1_matrix;
static int S1_offset1, S1_offset2;

// Precompute S1 for all anchor points that the quadruple loop can produce.
void init_S1_matrix(int n1, int n2) {

  int min_a = -n1, max_a = 2 * n1;
  int min_b = -n2, max_b = 2 * n2;

  S1_offset1 = -min_a;
  S1_offset2 = -min_b;

  S1_matrix.set_size(max_a - min_a + 1, max_b - min_b + 1);

  for (int a = min_a; a <= max_a; a++) {
    for (int b = min_b; b <= max_b; b++) {
      double s = 0.0;
      for (int t1 = 1; t1 <= n1; t1++)
        for (int t2 = 1; t2 <= n2; t2++)
          s += C_lookup(a - t1, b - t2);
      S1_matrix(a + S1_offset1, b + S1_offset2) = s;
    }
  }
}

// O(1) partial-sum lookup.
inline double S1_lookup(int a, int b) {
  return S1_matrix(a + S1_offset1, b + S1_offset2);
}


// ============================================================================
// SECTION 2: Derived quantities
// ============================================================================

// E = (1/N^2) sum_{d1,d2} (n1-|d1|)(n2-|d2|) C(d1,d2).
// This is the expectation of the squared sample mean (calE / N^2).
double compute_E(int n1, int n2) {
  double E = 0.0;
  for (int d1 = -(n1 - 1); d1 <= n1 - 1; d1++)
    for (int d2 = -(n2 - 1); d2 <= n2 - 1; d2++)
      E += (n1 - std::abs(d1)) * (n2 - std::abs(d2)) * C_lookup(d1, d2);
  return E / (static_cast<double>(n1 * n2) * static_cast<double>(n1 * n2));
}

// Factored triple sum used in terms 12-15.
inline double triple_sum_factored(int a, int b, double E, double N) {
  return 3.0 * N * N * E * S1_lookup(a, b);
}


// ============================================================================
// SECTION 3: Initialise all lookup tables at once
// ============================================================================

// Master initialiser — call at the start of every exported function.
// Skips reinitialisation if parameters haven't changed.
static void init_lookups(int n1, int n2, double beta,
                         double lambda1, double lambda2,
                         double alpha1, double alpha2,
                         double sigma) {
  if (n1 == cached_n1 && n2 == cached_n2 &&
      beta == cached_beta &&
      lambda1 == cached_l1 && lambda2 == cached_l2 &&
      alpha1 == cached_a1 && alpha2 == cached_a2 &&
      sigma == cached_sigma) {
    return;
  }
  init_C_matrix(n1, n2, beta, lambda1, lambda2, alpha1, alpha2, sigma);
  init_S1_matrix(n1, n2);
}


// ============================================================================
// SECTION 4: Internal engines (not exported to R)
// ============================================================================

// --- Cov(Chat(h), Chat(g)) with nonzero mean ---------------------------------

// Full 16-term Isserlis expansion for a single (h, g) pair.
// Assumes lookup tables are already initialised.
double compute_cov_entry_internal(int n1, int n2,
                                   int h1, int h2,
                                   int h3, int h4,
                                   double mu, double E, double N) {

  double N_h = (double)((n1 - std::abs(h1)) * (n2 - std::abs(h2)));
  double N_g = (double)((n1 - std::abs(h3)) * (n2 - std::abs(h4)));

  // Overlap bounds for lag h
  int T11 = (h1 >= 0) ? 1 : 1 - h1;
  int T12 = (h1 >= 0) ? n1 - h1 : n1;
  int T21 = (h2 >= 0) ? 1 : 1 - h2;
  int T22 = (h2 >= 0) ? n2 - h2 : n2;

  // Overlap bounds for lag g
  int S11 = (h3 >= 0) ? 1 : 1 - h3;
  int S12 = (h3 >= 0) ? n1 - h3 : n1;
  int S21 = (h4 >= 0) ? 1 : 1 - h4;
  int S22 = (h4 >= 0) ? n2 - h4 : n2;

  double mu2 = mu * mu;
  double mu4 = mu2 * mu2;

  double C_h1h2 = C_lookup(h1, h2);
  double C_h3h4 = C_lookup(h3, h4);

  // Precompute B_h, D_h (partial sums over the h-window)
  double B_h = 0.0, D_h = 0.0;
  for (int i = T11; i <= T12; i++)
    for (int j = T21; j <= T22; j++) {
      B_h += S1_lookup(i + h1, j + h2);
      D_h += S1_lookup(i, j);
    }
  B_h /= (N_h * N);
  D_h /= (N_h * N);

  // Precompute B_g, D_g (partial sums over the g-window)
  double B_g = 0.0, D_g = 0.0;
  for (int k = S11; k <= S12; k++)
    for (int ell = S21; ell <= S22; ell++) {
      B_g += S1_lookup(k + h3, ell + h4);
      D_g += S1_lookup(k, ell);
    }
  B_g /= (N_g * N);
  D_g /= (N_g * N);

  double result = 0.0;

  // Quadruple loop over (i, j, k, ell)
  for (int i = T11; i <= T12; i++) {
    for (int j = T21; j <= T22; j++) {
      for (int k = S11; k <= S12; k++) {
        for (int ell = S21; ell <= S22; ell++) {

          int di = k - i;
          int dj = ell - j;

          double C_di_dj = C_lookup(di, dj);
          double C_di_h1 = C_lookup(di - h1, dj - h2);
          double C_di_h3 = C_lookup(di + h3, dj + h4);
          double C_digh  = C_lookup(di + h3 - h1, dj + h4 - h2);

          double S1_i_j    = S1_lookup(i, j);
          double S1_ih_jh  = S1_lookup(i + h1, j + h2);
          double S1_k_ell  = S1_lookup(k, ell);
          double S1_kh_eh  = S1_lookup(k + h3, ell + h4);

          // Term 1 (+)
          result += mu4;
          result += mu2 * (C_h1h2 + C_h3h4 + C_di_dj + C_digh + C_di_h1 + C_di_h3);
          result += C_digh * C_di_dj + C_di_h1 * C_di_h3;

          // Term 2 (-)
          result -= mu4;
          result -= (mu2 / N) * (S1_i_j + S1_kh_eh + S1_k_ell);
          result -= mu2 * (C_di_h3 + C_di_dj + C_h3h4);
          result -= (1.0 / N) * (S1_i_j * C_h3h4 + S1_kh_eh * C_di_dj + S1_k_ell * C_di_h3);

          // Term 3 (-)
          result -= mu4;
          result -= (mu2 / N) * (S1_ih_jh + S1_kh_eh + S1_k_ell);
          result -= mu2 * (C_digh + C_di_h1 + C_h3h4);
          result -= (1.0 / N) * (S1_ih_jh * C_h3h4 + S1_kh_eh * C_di_h1 + S1_k_ell * C_digh);

          // Term 4 (-)
          result -= mu4;
          result -= (mu2 / N) * (S1_ih_jh + S1_i_j + S1_k_ell);
          result -= mu2 * (C_h1h2 + C_di_h1 + C_di_dj);
          result -= (1.0 / N) * (S1_ih_jh * C_di_dj + S1_i_j * C_di_h1 + S1_k_ell * C_h1h2);

          // Term 5 (-)
          result -= mu4;
          result -= (mu2 / N) * (S1_ih_jh + S1_i_j + S1_kh_eh);
          result -= mu2 * (C_h1h2 + C_digh + C_di_h3);
          result -= (1.0 / N) * (S1_ih_jh * C_di_h3 + S1_i_j * C_digh + S1_kh_eh * C_h1h2);

          // Term 6 (+)
          result += mu4 + mu2 * C_h3h4 + (mu2 + C_h3h4) * E;
          result += (2.0 * mu2 / N) * (S1_kh_eh + S1_k_ell);
          result += (2.0 / (N * N)) * S1_kh_eh * S1_k_ell;

          // Term 7 (+)
          result += mu4 + mu2 * C_di_dj + (mu2 + C_di_dj) * E;
          result += (2.0 * mu2 / N) * (S1_i_j + S1_k_ell);
          result += (2.0 / (N * N)) * S1_i_j * S1_k_ell;

          // Term 8 (+)
          result += mu4 + mu2 * C_di_h1 + (mu2 + C_di_h1) * E;
          result += (2.0 * mu2 / N) * (S1_ih_jh + S1_k_ell);
          result += (2.0 / (N * N)) * S1_ih_jh * S1_k_ell;

          // Term 9 (+)
          result += mu4 + mu2 * C_h1h2 + (mu2 + C_h1h2) * E;
          result += (2.0 * mu2 / N) * (S1_ih_jh + S1_i_j);
          result += (2.0 / (N * N)) * S1_ih_jh * S1_i_j;

          // Term 10 (+)
          result += mu4 + mu2 * C_digh + (mu2 + C_digh) * E;
          result += (2.0 * mu2 / N) * (S1_ih_jh + S1_kh_eh);
          result += (2.0 / (N * N)) * S1_ih_jh * S1_kh_eh;

          // Term 11 (+)
          result += mu4 + mu2 * C_di_h3 + (mu2 + C_di_h3) * E;
          result += (2.0 * mu2 / N) * (S1_i_j + S1_kh_eh);
          result += (2.0 / (N * N)) * S1_i_j * S1_kh_eh;

          // Term 12 (-)
          result -= mu4 + 3.0 * mu2 * E;
          result -= (3.0 * mu2 / N) * S1_ih_jh;
          result -= triple_sum_factored(i + h1, j + h2, E, N) / (N * N * N);

          // Term 13 (-)
          result -= mu4 + 3.0 * mu2 * E;
          result -= (3.0 * mu2 / N) * S1_i_j;
          result -= triple_sum_factored(i, j, E, N) / (N * N * N);

          // Term 14 (-)
          result -= mu4 + 3.0 * mu2 * E;
          result -= (3.0 * mu2 / N) * S1_kh_eh;
          result -= triple_sum_factored(k + h3, ell + h4, E, N) / (N * N * N);

          // Term 15 (-)
          result -= mu4 + 3.0 * mu2 * E;
          result -= (3.0 * mu2 / N) * S1_k_ell;
          result -= triple_sum_factored(k, ell, E, N) / (N * N * N);

          // Term 16 (+)
          result += mu4 + 6.0 * mu2 * E + 3.0 * E * E;

        } // ell
      } // k
    } // j
  } // i

  // Normalise by N_h * N_g
  result /= (N_h * N_g);

  // Correction from E[Chat(h)] * E[Chat(g)] subtraction
  result += C_h1h2 * B_g + C_h1h2 * D_g - C_h1h2 * E;
  result += C_h3h4 * B_h - B_h * B_g - B_h * D_g + B_h * E;
  result += C_h3h4 * D_h - D_h * B_g - D_h * D_g + D_h * E;
  result -= C_h3h4 * E   - E * B_g   - E * D_g   + E * E;

  return result;
}


// --- Cov(Chat(h), Chat(g)) with known zero mean ------------------------------

// Simplified formula when mu = 0: only the Isserlis product terms survive.
double compute_cov_entry_zeromean_internal(int n1, int n2,
                                            int h1, int h2,
                                            int h3, int h4) {

  double N_h = (double)((n1 - std::abs(h1)) * (n2 - std::abs(h2)));
  double N_g = (double)((n1 - std::abs(h3)) * (n2 - std::abs(h4)));

  int T11 = (h1 >= 0) ? 1 : 1 - h1;
  int T12 = (h1 >= 0) ? n1 - h1 : n1;
  int T21 = (h2 >= 0) ? 1 : 1 - h2;
  int T22 = (h2 >= 0) ? n2 - h2 : n2;

  int S11 = (h3 >= 0) ? 1 : 1 - h3;
  int S12 = (h3 >= 0) ? n1 - h3 : n1;
  int S21 = (h4 >= 0) ? 1 : 1 - h4;
  int S22 = (h4 >= 0) ? n2 - h4 : n2;

  double result = 0.0;
  for (int i = T11; i <= T12; i++)
    for (int j = T21; j <= T22; j++)
      for (int k = S11; k <= S12; k++)
        for (int ell = S21; ell <= S22; ell++) {
          int di = k - i, dj = ell - j;
          result += C_lookup(di + h3 - h1, dj + h4 - h2) * C_lookup(di, dj);
          result += C_lookup(di - h1, dj - h2) * C_lookup(di + h3, dj + h4);
        }

  return result / (N_h * N_g);
}


// ============================================================================
// SECTION 5: R-exported single-entry functions
// ============================================================================

// --- Cov(Chat(h), Chat(g)) with nonzero mean (scalar lags) -------------------

// Full Isserlis expansion for one pair of lag vectors.
// [[Rcpp::export]]
double cov_sample_autocov_opt(int n1, int n2,
                               int h1, int h2,
                               int h3, int h4,
                               double mu,
                               double beta,
                               NumericVector lambdas,
                               NumericVector alphas,
                               double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  double E = compute_E(n1, n2);

  return compute_cov_entry_internal(n1, n2, h1, h2, h3, h4, mu, E, N);
}


// --- Zero-mean version -------------------------------------------------------

// Cov(Chat(h), Chat(g)) when the process mean is known to be zero.
// [[Rcpp::export]]
double cov_sample_autocov_zeromean(int n1, int n2,
                                    int h1, int h2,
                                    int h3, int h4,
                                    double beta,
                                    NumericVector lambdas,
                                    NumericVector alphas,
                                    double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  return compute_cov_entry_zeromean_internal(n1, n2, h1, h2, h3, h4);
}


// --- Vector-lag convenience wrappers -----------------------------------------

// Cov(Chat(h), Chat(g)) taking lag vectors h = c(h1,h2) and g = c(g1,g2).
// [[Rcpp::export]]
double cov_Chat(int n1, int n2,
                NumericVector h, NumericVector g,
                double mu, double beta,
                NumericVector lambdas, NumericVector alphas,
                double sigma) {
  return cov_sample_autocov_opt(n1, n2,
                                 (int)h[0], (int)h[1],
                                 (int)g[0], (int)g[1],
                                 mu, beta, lambdas, alphas, sigma);
}

// Var(Chat(h)) = Cov(Chat(h), Chat(h)).
// [[Rcpp::export]]
double var_Chat(int n1, int n2,
                NumericVector h,
                double mu, double beta,
                NumericVector lambdas, NumericVector alphas,
                double sigma) {
  return cov_sample_autocov_opt(n1, n2,
                                 (int)h[0], (int)h[1],
                                 (int)h[0], (int)h[1],
                                 mu, beta, lambdas, alphas, sigma);
}


// --- Bias of Chat(h) --------------------------------------------------------

// Bias = E[Chat(h)] - C(h) = -B_h - D_h + E.
// [[Rcpp::export]]
double bias_Chat(int n1, int n2,
                 int h1, int h2,
                 double beta,
                 NumericVector lambdas, NumericVector alphas,
                 double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N   = (double)(n1 * n2);
  double N_h = (double)((n1 - std::abs(h1)) * (n2 - std::abs(h2)));
  double E   = compute_E(n1, n2);

  int T11 = (h1 >= 0) ? 1 : 1 - h1;
  int T12 = (h1 >= 0) ? n1 - h1 : n1;
  int T21 = (h2 >= 0) ? 1 : 1 - h2;
  int T22 = (h2 >= 0) ? n2 - h2 : n2;

  double B_h = 0.0, D_h = 0.0;
  for (int i = T11; i <= T12; i++)
    for (int j = T21; j <= T22; j++) {
      B_h += S1_lookup(i + h1, j + h2);
      D_h += S1_lookup(i, j);
    }
  B_h /= (N_h * N);
  D_h /= (N_h * N);

  return -B_h - D_h + E;
}


// ============================================================================
// SECTION 6: Omega matrix builders (covariance of vec(Gamma_hat))
// ============================================================================

// --- Truncated Omega_m for lags in {-m,...,m}^2 ------------------------------

// Build the (2m+1)^2 x (2m+1)^2 covariance matrix of vec(Gamma_hat_m).
// [[Rcpp::export]]
NumericMatrix build_Omega_m(int n1, int n2, int m,
                             double mu, double beta,
                             NumericVector lambdas,
                             NumericVector alphas,
                             double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  double E = compute_E(n1, n2);

  int L = 2 * m + 1;
  int p = L * L;
  NumericMatrix Omega(p, p);

  for (int k1 = 0; k1 < p; k1++) {
    int h1 = (k1 % L) - m;
    int h2 = (k1 / L) - m;
    for (int k2 = k1; k2 < p; k2++) {
      int h3 = (k2 % L) - m;
      int h4 = (k2 / L) - m;

      double val = compute_cov_entry_internal(n1, n2, h1, h2, h3, h4, mu, E, N);
      Omega(k1, k2) = val;
      Omega(k2, k1) = val;
    }
  }

  // Dimnames
  CharacterVector dn(p);
  for (int k = 0; k < p; k++)
    dn[k] = "C(" + std::to_string((k % L) - m) + "," +
                    std::to_string((k / L) - m) + ")";
  Omega.attr("dimnames") = List::create(dn, dn);

  return Omega;
}


// --- Full Omega for all lags {-(n1-1),...,n1-1} x {-(n2-1),...,n2-1} ---------

// Build the L1*L2 x L1*L2 covariance matrix of vec(Gamma_hat).
// [[Rcpp::export]]
NumericMatrix build_Omega_full(int n1, int n2,
                                double mu, double beta,
                                NumericVector lambdas,
                                NumericVector alphas,
                                double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  double E = compute_E(n1, n2);

  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;
  int o1 = n1 - 1, o2 = n2 - 1;
  int p = L1 * L2;

  NumericMatrix Omega(p, p);

  for (int k1 = 0; k1 < p; k1++) {
    int h1 = (k1 % L1) - o1;
    int h2 = (k1 / L1) - o2;
    for (int k2 = k1; k2 < p; k2++) {
      int h3 = (k2 % L1) - o1;
      int h4 = (k2 / L1) - o2;

      double val = compute_cov_entry_internal(n1, n2, h1, h2, h3, h4, mu, E, N);
      Omega(k1, k2) = val;
      Omega(k2, k1) = val;
    }
  }

  CharacterVector dn(p);
  for (int k = 0; k < p; k++)
    dn[k] = "C(" + std::to_string((k % L1) - o1) + "," +
                    std::to_string((k / L1) - o2) + ")";
  Omega.attr("dimnames") = List::create(dn, dn);

  return Omega;
}


// ============================================================================
// SECTION 7: Bias vector builders
// ============================================================================

// --- Truncated bias vector b_m = E[vec(Gamma_hat_m)] - vec(Gamma_m) ----------

// Build the (2m+1)^2 bias vector for lags in {-m,...,m}^2.
// [[Rcpp::export]]
NumericVector build_bias_m(int n1, int n2, int m,
                            double beta,
                            NumericVector lambdas,
                            NumericVector alphas,
                            double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  double E = compute_E(n1, n2);

  int L = 2 * m + 1;
  int p = L * L;
  NumericVector bvec(p);

  for (int k = 0; k < p; k++) {
    int h1 = (k % L) - m;
    int h2 = (k / L) - m;
    double N_h = (double)((n1 - std::abs(h1)) * (n2 - std::abs(h2)));

    int T11 = (h1 >= 0) ? 1 : 1 - h1;
    int T12 = (h1 >= 0) ? n1 - h1 : n1;
    int T21 = (h2 >= 0) ? 1 : 1 - h2;
    int T22 = (h2 >= 0) ? n2 - h2 : n2;

    double B_h = 0.0, D_h = 0.0;
    for (int i = T11; i <= T12; i++)
      for (int j = T21; j <= T22; j++) {
        B_h += S1_lookup(i + h1, j + h2);
        D_h += S1_lookup(i, j);
      }
    B_h /= (N_h * N);
    D_h /= (N_h * N);

    bvec[k] = -B_h - D_h + E;
  }

  CharacterVector names(p);
  for (int k = 0; k < p; k++)
    names[k] = "C(" + std::to_string((k % L) - m) + "," +
                       std::to_string((k / L) - m) + ")";
  bvec.attr("names") = names;

  return bvec;
}


// --- Full bias vector for all lags -------------------------------------------

// Build the L1*L2 bias vector for all lags.
// [[Rcpp::export]]
NumericVector build_bias_full(int n1, int n2,
                               double beta,
                               NumericVector lambdas,
                               NumericVector alphas,
                               double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  double E = compute_E(n1, n2);

  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;
  int o1 = n1 - 1, o2 = n2 - 1;
  int p = L1 * L2;

  NumericVector bvec(p);

  for (int k = 0; k < p; k++) {
    int h1 = (k % L1) - o1;
    int h2 = (k / L1) - o2;
    double N_h = (double)((n1 - std::abs(h1)) * (n2 - std::abs(h2)));

    int T11 = (h1 >= 0) ? 1 : 1 - h1;
    int T12 = (h1 >= 0) ? n1 - h1 : n1;
    int T21 = (h2 >= 0) ? 1 : 1 - h2;
    int T22 = (h2 >= 0) ? n2 - h2 : n2;

    double B_h = 0.0, D_h = 0.0;
    for (int i = T11; i <= T12; i++)
      for (int j = T21; j <= T22; j++) {
        B_h += S1_lookup(i + h1, j + h2);
        D_h += S1_lookup(i, j);
      }
    B_h /= (N_h * N);
    D_h /= (N_h * N);

    bvec[k] = -B_h - D_h + E;
  }

  CharacterVector names(p);
  for (int k = 0; k < p; k++)
    names[k] = "C(" + std::to_string((k % L1) - o1) + "," +
                       std::to_string((k / L1) - o2) + ")";
  bvec.attr("names") = names;

  return bvec;
}


// ============================================================================
// SECTION 8: MSE matrix builders (MSE = Omega + b * b')
// ============================================================================

// --- Truncated MSE_m ---------------------------------------------------------

// Build the (2m+1)^2 x (2m+1)^2 MSE matrix for lags in {-m,...,m}^2.
// [[Rcpp::export]]
NumericMatrix build_MSE_m(int n1, int n2, int m,
                           double mu, double beta,
                           NumericVector lambdas,
                           NumericVector alphas,
                           double sigma,
                           bool add_dimnames = false) {

  NumericMatrix Omega = build_Omega_m(n1, n2, m, mu, beta, lambdas, alphas, sigma);
  NumericVector bias  = build_bias_m(n1, n2, m, beta, lambdas, alphas, sigma);

  int p = bias.size();
  NumericMatrix MSE(p, p);

  for (int i = 0; i < p; i++)
    for (int j = 0; j < p; j++)
      MSE(i, j) = Omega(i, j) + bias[i] * bias[j];

  if (add_dimnames) {
    int L = 2 * m + 1;
    CharacterVector dn(p);
    for (int k = 0; k < p; k++)
      dn[k] = "C(" + std::to_string((k % L) - m) + "," +
                      std::to_string((k / L) - m) + ")";
    MSE.attr("dimnames") = List::create(dn, dn);
  }
  return MSE;
}


// --- Full MSE for all lags ---------------------------------------------------

// Build the L1*L2 x L1*L2 MSE matrix for all lags.
// [[Rcpp::export]]
NumericMatrix build_MSE_full(int n1, int n2,
                              double mu, double beta,
                              NumericVector lambdas,
                              NumericVector alphas,
                              double sigma,
                              bool add_dimnames = false) {

  NumericMatrix Omega = build_Omega_full(n1, n2, mu, beta, lambdas, alphas, sigma);
  NumericVector bias  = build_bias_full(n1, n2, beta, lambdas, alphas, sigma);

  int p = bias.size();
  NumericMatrix MSE(p, p);

  for (int i = 0; i < p; i++)
    for (int j = 0; j < p; j++)
      MSE(i, j) = Omega(i, j) + bias[i] * bias[j];

  if (add_dimnames) {
    int L1 = 2 * n1 - 1;
    int o1 = n1 - 1, o2 = n2 - 1;
    CharacterVector dn(p);
    for (int k = 0; k < p; k++)
      dn[k] = "C(" + std::to_string((k % L1) - o1) + "," +
                      std::to_string((k / L1) - o2) + ")";
    MSE.attr("dimnames") = List::create(dn, dn);
  }
  return MSE;
}


// ============================================================================
// SECTION 9: Expand truncated MSE_m into the full lag grid
// ============================================================================

// --- Direct expansion from (n1, n2, m) parameters ----------------------------

// Embed MSE_m into the full L1*L2 x L1*L2 matrix.
// Inner block: MSE_m.  Outer entries: bias(h) * bias(g).
// [[Rcpp::export]]
NumericMatrix expand_MSE_m_to_full(int n1, int n2, int m,
                                    double mu, double beta,
                                    NumericVector lambdas,
                                    NumericVector alphas,
                                    double sigma,
                                    bool add_dimnames = false) {

  NumericMatrix MSE_m   = build_MSE_m(n1, n2, m, mu, beta, lambdas, alphas, sigma, false);
  NumericVector bias_f  = build_bias_full(n1, n2, beta, lambdas, alphas, sigma);

  int L1 = 2 * n1 - 1, L2 = 2 * n2 - 1;
  int o1 = n1 - 1, o2 = n2 - 1;
  int p_full = L1 * L2;
  int L_m = 2 * m + 1;

  NumericMatrix MSE_exp(p_full, p_full);

  for (int k1 = 0; k1 < p_full; k1++) {
    int h1 = (k1 % L1) - o1, h2 = (k1 / L1) - o2;
    bool h_in = (std::abs(h1) <= m) && (std::abs(h2) <= m);

    for (int k2 = k1; k2 < p_full; k2++) {
      int h3 = (k2 % L1) - o1, h4 = (k2 / L1) - o2;
      bool g_in = (std::abs(h3) <= m) && (std::abs(h4) <= m);

      double val;
      if (h_in && g_in) {
        int k1_m = (h2 + m) * L_m + (h1 + m);
        int k2_m = (h4 + m) * L_m + (h3 + m);
        val = MSE_m(k1_m, k2_m);
      } else {
        val = bias_f[k1] * bias_f[k2];
      }

      MSE_exp(k1, k2) = val;
      MSE_exp(k2, k1) = val;
    }
  }

  if (add_dimnames) {
    CharacterVector dn(p_full);
    for (int k = 0; k < p_full; k++)
      dn[k] = "C(" + std::to_string((k % L1) - o1) + "," +
                      std::to_string((k / L1) - o2) + ")";
    MSE_exp.attr("dimnames") = List::create(dn, dn);
  }
  return MSE_exp;
}


// --- Expansion from pre-transformed components --------------------------------

// Like above, but takes already-transformed Omega_m and bias_m as inputs.
// Useful when W*_m^{-1} has already been applied on the R side.
// [[Rcpp::export]]
NumericMatrix expand_MSE_m_to_full_from_components(
    NumericMatrix Omega_m_transformed,
    NumericVector bias_m_transformed,
    int n1, int n2, int m,
    double beta,
    NumericVector lambdas,
    NumericVector alphas,
    double sigma,
    bool add_dimnames = false) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  int L1 = 2 * n1 - 1, L2 = 2 * n2 - 1;
  int o1 = n1 - 1, o2 = n2 - 1;
  int p_full = L1 * L2;
  int L_m = 2 * m + 1;

  NumericMatrix MSE_exp(p_full, p_full);

  for (int k1 = 0; k1 < p_full; k1++) {
    int h1 = (k1 % L1) - o1, h2 = (k1 / L1) - o2;
    bool h_in = (std::abs(h1) <= m) && (std::abs(h2) <= m);
    double C_h = C_lookup(h1, h2);

    for (int k2 = k1; k2 < p_full; k2++) {
      int h3 = (k2 % L1) - o1, h4 = (k2 / L1) - o2;
      bool g_in = (std::abs(h3) <= m) && (std::abs(h4) <= m);

      double val;
      if (h_in && g_in) {
        int k1_m = (h2 + m) * L_m + (h1 + m);
        int k2_m = (h4 + m) * L_m + (h3 + m);
        val = Omega_m_transformed(k1_m, k2_m) +
              bias_m_transformed[k1_m] * bias_m_transformed[k2_m];
      } else {
        val = C_h * C_lookup(h3, h4);
      }

      MSE_exp(k1, k2) = val;
      MSE_exp(k2, k1) = val;
    }
  }

  if (add_dimnames) {
    CharacterVector dn(p_full);
    for (int k = 0; k < p_full; k++)
      dn[k] = "C(" + std::to_string((k % L1) - o1) + "," +
                      std::to_string((k / L1) - o2) + ")";
    MSE_exp.attr("dimnames") = List::create(dn, dn);
  }
  return MSE_exp;
}


// ============================================================================
// SECTION 10: Compact (mu-free) covariance formula
//
// This implements the simplified closed-form expression for
//   Cov(Chat(h), Chat(g))
// that does NOT depend on the true mean mu. It collects the 16 Isserlis
// terms into 6 groups using the four anchor points
//   p1 = (i+h1, j+h2),  p2 = (i, j),
//   p3 = (k+h3, l+h4),  p4 = (k, l)
// and the operators
//   Phi(p) = sum_{t1,t2} C(p - t),   calE = sum_{t,s} C(t - s).
//
// The formula is:
//   (1/N_h N_g) sum_{i,j,k,l} {
//     C(d+g-h)*C(d) + C(d-h)*C(d+g)                         [line 1]
//     - (1/N) [ Phi(p2)(C(d)+C(d+g)) + Phi(p1)(C(d-h)+C(d+g-h))
//              + Phi(p4)(C(d)+C(d-h)) + Phi(p3)(C(d+g)+C(d+g-h)) ]  [line 2-3]
//     + (calE/N^2) [ C(d) + C(d-h) + C(d+g-h) + C(d+g) ]   [line 4]
//     + (1/N^2) [ 2*Phi(p2)*Phi(p1) + Phi(p2)*Phi(p4) + Phi(p2)*Phi(p3)
//               + Phi(p1)*Phi(p4) + Phi(p1)*Phi(p3) + 2*Phi(p4)*Phi(p3) ]  [line 5-6]
//     - (2*calE/N^3) [ Phi(p1) + Phi(p2) + Phi(p3) + Phi(p4) ]  [line 7]
//     + 2*calE^2/N^4                                          [line 8]
//   }
//
// Notably, the result does not depend on the true mean mu.
// Use these functions to verify against the 16-term version (Sections 4-5).
// ============================================================================


// --- Internal engine for compact formula -------------------------------------

// Compact Cov(Chat(h), Chat(g)) — does not depend on mu.
// Assumes lookup tables are already initialised.
// calE is the unnormalised double sum: calE = N^2 * E.
double compute_cov_entry_compact_internal(int n1, int n2,
                                           int h1, int h2,
                                           int h3, int h4,
                                           double calE, double N) {

  double N_h = (double)((n1 - std::abs(h1)) * (n2 - std::abs(h2)));
  double N_g = (double)((n1 - std::abs(h3)) * (n2 - std::abs(h4)));

  // Overlap bounds for lag h
  int T11 = (h1 >= 0) ? 1 : 1 - h1;
  int T12 = (h1 >= 0) ? n1 - h1 : n1;
  int T21 = (h2 >= 0) ? 1 : 1 - h2;
  int T22 = (h2 >= 0) ? n2 - h2 : n2;

  // Overlap bounds for lag g
  int S11 = (h3 >= 0) ? 1 : 1 - h3;
  int S12 = (h3 >= 0) ? n1 - h3 : n1;
  int S21 = (h4 >= 0) ? 1 : 1 - h4;
  int S22 = (h4 >= 0) ? n2 - h4 : n2;

  double N2 = N * N;
  double N3 = N2 * N;
  double N4 = N2 * N2;

  double result = 0.0;

  for (int i = T11; i <= T12; i++) {
    for (int j = T21; j <= T22; j++) {
      for (int k = S11; k <= S12; k++) {
        for (int ell = S21; ell <= S22; ell++) {

          int di = k - i;
          int dj = ell - j;

          // Four covariance values at difference d with lag shifts
          double C_d       = C_lookup(di, dj);
          double C_d_mh    = C_lookup(di - h1, dj - h2);
          double C_d_pg    = C_lookup(di + h3, dj + h4);
          double C_d_pgmh  = C_lookup(di + h3 - h1, dj + h4 - h2);

          // Four Phi values at the anchor points
          double Phi_p1 = S1_lookup(i + h1, j + h2);   // Phi(i+h1, j+h2)
          double Phi_p2 = S1_lookup(i, j);              // Phi(i, j)
          double Phi_p3 = S1_lookup(k + h3, ell + h4);  // Phi(k+h3, l+h4)
          double Phi_p4 = S1_lookup(k, ell);             // Phi(k, l)

          // Line 1: Isserlis product terms
          double val = C_d_pgmh * C_d + C_d_mh * C_d_pg;

          // Lines 2-3: -(1/N) * [Phi * C cross terms]
          val -= (1.0 / N) * (
            Phi_p2 * (C_d + C_d_pg) +
            Phi_p1 * (C_d_mh + C_d_pgmh) +
            Phi_p4 * (C_d + C_d_mh) +
            Phi_p3 * (C_d_pg + C_d_pgmh)
          );

          // Line 4: +(calE/N^2) * [sum of four C values]
          val += (calE / N2) * (C_d + C_d_mh + C_d_pgmh + C_d_pg);

          // Lines 5-6: +(1/N^2) * [Phi * Phi cross terms]
          val += (1.0 / N2) * (
            2.0 * Phi_p2 * Phi_p1 +
                  Phi_p2 * Phi_p4 +
                  Phi_p2 * Phi_p3 +
                  Phi_p1 * Phi_p4 +
                  Phi_p1 * Phi_p3 +
            2.0 * Phi_p4 * Phi_p3
          );

          // Line 7: -(2*calE/N^3) * [sum of four Phi values]
          val -= (2.0 * calE / N3) * (Phi_p1 + Phi_p2 + Phi_p3 + Phi_p4);

          // Line 8: +2*calE^2/N^4
          val += 2.0 * calE * calE / N4;

          result += val;

        } // ell
      } // k
    } // j
  } // i

  return result / (N_h * N_g);
}


// --- R-exported single entry --------------------------------------------------

// Compact Cov(Chat(h), Chat(g)) — does not depend on mu.
// [[Rcpp::export]]
double cov_Chat_compact(int n1, int n2,
                         NumericVector h, NumericVector g,
                         double beta,
                         NumericVector lambdas, NumericVector alphas,
                         double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  // calE = N^2 * E, where compute_E already divides by N^2
  double calE = N * N * compute_E(n1, n2);

  return compute_cov_entry_compact_internal(n1, n2,
                                             (int)h[0], (int)h[1],
                                             (int)g[0], (int)g[1],
                                             calE, N);
}

// Var(Chat(h)) using the compact formula.
// [[Rcpp::export]]
double var_Chat_compact(int n1, int n2,
                         NumericVector h,
                         double beta,
                         NumericVector lambdas, NumericVector alphas,
                         double sigma) {

  return cov_Chat_compact(n1, n2, h, h, beta, lambdas, alphas, sigma);
}


// --- Compact Omega_m builder --------------------------------------------------

// Build (2m+1)^2 x (2m+1)^2 Omega using the compact formula (mu-free).
// [[Rcpp::export]]
NumericMatrix build_Omega_m_compact(int n1, int n2, int m,
                                     double beta,
                                     NumericVector lambdas,
                                     NumericVector alphas,
                                     double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  double calE = N * N * compute_E(n1, n2);

  int L = 2 * m + 1;
  int p = L * L;
  NumericMatrix Omega(p, p);

  for (int k1 = 0; k1 < p; k1++) {
    int h1 = (k1 % L) - m;
    int h2 = (k1 / L) - m;
    for (int k2 = k1; k2 < p; k2++) {
      int h3 = (k2 % L) - m;
      int h4 = (k2 / L) - m;

      double val = compute_cov_entry_compact_internal(n1, n2, h1, h2, h3, h4, calE, N);
      Omega(k1, k2) = val;
      Omega(k2, k1) = val;
    }
  }

  CharacterVector dn(p);
  for (int k = 0; k < p; k++)
    dn[k] = "C(" + std::to_string((k % L) - m) + "," +
                    std::to_string((k / L) - m) + ")";
  Omega.attr("dimnames") = List::create(dn, dn);

  return Omega;
}


// --- Compact Omega_full builder -----------------------------------------------

// Build L1*L2 x L1*L2 Omega using the compact formula (mu-free).
// [[Rcpp::export]]
NumericMatrix build_Omega_full_compact(int n1, int n2,
                                        double beta,
                                        NumericVector lambdas,
                                        NumericVector alphas,
                                        double sigma) {

  init_lookups(n1, n2, beta, lambdas[0], lambdas[1],
               alphas[0], alphas[1], sigma);

  double N = (double)(n1 * n2);
  double calE = N * N * compute_E(n1, n2);

  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;
  int o1 = n1 - 1, o2 = n2 - 1;
  int p = L1 * L2;

  NumericMatrix Omega(p, p);

  for (int k1 = 0; k1 < p; k1++) {
    int h1 = (k1 % L1) - o1;
    int h2 = (k1 / L1) - o2;
    for (int k2 = k1; k2 < p; k2++) {
      int h3 = (k2 % L1) - o1;
      int h4 = (k2 / L1) - o2;

      double val = compute_cov_entry_compact_internal(n1, n2, h1, h2, h3, h4, calE, N);
      Omega(k1, k2) = val;
      Omega(k2, k1) = val;
    }
  }

  CharacterVector dn(p);
  for (int k = 0; k < p; k++)
    dn[k] = "C(" + std::to_string((k % L1) - o1) + "," +
                    std::to_string((k / L1) - o2) + ")";
  Omega.attr("dimnames") = List::create(dn, dn);

  return Omega;
}


// --- Compact MSE_m builder (MSE = Omega_compact + b * b') --------------------

// Build (2m+1)^2 x (2m+1)^2 MSE using the compact (mu-free) Omega.
// [[Rcpp::export]]
NumericMatrix build_MSE_m_compact(int n1, int n2, int m,
                                   double beta,
                                   NumericVector lambdas,
                                   NumericVector alphas,
                                   double sigma,
                                   bool add_dimnames = false) {

  NumericMatrix Omega = build_Omega_m_compact(n1, n2, m, beta, lambdas, alphas, sigma);
  NumericVector bias  = build_bias_m(n1, n2, m, beta, lambdas, alphas, sigma);

  int p = bias.size();
  NumericMatrix MSE(p, p);

  for (int i = 0; i < p; i++)
    for (int j = 0; j < p; j++)
      MSE(i, j) = Omega(i, j) + bias[i] * bias[j];

  if (add_dimnames) {
    int L = 2 * m + 1;
    CharacterVector dn(p);
    for (int k = 0; k < p; k++)
      dn[k] = "C(" + std::to_string((k % L) - m) + "," +
                      std::to_string((k / L) - m) + ")";
    MSE.attr("dimnames") = List::create(dn, dn);
  }
  return MSE;
}


// --- Compact MSE_full builder ------------------------------------------------

// Build L1*L2 x L1*L2 MSE using the compact (mu-free) Omega.
// [[Rcpp::export]]
NumericMatrix build_MSE_full_compact(int n1, int n2,
                                      double beta,
                                      NumericVector lambdas,
                                      NumericVector alphas,
                                      double sigma,
                                      bool add_dimnames = false) {

  NumericMatrix Omega = build_Omega_full_compact(n1, n2, beta, lambdas, alphas, sigma);
  NumericVector bias  = build_bias_full(n1, n2, beta, lambdas, alphas, sigma);

  int p = bias.size();
  NumericMatrix MSE(p, p);

  for (int i = 0; i < p; i++)
    for (int j = 0; j < p; j++)
      MSE(i, j) = Omega(i, j) + bias[i] * bias[j];

  if (add_dimnames) {
    int L1 = 2 * n1 - 1;
    int o1 = n1 - 1, o2 = n2 - 1;
    CharacterVector dn(p);
    for (int k = 0; k < p; k++)
      dn[k] = "C(" + std::to_string((k % L1) - o1) + "," +
                      std::to_string((k / L1) - o2) + ")";
    MSE.attr("dimnames") = List::create(dn, dn);
  }
  return MSE;
}