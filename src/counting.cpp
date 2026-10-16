//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// Code to build the counting factors A, B, D and the
// counting matrices used to construct the weight matrix W_n and W*_n.
//-----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "counting.h"
#include <algorithm>
#include <cmath>

//-----------------------------------------------------------------------------
// Scalar counting factors A, B, D
//-----------------------------------------------------------------------------

// Count of (i, t) pairs producing lag r in the A-term.
// Parameters:
//   n:    size of the spatial dimension.
//   h, r: lag and summation lag.
// Returns:
//   max(0, n - max(|r|, |h|, |r+h|)).
// [[Rcpp::export]]
int A_count(int n, int h, int r) {
  int val = n - std::max({std::abs(r), std::abs(h), std::abs(r + h)});
  int A_nrh = val > 0 ? val : 0;
  return A_nrh;
}

// Count of (i, t) pairs producing lag r in the B-term.
// Parameters:
//   n:    size of the spatial dimension.
//   h, r: lag and summation lag.
// Returns:
//   max(0, n - max(|r|, |h|, |r-h|)).
// [[Rcpp::export]]
int B_count(int n, int h, int r) {
  int val = n - std::max({std::abs(r), std::abs(h), std::abs(r - h)});
  int B_nrh = val > 0 ? val : 0;
  return B_nrh;
}

// Count of (t1, t2) pairs producing lag r in the D-term.
// Parameters:
//   n: size of the spatial dimension.
//   r: summation lag.
// Returns: 
//   max(0, n - |r|).
// [[Rcpp::export]]
int D_count(int n, int r) {
  int val = n - std::abs(r);
  int D_r = val > 0 ? val : 0;
  return D_r;
}


//-----------------------------------------------------------------------------
// COUNTING MATRICES A, B, D, DELTA AND THEIR NORMALIZED VERSIONS
//-----------------------------------------------------------------------------

// Build the 1D counting matrices over the lag grid (rows = h, cols = r).
// Parameters:
//   n_k: size of the spatial dimension.
// Returns:
//   A, B, D:                L x L raw counts; I: L x L identity.
//   A_norm, B_norm, D_norm: N_h normalised counts.
//   Delta, Delta_mat:       Bartlett window (n_k-|h|)/n_k
//   L:                      lag-grid size, L = 2*n_k - 1.
// [[Rcpp::export]]
List build_count_matrices(int n_k) {

  // Size of the lag grid
  int L = 2 * n_k - 1;

  // Lag grid
  IntegerVector lags(L);
  for (int i = 0; i < L; i++)
    lags[i] = i - (n_k - 1);

  // Raw counting matrices
  NumericMatrix A(L, L), B(L, L), D(L, L), I(L, L);

  for (int i = 0; i < L; i++) I(i, i) = 1.0;
  for (int i = 0; i < L; i++) {
    int h = lags[i];
    for (int j = 0; j < L; j++) {
      int r = lags[j];

      A(i, j) = A_count(n_k, h, r);
      B(i, j) = B_count(n_k, h, r);
      D(i, j) = D_count(n_k, r);
    }
  }

  // Normalised matrices
  double N = static_cast<double>(n_k) * static_cast<double>(n_k);
  NumericMatrix A_norm(L, L), B_norm(L, L), D_norm(L, L);

  NumericVector row_den(L);
    for (int i = 0; i < L; i++) {
      //n_k * (n_k - |h|)
      row_den[i] = static_cast<double>(n_k) *
                   (static_cast<double>(n_k) - std::abs(lags[i]));
    }
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        A_norm(i, j) = A(i, j) / row_den[i];
        B_norm(i, j) = B(i, j) / row_den[i];
        D_norm(i, j) = D(i, j) / N;
      }
    }
  // Vector Bartlett window Delta(h) = (n_k - |h|)/n_k 
  NumericVector Delta(L);
  for (int i = 0; i < L; i++)
    Delta[i] = (static_cast<double>(n_k) - std::abs(lags[i])) / static_cast<double>(n_k);

  NumericMatrix Delta_mat(L, L);
  for (int i = 0; i < L; i++)
    Delta_mat(i, i) = Delta[i];

  return List::create(
    Named("A")         = A,
    Named("B")         = B,
    Named("D")         = D,
    Named("I")         = I,
    Named("A_norm")    = A_norm,
    Named("B_norm")    = B_norm,
    Named("D_norm")    = D_norm,
    Named("Delta")     = Delta,
    Named("Delta_mat") = Delta_mat,
    Named("L")         = L
  );
}


//-----------------------------------------------------------------------------
//  SELECTION MATRICES FOR THE BIAS-CORRECTED ESTIMATOR
//-----------------------------------------------------------------------------

// P_m: picks the central lags {-m,...,m} out of the full grid; row p has a 1
// at q = ceil(L/2) - m + p - 1.
// Parameters:
//   L: size of the full lag grid (L = 2*n_k - 1).
//   m: truncation cut-off (0 <= m <= n_k-1).
// Returns:
//   the M x L selection matrix, M = 2m+1.
// [[Rcpp::export]]
arma::mat build_Pm(int L, int m) {
  int M = 2 * m + 1;  // number of selected lags
  arma::mat P = arma::zeros<arma::mat>(M, L);
  int center = static_cast<int>(std::ceil(L / 2.0)); // 1-based center (lag 0)
  for (int p = 1; p <= M; ++p) {                   // 1-based row index
    int q = center - m + p - 1;                    // 1-based column, per eq.
    // back to 0-based for storage reasons
    P(p - 1, q - 1) = 1.0;                          
  }
  return P;
}

// P^c_m: complement of P_m picks the OUTER lags, i.e. the full grid minus
// the central block {-m,...,m}. Row p has a single 1 at q = p or q = p + M
// Parameters:
//   L: size of the full lag grid (L = 2*n_k - 1).
//   m: truncation cut-off; the central block is the part REMOVED.
// Returns:
//   the K x L selection matrix, K = L - (2m+1) = 2*(n_k-m-1).
// [[Rcpp::export]]
arma::mat build_Pmc(int L, int m) {
  int n_k = (L + 1) / 2;                 // recover grid size from L = 2n_k - 1
  int M   = 2 * m + 1;                   // width of the removed central block
  int K   = L - M;                       // number of outer lags = 2*(n_k-m-1)
  arma::mat P = arma::zeros<arma::mat>(K, L);
  for (int p = 1; p <= K; ++p) {         // 1-based row index
    int q;
    if (p <= n_k - m - 1) {
      q = p;                             // left block: 
    } else {
      q = p + M;                         // right block:
    }
    P(p - 1, q - 1) = 1.0;               // back to 0-based for storage
  }
  return P;
}

//-----------------------------------------------------------------------------
// S_m: Selection matrix picking the lags {-m1,...,m1} x {-m2,...,m2}.
// Parameters:
//   L1, L2: lag-grid sizes, Li = 2*ni - 1.
//   m1:     cut-off for dimension 1.
//   m2:     cut-off for dimension 2. Default -1 => m2 = m1.
// Returns:
//   S_m = kron(P2, P1), of size (M1*M2) x (L1*L2), Mi = 2*mi + 1.
// [[Rcpp::export]]
arma::mat build_Sm(int L1, int L2, int m1, int m2 = -1) {
  if (m2 < 0) m2 = m1;                   // default: same cut-off in both dims
  arma::mat P1 = build_Pm(L1, m1);
  arma::mat P2 = build_Pm(L2, m2);
  return arma::kron(P2, P1);
}

//-----------------------------------------------------------------------------
// WEIGHT MATRICES for the bias-corrected estimator.

// Parameters:
//   mat1, mat2:   outputs of build_count_matrices(n1) and (n2).
//   m1, m2:       truncation cut-offs per dimension. Default -1 => m2 = m1.
//   Nh_normalize: TRUE = standard 1/N_h; FALSE = constant 1/N.
// Returns:
//   W, W_star:     full unnormalised / normalised weight matrices.
//   W_m, W_star_m: their truncations to the central window.
// [[Rcpp::export]]
List build_weights(List mat1, List mat2, int m1, int m2 = -1,
               bool Nh_normalize = true) {

  if (m2 < 0) m2 = m1; // default: same cut-off in both dims

  // Raw counts, normalised counts and grid sizes from the R lists
  arma::mat A1  = as<arma::mat>(mat1["A"]);
  arma::mat B1  = as<arma::mat>(mat1["B"]);
  arma::mat D1  = as<arma::mat>(mat1["D"]);
  arma::mat A1n = as<arma::mat>(mat1["A_norm"]);
  arma::mat B1n = as<arma::mat>(mat1["B_norm"]);
  arma::mat D1n = as<arma::mat>(mat1["D_norm"]);
  arma::mat Del1 = as<arma::mat>(mat1["Delta_mat"]);
  int L1        = as<int>(mat1["L"]);

  arma::mat A2  = as<arma::mat>(mat2["A"]);
  arma::mat B2  = as<arma::mat>(mat2["B"]);
  arma::mat D2  = as<arma::mat>(mat2["D"]);
  arma::mat A2n = as<arma::mat>(mat2["A_norm"]);
  arma::mat B2n = as<arma::mat>(mat2["B_norm"]);
  arma::mat D2n = as<arma::mat>(mat2["D_norm"]);
  arma::mat Del2 = as<arma::mat>(mat2["Delta_mat"]);
  int L2        = as<int>(mat2["L"]);

  int L12 = L1 * L2;                        // full lag-grid dimension
  arma::mat I = arma::eye<arma::mat>(L12, L12);

  // Full UNNORMALIZED weight matrix
  arma::mat W = I - arma::kron(A2, A1) - arma::kron(B2, B1) + arma::kron(D2, D1);

  if (!Nh_normalize) {
    // = diag(Delta2 (x) Delta1)
      arma::mat Delta_kron = arma::kron(Del2, Del1);   
      W = Delta_kron * W;
  }

  // Full NORMALIZED weight matrix
  arma::mat W_star = I - arma::kron(A2n, A1n) - arma::kron(B2n, B1n) + arma::kron(D2n, D1n);
  if (!Nh_normalize) {
    // = diag(Delta2 (x) Delta1)
      arma::mat Delta_kron = arma::kron(Del2, Del1);   
      W_star = Delta_kron * W_star;
  }

  arma::mat Sm = build_Sm(L1, L2, m1, m2);   // (M1*M2) x (L1*L2)
  arma::mat W_m = Sm * W * Sm.t();
  arma::mat W_star_m = Sm * W_star * Sm.t();
  

  return List::create(
    Named("W")        = Rcpp::wrap(W),
    Named("W_star")   = Rcpp::wrap(W_star),
    Named("W_m")      = Rcpp::wrap(W_m),
    Named("W_star_m") = Rcpp::wrap(W_star_m)
  );
}

