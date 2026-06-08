//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// Code to build the counting factors A, B, D and the
// counting matrices used to construct the weight matrix W*_n.
//
// Purely combinatorial — no dependence on the covariance model.
//-----------------------------------------------------------------------------

#include "counting.h"
#include <algorithm>
#include <cmath>


//-----------------------------------------------------------------------------
// Scalar counting factors A, B, D
//-----------------------------------------------------------------------------

// Count of (i, t) pairs producing lag r in the A-term).
// [[Rcpp::export]]
int A_count(int n, int h, int r) {
  int val = n - std::max({std::abs(r), std::abs(h), std::abs(r + h)});
  return val > 0 ? val : 0; // condition ? if true : if false
}

// Count of (i, t) pairs producing lag r in the B-term .
// [[Rcpp::export]]
int B_count(int n, int h, int r) {
  int val = n - std::max({std::abs(r), std::abs(h), std::abs(r - h)});
  return val > 0 ? val : 0; // condition ? if true : if false
}

// Count of (t1, t2) pairs producing lag r in the D-term.
// [[Rcpp::export]]
int D_count(int n, int r) {
  int val = n - std::abs(r);
  return val > 0 ? val : 0; // condition ? if true : if false
}


//-----------------------------------------------------------------------------
// Full counting matrices and their normalised versions
//-----------------------------------------------------------------------------

// Build L x L matrices A, B, D and their normalised versions A*, B*, D*
// for one spatial dimension of size n_k.
//
// L = 2*n_k - 1, with lag grid {-(n_k-1), ..., 0, ..., n_k-1}.
//
// Nh_normalize = true:  A*(h,r) = A(h,r) / (n_k * (n_k - |h|)),
//                        B*(h,r) = B(h,r) / (n_k * (n_k - |h|)),
//                        D*(h,r) = D(r) / n_k^2.
// Nh_normalize = false: all divided by n_k^2.
// [[Rcpp::export]]
List build_count_matrices(int n_k, bool Nh_normalize) {
  // Builds a list containing:
  //   A, B, D: L x L matrices of raw counts for the A, B, D terms.
  //   A_norm, B_norm, D_norm: L x L matrices of normalised counts for the A, B, D terms.
  //   L: the size of the lag grid (L = 2*n_k - 1).
  // Parameters:
  //   n_k: size of the spatial dimension (number of grid points).
  //   Nh_normalize: if true, normalise by N_h; if false, normalise by N = n_k^2.

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

  if (Nh_normalize) {
    NumericVector row_den(L);
    for (int i = 0; i < L; i++) {
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
  } else { // We normalise by N = n_k^2 for all entries
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        A_norm(i, j) = A(i, j) / N;
        B_norm(i, j) = B(i, j) / N;
        D_norm(i, j) = D(i, j) / N;
      }
    }
  }

  return List::create(
    Named("A")      = A,
    Named("B")      = B,
    Named("D")      = D,
    Named("I")      = I,
    Named("A_norm") = A_norm,
    Named("B_norm") = B_norm,
    Named("D_norm") = D_norm,
    Named("L")      = L
  );
}
