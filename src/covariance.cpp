// ============================================================================
// Author: Antonio Fioravanti
// True autocovariance model for the modified exponential.
// Provides R-callable functions to evaluate C(h) at a single lag, build the
// full L1 x L2 matrix of true autocovariances, build the covariance matrix
// on a grid.
// ============================================================================

#include <Rcpp.h>
#include "covariance.h"
using namespace Rcpp;
//-----------------------------------------------------------------------------

// Function for R 
// [[Rcpp::export]]
double eval_true_cov(NumericVector hvec,
                        double beta,
                        NumericVector lambdas,
                        NumericVector alphas,
                        double sigma) {
// Evaluates the true covariance given by:
//   C(h1,h2) = sigma2 * exp(-(|h1|^alpha1 / lambda1
//                              + |h2|^alpha2 / lambda2
//                              + beta * |h1 - h2| / lambda1))
// at the lag vector hvec = (h1, h2).
  if (hvec.size() != 2)
    stop("hvec must have length 2");
  if (lambdas.size() != 2)
    stop("lambdas must have length 2");
  if (alphas.size() != 2)
    stop("alphas must have length 2");

  return true_cov(hvec[0], hvec[1],
                  beta,
                  lambdas[0], lambdas[1],
                  alphas[0], alphas[1],
                  sigma);
}

// ----------------------------------------------------------------------------
//  Full covariance matrix on a lattice grid

// Builds the N x N covariance matrix Sigma for the grid locations in the grid.
// Parameters:
//   grid: a DataFrame with columns t1 and t2 for the coordinates of the N
//         grid points.
//   sigma, alpha1, alpha2, lambda1, lambda2, beta: parameters of the modified
//   exponential model from eval_true_cov.
// [[Rcpp::export]]
NumericMatrix ModifiedExponentialCovariance(
    DataFrame grid,
    double sigma   = 1.0,
    double alpha1  = 1.0,
    double alpha2  = 1.0,
    double lambda1 = 1.0,
    double lambda2 = 1.0,
    double beta    = 0.0) {

  NumericVector t1 = grid["t1"];
  NumericVector t2 = grid["t2"];
  int n = t1.size();
  if (t2.size() != n)
    stop("grid$t1 and grid$t2 must have the same length");

  NumericMatrix cov(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cov(i, j) = true_cov(t1[i] - t1[j], t2[i] - t2[j],
                                beta, lambda1, lambda2,
                                alpha1, alpha2, sigma);
    }
  }
  return cov;
}


// ----------------------------------------------------------------------------
// Lag indexed matrix Gamma of true autocovariances


// [[Rcpp::export]]
NumericMatrix build_true_Gamma(NumericVector nvec,
                               double beta,
                               NumericVector lambdas,
                               NumericVector alphas,
                               double sigma) {
// Builds the (2*n1-1) x (2*n2-1) matrix Gamma (h1, h2) for all lags h1, h2.
// Rows indexed by h1 in {-(n1-1), ..., n1-1},
// columns by h2 in {-(n2-1), ..., n2-1}.
// Parameters:
//   nvec = (n1, n2) for the grid size (n1 x n2).
//   beta, lambdas, alphas, sigma: parameters of the modified exponential model
//   from eval_true_cov.

  if (nvec.size() != 2)
    stop("nvec must have length 2");
  if (lambdas.size() != 2)
    stop("lambdas must have length 2");
  if (alphas.size() != 2)
    stop("alphas must have length 2");

  int n1 = (int)nvec[0], n2 = (int)nvec[1];
  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;

  NumericMatrix Ctrue(L1, L2);

  for (int i = 0; i < L1; ++i) {
    int h1 = i - (n1 - 1);
    for (int j = 0; j < L2; ++j) {
      int h2 = j - (n2 - 1);
      Ctrue(i, j) = true_cov(h1, h2, beta,
                              lambdas[0], lambdas[1],
                              alphas[0], alphas[1],
                              sigma);
    }
  }

  // Dimnames
  CharacterVector rown(L1), coln(L2);
  for (int i = 0; i < L1; ++i)
    rown[i] = "h1=" + std::to_string(i - (n1 - 1));
  for (int j = 0; j < L2; ++j)
    coln[j] = "h2=" + std::to_string(j - (n2 - 1));
  Ctrue.attr("dimnames") = List::create(rown, coln);

  return Ctrue;
}
