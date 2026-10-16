//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// True autocovariance models, selected by model = "modexp" | "gencauchy"
// Reference: https://rdrr.io/cran/geoR/man/cov.spatial.html
//
//   d1  = (|h1| / lambda1)^alpha1,   d2 = (|h2| / lambda2)^alpha2
//   d12 = (|h1 - h2| / sqrt(lambda1 * lambda2))^min(alpha1, alpha2)
//
//   modexp     beta = 0:  C = sigma^2 * exp(-d1) * exp(-d2)
//              beta > 0:  C = sigma^2 * exp(-(d1 + d2 + beta * d12))
//   gencauchy  beta = 0:  C = sigma^2 * (1 + d1)^(-delta/alpha1) * (1 + d2)^(-delta/alpha2)
//              beta > 0:  C = sigma^2 * (1 + d1 + d2 + beta * d12)^(-delta/ahar)
//              abar = 2 * alpha1 * alpha2 / (alpha1 + alpha2)   (harmonic mean)
//   where:
//   lambda: range parameter (0, inf)
//   alpha:  smoothness (or shape) parameter (0,2]
//   beta.   separability parameter beta=0 (separable) vs beta >0 (nonseparable)
//   delta: (gencauchy only) long memory (0,inf). Long memory for <=1.
//-----------------------------------------------------------------------------

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <algorithm>
using namespace Rcpp;

//-----------------------------------------------------------------------------
// Model id and parameter bundle
//-----------------------------------------------------------------------------
enum CovModel { COV_MODEXP = 0, COV_GENCAUCHY = 1 };

struct CovParams {
  double beta      = 0.0;
  double lambda1   = 1.0, lambda2 = 1.0;
  double alpha1    = 1.0, alpha2  = 1.0;
  double sigma     = 1.0;
  double delta     = 1.0;   // gencauchy only
  int    model     = COV_MODEXP;
};

static int parse_cov_model(std::string s) {
  // convert string into lower case
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (s == "modexp")    return COV_MODEXP;
  if (s == "gencauchy") return COV_GENCAUCHY;
  stop("unknown covariance model '" + s +
       "'; expected \"modexp\" or \"gencauchy\"");
  return COV_MODEXP; // never reached
}

static void validate_cov_params(const CovParams& p) {
  if (p.lambda1 <= 0.0 || p.lambda2 <= 0.0)
    stop("lambdas must be > 0");
  if (p.alpha1 <= 0.0 || p.alpha1 > 2.0 || p.alpha2 <= 0.0 || p.alpha2 > 2.0)
    stop("alphas must be in (0, 2] for positive-definiteness");
  if (p.beta < 0.0)
    stop("beta must be >= 0");
  if (p.model == COV_GENCAUCHY && p.delta <= 0.0)
    stop("delta must be > 0 for the gencauchy model");
}

// R arguments -> validated parameter bundle
static CovParams make_params(double beta,
                             NumericVector lambdas,
                             NumericVector alphas,
                             double sigma,
                             std::string model,
                             double delta) {
  if (lambdas.size() != 2) stop("lambdas must have length 2");
  if (alphas.size() != 2)  stop("alphas must have length 2");

  CovParams p;
  p.model   = parse_cov_model(model);
  p.sigma   = sigma;
  p.alpha1  = alphas[0];  p.alpha2  = alphas[1];
  p.lambda1 = lambdas[0]; p.lambda2 = lambdas[1];
  p.beta    = beta;
  p.delta   = delta;
  validate_cov_params(p);
  return p;
}

//-----------------------------------------------------------------------------
// The covariance functions
//-----------------------------------------------------------------------------

// True covariance at one lag (h1, h2)
static double true_cov(double h1, double h2, const CovParams& p) {
  double sigma2 = p.sigma * p.sigma;
  double d1 = std::pow(std::abs(h1) / p.lambda1, p.alpha1);
  double d2 = std::pow(std::abs(h2) / p.lambda2, p.alpha2);

  // coupling term, only used when beta > 0
  double a12 = std::min(p.alpha1, p.alpha2);
  double d12 = std::pow(std::abs(h1 - h2) / std::sqrt(p.lambda1 * p.lambda2), a12);

  if (p.model == COV_MODEXP) {
    if (p.beta == 0.0) {
      // separable
      return sigma2 * std::exp(-d1) * std::exp(-d2);
    } else {
      // non-separable
      return sigma2 * std::exp(-(d1 + d2 + p.beta * d12));
    }
  } else if (p.model == COV_GENCAUCHY) {
    if (p.beta == 0.0) {
      // separable: each dimension with its own alpha
      return sigma2 * std::pow(1.0 + d1, -p.delta / p.alpha1)
                    * std::pow(1.0 + d2, -p.delta / p.alpha2);
    } else {
      // non-separable: harmonic mean of the alphas in the exponent
      double ahar = 2.0 * p.alpha1 * p.alpha2 / (p.alpha1 + p.alpha2);
      return sigma2 * std::pow(1.0 + d1 + d2 + p.beta * d12, -p.delta / ahar);
    }
  }

  // Potentially extend to other covariances from 
  //  https://rdrr.io/cran/geoR/man/cov.spatial.html
  return NA_REAL;
}

//-----------------------------------------------------------------------------
// Function for R
//-----------------------------------------------------------------------------

// [[Rcpp::export]]
double eval_true_cov(NumericVector hvec,
                     double beta,
                     NumericVector lambdas,
                     NumericVector alphas,
                     double sigma,
                     std::string model = "modexp",
                     double delta = 1.0) {
// Evaluates the true covariance at the lag vector hvec = (h1, h2) under the
// model selected by `model`: "modexp" or "gencauchy".
  if (hvec.size() != 2)
    stop("hvec must have length 2");

  CovParams p = make_params(beta, lambdas, alphas, sigma, model, delta);

  return true_cov(hvec[0], hvec[1], p);
}


// Lag indexed matrix Gamma of true autocovariances

// [[Rcpp::export]]
NumericMatrix build_true_Gamma(NumericVector nvec,
                               double beta,
                               NumericVector lambdas,
                               NumericVector alphas,
                               double sigma,
                               std::string model = "modexp",
                               double delta = 1.0) {
// Builds the (2*n1-1) x (2*n2-1) matrix Gamma(h1, h2) for all lags h1, h2.
// Rows indexed by h1; Columns by h2 
  if (nvec.size() != 2)
    stop("nvec must have length 2");

  CovParams p = make_params(beta, lambdas, alphas, sigma, model, delta);

  int n1 = (int)nvec[0], n2 = (int)nvec[1];
  int L1 = 2 * n1 - 1;
  int L2 = 2 * n2 - 1;

  NumericMatrix Ctrue(L1, L2);

  for (int i = 0; i < L1; ++i) {
    int h1 = i - (n1 - 1);
    for (int j = 0; j < L2; ++j) {
      int h2 = j - (n2 - 1);
      Ctrue(i, j) = true_cov(h1, h2, p);
    }
  }

  CharacterVector rown(L1), coln(L2);
  for (int i = 0; i < L1; ++i)
    rown[i] = "h1=" + std::to_string(i - (n1 - 1));
  for (int j = 0; j < L2; ++j)
    coln[j] = "h2=" + std::to_string(j - (n2 - 1));
  Ctrue.attr("dimnames") = List::create(rown, coln);

  return Ctrue;
}

//-----------------------------------------------------------------------------
// True variance-covariance matrix Sigma of the n1 x n2 lattice
//-----------------------------------------------------------------------------

// N x N covariance matrix (N = n1*n2) for the regular lattice.
//   Sites ordered column-major 

// [[Rcpp::export]]
NumericMatrix build_true_Sigma(NumericVector nvec,
                               double beta,
                               NumericVector lambdas,
                               NumericVector alphas,
                               double sigma,
                               std::string model = "modexp",
                               double delta = 1.0) {

  NumericMatrix Gamma = build_true_Gamma(nvec, beta, lambdas, alphas, sigma,
                                         model, delta);

  int n1 = (int)nvec[0], n2 = (int)nvec[1];
  int N  = n1 * n2;

  NumericMatrix Sigma(N, N);

  for (int s2 = 0; s2 < n2; ++s2) {
    for (int s1 = 0; s1 < n1; ++s1) {
      int s_idx = s2 * n1 + s1;
      for (int t2 = 0; t2 < n2; ++t2) {
        for (int t1 = 0; t1 < n1; ++t1) {
          int t_idx = t2 * n1 + t1;
          Sigma(s_idx, t_idx) = Gamma((s1 - t1) + (n1 - 1),
                                      (s2 - t2) + (n2 - 1));
        }
      }
    }
  }

  return Sigma;
}
