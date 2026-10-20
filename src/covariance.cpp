//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// True autocovariance models, selected by
//   model = "modexp" | "gencauchy"                       
//         | "modexp2" | "gencauchy2" | "gneiting"  
// Reference: https://rdrr.io/cran/geoR/man/cov.spatial.html
//
//   d1  = (|h1| / lambda1)^alpha1,   d2 = (|h2| / lambda2)^alpha2
//
//   d12 = (|h1 - h2| / sqrt(lambda1 * lambda2))^min(alpha1, alpha2)
//   modexp     beta = 0:  C = sigma^2 * exp(-d1) * exp(-d2)
//              beta > 0:  C = sigma^2 * exp(-(d1 + d2 + beta * d12))
//   gencauchy  beta = 0:  C = sigma^2 * (1 + d1)^(-delta/alpha1) * (1 + d2)^(-delta/alpha2)
//              beta > 0:  C = sigma^2 * (1 + d1 + d2 + beta * d12)^(-delta/ahar)
//   NOTE: in these models beta also changes the margins and breaks axial
//   symmetry, so beta is NOT a pure separability parameter.
//
//  NEW MODELS:
//   Margins C(h1,0), C(0,h2) do not depend on beta; beta = 0 is separable.
//   Dq  = (d1^q + d2^q)^(1/q),  q = 2 / max(alpha1, alpha2)
//   modexp2     C = sigma^2 * exp(-[(1-beta)(d1 + d2) + beta * Dq])
//   gencauchy2  C = sigma^2 * [(1+d1)(1+d2)]^(-(1-beta) nu) * (1+Dq)^(-beta nu),
//                  nu = delta/alpha, requires alpha1 == alpha2;
//                  short memory iff delta*(2-beta) > 2 (beta couples memory);
//                  long memory for all beta in [0,1] iff delta <= 1
//   gneiting       Gneiting (2002), adapted to spatial only
//                  C = sigma^2 * (1 + d2)^(-delta) * exp(-d1 / (1 + d2)^(beta*alpha1/2))
//                  adapted such that Gneiting's:
//                  c = lambda1^(-alpha1)
//                  a = lambda2^(-alpha2),
//                  tau = delta.
//                  alpha_G = alpha2/2,
//                  gamma_G = alpha1/2, 
//                  Requires:
//                   delta > 0 and delta >= beta/2;
//                  short memory iff alpha2*(delta-beta/2) > 1.
//                  Calibrated to old modexp (lambda): 
//                  alphas = (1, 2),  lambdas = (lambda, lambda), delta = 1.5.
//
//   where:
//   lambda: range parameter (0, inf)
//   alpha:  shape parameter (0, 2]
//   beta:   separability parameter (old: [0, inf), new: [0, 1]); 0 = separable
//   delta:  gencauchy: tail exponent (0, inf); gneiting: tau
//-----------------------------------------------------------------------------

#include <Rcpp.h>
#include <cmath>
#include <cctype>      // std::tolower (needed for portable package builds)
#include <string>
#include <algorithm>
using namespace Rcpp;

//-----------------------------------------------------------------------------
// Model id and parameter bundle
//-----------------------------------------------------------------------------
enum CovModel {
  COV_MODEXP        = 0,   // legacy
  COV_GENCAUCHY     = 1,   // legacy (not used in paper)
  COV_MODEXP2       = 2,   // NEW v0.8.0
  COV_GENCAUCHY2    = 3,   // NEW v0.8.0
  COV_GNEITING      = 4    // NEW v0.8.0
};

struct CovParams {
  double beta      = 0.0;
  double lambda1   = 1.0, lambda2 = 1.0;
  double alpha1    = 1.0, alpha2  = 1.0;
  double sigma     = 1.0;
  double delta     = 1.0;   // gencauchy*: tail exponent; gneiting: tau
  int    model     = COV_MODEXP;
};

static int parse_cov_model(std::string s) {
  // convert string into lower case
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (s == "modexp")        return COV_MODEXP;
  if (s == "gencauchy")     return COV_GENCAUCHY;
  if (s == "modexp2")    return COV_MODEXP2;      // NEW
  if (s == "gencauchy2") return COV_GENCAUCHY2;   // NEW
  if (s == "gneiting")      return COV_GNEITING;        // NEW
  stop("unknown covariance model '" + s + "'; expected \"modexp\", "
       "\"gencauchy\", \"modexp2\", \"gencauchy2\" or \"gneiting\"");
  return COV_MODEXP; // never reached
}

static void validate_cov_params(const CovParams& p) {
  if (p.lambda1 <= 0.0 || p.lambda2 <= 0.0)
    stop("lambdas must be > 0");
  if (p.alpha1 <= 0.0 || p.alpha1 > 2.0 || p.alpha2 <= 0.0 || p.alpha2 > 2.0)
    stop("alphas must be in (0, 2] for positive-definiteness");
  if (p.beta < 0.0)
    stop("beta must be >= 0");
  if ((p.model == COV_GENCAUCHY || p.model == COV_GENCAUCHY2 ||
       p.model == COV_GNEITING) && p.delta <= 0.0)
    stop("delta must be > 0 for the gencauchy and gneiting models");

  // ---- NEW: checks for the paper-2 models ----
  const bool new_model = (p.model == COV_MODEXP2 ||
                          p.model == COV_GENCAUCHY2 ||
                          p.model == COV_GNEITING);
  if (new_model && p.beta > 1.0)
    stop("beta must be in [0, 1] for modexp2, gencauchy2 and gneiting");

  if ((p.model == COV_MODEXP2 || p.model == COV_GENCAUCHY2) &&
      p.beta > 0.0 && std::max(p.alpha1, p.alpha2) >= 2.0)
    Rcpp::warning("max(alpha) = 2: q = 1, so Dq = d1 + d2 and beta has no "
                  "effect (the model reduces to the separable one)");

  if (p.model == COV_GENCAUCHY2) {
    if (p.alpha1 != p.alpha2)
      stop("gencauchy2 requires alpha1 == alpha2");
    if (p.delta > 1.0)
      Rcpp::warning("gencauchy2: beta also couples the memory: short memory iff "
                    "delta*(2-beta) > 2. Use delta <= 1 for long memory at every "
                    "beta in [0, 1] (separable and non-separable alike)");
  }

  if (p.model == COV_GNEITING) {
    if (p.delta < p.beta / 2.0)
      stop("gneiting requires tau (= delta) >= beta/2");
    if (p.alpha2 * (p.delta - p.beta / 2.0) <= 1.0)
      Rcpp::warning("gneiting: alpha2*(tau - beta/2) <= 1, covariance not "
                    "summable (long memory); e.g. use alphas = (1, 2), delta = 1.5");
  }
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

  // ---------------- paper 1 models (unchanged) ----------------
  if (p.model == COV_MODEXP || p.model == COV_GENCAUCHY) {
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
    } else {
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
  }

  // ---------------- NEW: paper 2 models ----------------
  if (p.model == COV_MODEXP2 || p.model == COV_GENCAUCHY2) {
    // non-additive variogram with the same margins as d1 + d2
    const double q  = 2.0 / std::max(p.alpha1, p.alpha2);          // q >= 1
    const double Dq = std::pow(std::pow(d1, q) + std::pow(d2, q), 1.0 / q);

    if (p.model == COV_MODEXP2) {
      return sigma2 * std::exp(-((1.0 - p.beta) * (d1 + d2) + p.beta * Dq));
    } else {
      const double nu = p.delta / p.alpha1;                          // alpha1 == alpha2
      return sigma2 * std::pow((1.0 + d1) * (1.0 + d2), -(1.0 - p.beta) * nu)
                    * std::pow(1.0 + Dq, -p.beta * nu);
    }
  }

  if (p.model == COV_GNEITING) {
    // Gneiting (2002): psi = 1 + a|h2|^(2 alpha_G) = 1 + d2,
    // c|h1|^(2 gamma_G) = d1, exponent beta*gamma_G = beta*alpha1/2, tau = delta
    const double psi = 1.0 + d2;
    return sigma2 * std::pow(psi, -p.delta)
                  * std::exp(-d1 / std::pow(psi, 0.5 * p.beta * p.alpha1));
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
// model selected by `model` (see header).
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