//-----------------------------------------------------------------------------
// Author: Antonio Fioravanti
// Data-driven selection of the truncation cut-off m for the bias-corrected
// lattice covariance estimator.
//
// Convention: h1 = vertical lag (rows), h2 = horizontal lag (columns).
//-----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include "naive.h"
#include "bc_estimator.h"
#include "select_m.h"

using namespace Rcpp;

//-----------------------------------------------------------------------------
// WLS FIT OF THE RANGE FROM A CORRELOGRAM
//-----------------------------------------------------------------------------
// Fit lambda in the exponential correlogram rho(h) = exp(-h / lambda) by
// weighted least squares through the origin on the log scale:
//
//   log rho(h) = -(1 / lambda) * h,      h = 1, ..., K
//   1 / lambda_hat = - sum_h w_h * h * log rho(h) / sum_h w_h * h^2
//
// with weights w_h proportional to the number of contributing pairs along the
// axis, w_h = n_k - h, so that short (better estimated) lags dominate. Lags with
// a non-positive or non-finite correlation are dropped (log undefined); if none
// survives, or the fitted slope is non-positive, NA is returned.
//
// Parameters:
//   rho: correlogram values rho(1), ..., rho(K).
//   n_k: axis length (only used to build the weights).
static double fit_lambda_wls(const std::vector<double>& rho, int n_k) {

  double num = 0.0, den = 0.0;
  int used = 0;

  for (std::size_t idx = 0; idx < rho.size(); ++idx) {
    // we need lags h = 1, ..., K, but idx is 0-based
    int    h = static_cast<int>(idx) + 1;
    double r = rho[idx];
    if (!std::isfinite(r) || r <= 0.0) continue;      // log undefined
    // weight w_h = n_k - h, but clamp to at least 1 to avoid zero weight
    double w = static_cast<double>(std::max(1, n_k - h));
    num += w * h * std::log(r);
    den += w * h * h;
    ++used;
  }
  // if no usable lags or the fitted slope is non-positive, return NA
  if (used == 0 || den <= 0.0) return NA_REAL;
  // fitted slope = -1 / lambda_hat
  double inv_lambda = -num / den;   
  // if the fitted slope is non-positive, return NA
  if (!std::isfinite(inv_lambda) || inv_lambda <= 0.0) return NA_REAL;
  
  // return the lambda hat
  return 1.0 / inv_lambda;
}

//-----------------------------------------------------------------------------
// BLOCK-HOLDOUT KRIGING CROSS-VALIDATION  (method = "cv")
//-----------------------------------------------------------------------------
// The cut-off is scored by how well the covariance it produces PREDICTS. For a
// candidate m: estimate Sigma_hat_m on the whole lattice, then hide a square
// block of sites, ordinary-krige them from the surrounding ring, and take the
// RMSE. Averaging over a few shifted blocks gives the CV score; argmin wins.
//
// The folds are the central block and four copies of it shifted by +-`shift`
// along both axes, so the scored blocks sit in slightly different parts of the
// field without ever leaving the lattice.
//
// Cost: each fold solves a dense (N - b^2) system, N = n1*n2, so the CV is
// O(|grid| * N^3) -- much heavier than the plug-in rule. Shrink cv_block or the
// grid on large lattices.

// One fold: the held-out block with top-left corner (r1, c1) and the ring around
// it, as 0-based column-major site indices into the n1 x n2 lattice.

// define the data structure named CvFold that contains two members: 
// pred and obs, both of type arma::uvec (Armadillo’s vector of unsigned 
// integer indices.)
struct CvFold { arma::uvec pred, obs; };

// ----
static CvFold make_fold(int n1, int n2, int r1, int c1, int b1, int b2) {
  // Splits a n1 x n2 lattice into a held-out block (pred) and the surrounding 
  // ring (obs), returning the column-major indices of the two sets.
   // r1, c1: top-left corner of the held-out block (0-based)
   // b1, b2: block size along each axis (height and width)
   // the held out block is the rectangle [r1, r1 + b1) x [c1, c1 + b2)
  std::vector<arma::uword> pred, obs;
  for (int c = 0; c < n2; ++c) {
    for (int r = 0; r < n1; ++r) {
      // convert the 2D index (r, c) into a 1D index in column-major order
      // i.e. indx(r, c) = c * n1 + r
      arma::uword idx = static_cast<arma::uword>(c * n1 + r);
      // if the location lies inside the selected block, add it to the 
      // pred vector, otherwise add it to the obs vector
      bool inside = (r >= r1 && r < r1 + b1 && c >= c1 && c < c1 + b2);
      // both pred and obs are  arma uword vectors, 
      // so we can use the push_back method to add the index to the
      (inside ? pred : obs).push_back(idx);
    }
  }

  CvFold f;
  f.pred = arma::uvec(pred);
  f.obs  = arma::uvec(obs);
  return f;
}
// ----


static std::vector<CvFold> build_cv_folds(int n1, int n2, int b, int shift) {
// It creates up to five square held-out blocks for cross-validation
    // n1, n2: lattice dimensions
    // b: block size (square) 
    // shift: diagonal offset of the four non-central folds

  // find the top left corner of the central block
  int s1 = (n1 - b) / 2, s2 = (n2 - b) / 2;
  // dr, dc: the diagonal shifts of the four directions
  const int dr[5] = {0, -1, -1, 1, 1};
  const int dc[5] = {0, -1, 1, -1, 1};

  std::vector<CvFold> folds;
  std::vector<std::pair<int,int> > seen;

  for (int k = 0; k < 5; ++k) {
    int r1 = std::min(std::max(s1 + dr[k] * shift, 0), n1 - b);
    int c1 = std::min(std::max(s2 + dc[k] * shift, 0), n2 - b);
    std::pair<int,int> corner(r1, c1);
    if (std::find(seen.begin(), seen.end(), corner) != seen.end()) continue;
    seen.push_back(corner);
    folds.push_back(make_fold(n1, n2, r1, c1, b, b));
  }
  return folds;
}

// Ordinary kriging of the held-out block from the ring, and its RMSE.
static double fold_rmse(const arma::mat& Sigma, const arma::vec& y,
                        const CvFold& f, double nugget) {

  arma::mat S_oo = Sigma(f.obs, f.obs);
  arma::mat S_po = Sigma(f.pred, f.obs);
  arma::vec y_obs = y(f.obs);
  int n_o = static_cast<int>(f.obs.n_elem);

  S_oo.diag() += nugget * arma::mean(Sigma.diag());   // jitter for stability

  // Solve once for both right-hand sides: the data and the constant 1.
  arma::mat rhs = arma::join_rows(y_obs, arma::ones<arma::vec>(n_o));
  arma::mat sol;
  if (!arma::solve(sol, S_oo, rhs)) return NA_REAL;

  double mu = arma::accu(sol.col(0)) / arma::accu(sol.col(1));   // OK mean
  arma::vec yhat = S_po * (sol.col(0) - mu * sol.col(1)) + mu;

  return std::sqrt(arma::mean(arma::square(yhat - y(f.pred))));
}

// CV score of one candidate cut-off: mean fold RMSE of Sigma_hat_(m1,m2).
static double cv_score_m(NumericMatrix X, int m1, int m2,
                         const std::vector<CvFold>& folds,
                         const arma::vec& y, double nugget, bool Nh_normalize) {

  int n1 = X.nrow(), n2 = X.ncol();
  if (m1 > n1 - 2 || m2 > n2 - 2) return NA_REAL;   // W*_m not admissible

  List fit   = bias_corrected_estimator(X, m1, m2, Nh_normalize);
  List Sig   = fit_to_Sigma(fit, "corrected");
  arma::mat Sigma = as<arma::mat>(Sig["Sigma_hat_m_O_psd"]);

  double total = 0.0; int used = 0;
  for (std::size_t k = 0; k < folds.size(); ++k) {
    double r = fold_rmse(Sigma, y, folds[k], nugget);
    if (std::isfinite(r)) { total += r; ++used; }
  }
  return used == 0 ? NA_REAL : total / used;
}

//-----------------------------------------------------------------------------
// DATA-DRIVEN SELECTION OF THE TRUNCATION CUT-OFF
//-----------------------------------------------------------------------------
// Select the truncation cut-off from the naive correlograms along the two
// coordinate axes and the grid size:
//
//   c0        = Chat(0,0)
//   rho_k(h)  = Chat(h * e_k) / c0        (e_1 = (1,0), e_2 = (0,1)); k = 1, 2
//   lambda_hat_k = WLS fit of rho_k(h) = exp(-h / lambda) over h = 1, ..., K
//   m_k       = round( kappa * (n_k * lambda_hat_k)^eta )
//
// With pooled = true a single cut-off is produced from the averaged correlogram
// rho(h) = (rho_1(h) + rho_2(h)) / 2 and n = min(n1, n2), and m1 = m2 = m.
//
// method = "cv" ignores all of the above and instead picks the cut-off that
// minimises the block-holdout kriging CV score (see the section above). With
// pooled = true only m1 = m2 is scanned (|grid| fits); with pooled = false the
// full m1 x m2 grid is scanned (|grid|^2 fits, and each fit is O(N^3)).
//
// Parameters:
//   X:            n1 x n2 data matrix.
//   method:       "rule" -> the plug-in rule - fitted exponential correlogram;
//                 "cv"   -> block-holdout kriging CV.
//   pooled:       pool the two axes into a single cut-off m1 = m2 (see above).
//   kappa, eta:   tuning constants of the rule ("rule" only).
//   K:            number of correlogram lags used in the WLS fit ("rule" only);
//                 if < 0 (the default) K = max(1, floor(min(n1, n2) / 2)). Per
//                 axis K is capped at n_k - 1.
//   m_max:        if >= 0, cap each m_k at min(m_max, n_k - 2); if < 0, ignore.
//   Nh_normalize: passed through to naive_sample_cov().
//   cv_grid:      candidate cut-offs scanned by the CV; empty (the default)
//                 means 1 ... min(min(n1,n2) / 2, upper cap).
//   cv_block:     side of the held-out square block; < 0 (the default) means
//                 max(2, min(n1, n2) / 3).
//   cv_shift:     diagonal offset of the four non-central folds, in sites.
//   cv_nugget:    relative jitter added to the kriging system's diagonal.
// Returns a list recording how the cut-off was reached. Always: m1, m2, n1, n2
// (int), pooled (bool), method (string). For "rule" additionally rho1, rho2 (the
// lag-one correlations), rho1_curve, rho2_curve, lambda1_hat, lambda2_hat, K,
// kappa, eta. For "cv" additionally cv_grid, cv_score (a vector when pooled, a
// |grid| x |grid| matrix of m1 by m2 scores otherwise), cv_block, cv_shift.
// [[Rcpp::export]]
List select_m(NumericMatrix X,
              std::string method = "rule",
              bool pooled = false,
              double kappa = 0.445,
              double eta = 0.5,
              int K = -1,
              int m_max = -1,
              bool Nh_normalize = true,
              IntegerVector cv_grid = IntegerVector::create(),
              int cv_block = -1,
              int cv_shift = 2,
              double cv_nugget = 1e-6) {

  int n1 = X.nrow();
  int n2 = X.ncol();

  if (n1 < 3)
    stop("select_m: dimension 1 has size n1 = %d < 3; no admissible cut-off "
         "exists (need n1 - 2 >= 1).", n1);
  if (n2 < 3)
    stop("select_m: dimension 2 has size n2 = %d < 3; no admissible cut-off "
         "exists (need n2 - 2 >= 1).", n2);

  if (method != "rule" && method != "cv")
    stop("select_m: method must be \"rule\" (alias \"param\") or \"cv\".");

  //---------------------------------------------------------------------------
  // method = "cv": block-holdout kriging cross-validation
  //---------------------------------------------------------------------------
  if (method == "cv") {
    
    // upper limit of admissimble m cut-off 
    int upper = (m_max >= 0) ? std::min(m_max, std::min(n1, n2) - 2)
                             : std::min(n1, n2) - 2;

    // Default grid: 1 ... min(n1,n2)/2, within the admissible range.
    std::vector<int> grid;
    if (cv_grid.size() == 0) {
      for (int m = 1; m <= std::min(std::min(n1, n2) / 2, upper); ++m)
      // push_back is a method of std::vector that appends an element to the
      // end of the vector - does not modify an existing entry; 
      // it increases the vector length by one and appends
        grid.push_back(m); // add m at the end of the vector grid
    } else {
      for (int i = 0; i < cv_grid.size(); ++i)
        if (cv_grid[i] >= 1 && cv_grid[i] <= upper) grid.push_back(cv_grid[i]);
    }
    if (grid.empty())
      stop("select_m: the cross-validation grid is empty after clamping to "
           "[1, %d]; widen cv_grid or relax m_max.", upper);

    int b = (cv_block >= 0) ? cv_block : std::max(2, std::min(n1, n2) / 3);
    b = std::min(b, std::min(n1, n2) - 1);           // leave a ring to krige from
    if (b < 1) stop("select_m: cv_block must leave at least one ring site.");

    std::vector<CvFold> folds = build_cv_folds(n1, n2, b, cv_shift);
    arma::vec y = arma::vectorise(as<arma::mat>(X));  // column-major, as Sigma

    int G = static_cast<int>(grid.size());
    int m1 = grid[0], m2 = grid[0];
    double best = R_PosInf;
    NumericVector score_vec;   // pooled
    NumericMatrix score_mat;   // full grid

    if (pooled) {
      score_vec = NumericVector(G);
      for (int i = 0; i < G; ++i) {
        score_vec[i] = cv_score_m(X, grid[i], grid[i], folds, y, cv_nugget,
                                  Nh_normalize);
        if (std::isfinite(score_vec[i]) && score_vec[i] < best) {
          best = score_vec[i]; m1 = m2 = grid[i];
        }
      }
    } else {
      score_mat = NumericMatrix(G, G);               // rows = m1, cols = m2
      for (int j = 0; j < G; ++j)
        for (int i = 0; i < G; ++i) {
          double s = cv_score_m(X, grid[i], grid[j], folds, y, cv_nugget,
                                Nh_normalize);
          score_mat(i, j) = s;
          if (std::isfinite(s) && s < best) { best = s; m1 = grid[i]; m2 = grid[j]; }
        }
    }

    if (!std::isfinite(best))
      Rcpp::warning("select_m: no cross-validation score was finite; falling "
                    "back to the smallest cut-off in the grid.");

    // A vector when pooled, a |grid| x |grid| matrix otherwise; wrapped here so
    // both arms of the conditional have the same (SEXP) type.
    SEXP cv_score = pooled ? Rcpp::wrap(score_vec) : Rcpp::wrap(score_mat);

    return List::create(
      Named("m1")       = m1,
      Named("m2")       = m2,
      Named("cv_grid")  = Rcpp::wrap(grid),
      Named("cv_score") = cv_score,
      Named("cv_block") = b,
      Named("cv_shift") = cv_shift,
      Named("n_folds")  = static_cast<int>(folds.size()),
      Named("n1")       = n1,
      Named("n2")       = n2,
      Named("pooled")   = pooled,
      Named("method")   = "cv"
    );
  }

  //---------------------------------------------------------------------------
  // method = "rule": plug-in cut-off from the correlograms
  //---------------------------------------------------------------------------
  // Number of correlogram lags used in the WLS fit.
  if (K < 0) K = std::max(1, std::min(n1, n2) / 2);

  double c0 = naive_sample_cov(X, NumericVector::create(0, 0), Nh_normalize);

  // Correlograms along each axis, rho_k(1), ..., rho_k(K_k).
  auto correlogram = [&](int axis, int n_k) {
    int K_k = std::max(1, std::min(K, n_k - 1));
    std::vector<double> rho(K_k, NA_REAL);
    if (c0 <= 0.0 || !std::isfinite(c0)) return rho;
    for (int h = 1; h <= K_k; ++h) {
      NumericVector hvec = (axis == 1) ? NumericVector::create(h, 0)
                                       : NumericVector::create(0, h);
      // 0-based index, but hvec is 1-based (lag h = 1, ..., K_k)
      rho[h - 1] = naive_sample_cov(X, hvec, Nh_normalize) / c0;
    }
    return rho;
  };

  std::vector<double> rho1 = correlogram(1, n1);
  std::vector<double> rho2 = correlogram(2, n2);

  // kappa * (n * lambda_hat)^eta, clamped to the admissible range [1, n - 2].
  auto cutoff = [&kappa, &eta, &m_max](double lambda_hat, int n_k,
                                       const char* dim) -> int {
    if (!std::isfinite(lambda_hat)) {
      Rcpp::warning("select_m: the correlogram of dimension %s gives no usable "
                    "range estimate (c0 <= 0 or no positive correlation); "
                    "falling back to m = 1.", dim);
      return 1;
    }
    // kappa * (n_k * lambda_hat)^eta, rounded to the nearest integer
    int m = static_cast<int>(std::lround(kappa * std::pow(n_k * lambda_hat, eta)));
    int upper = (m_max >= 0) ? std::min(m_max, n_k - 2) : (n_k - 2);
    return std::max(1, std::min(m, upper));
  };

  int m1, m2;
  double lambda1_hat, lambda2_hat;

  if (pooled) {
    int n = std::min(n1, n2); // choose the smaller axis
    // choose the smaller correlogram
    std::size_t Kp = std::min(rho1.size(), rho2.size()); 
    std::vector<double> rho_pool(Kp);
    // average the two correlograms, ignoring the longer one
    for (std::size_t i = 0; i < Kp; ++i) rho_pool[i] = 0.5 * (rho1[i] + rho2[i]);

    double lambda_hat = fit_lambda_wls(rho_pool, n);
    int m = cutoff(lambda_hat, n, "pooled");
    m1 = m2 = m;
    lambda1_hat = lambda2_hat = lambda_hat;
  } else {
    lambda1_hat = fit_lambda_wls(rho1, n1);
    lambda2_hat = fit_lambda_wls(rho2, n2);
    m1 = cutoff(lambda1_hat, n1, "1");
    m2 = cutoff(lambda2_hat, n2, "2");
  }

  return List::create(
    Named("m1")          = m1,
    Named("m2")          = m2,
    Named("rho1")        = rho1.empty() ? NA_REAL : rho1[0],
    Named("rho2")        = rho2.empty() ? NA_REAL : rho2[0],
    Named("rho1_curve")  = Rcpp::wrap(rho1),
    Named("rho2_curve")  = Rcpp::wrap(rho2),
    Named("lambda1_hat") = lambda1_hat,
    Named("lambda2_hat") = lambda2_hat,
    Named("n1")          = n1,
    Named("n2")          = n2,
    Named("K")           = K,
    Named("kappa")       = kappa,
    Named("eta")         = eta,
    Named("pooled")      = pooled,
    Named("method")      = method
  );
}
