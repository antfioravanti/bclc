#------------------------------------------------------------------------------
# mse.R
# R-level wrappers and utilities for the analytical (true-parameter) MSE
# computations implemented in true_cov_precompute.cpp.
#
# The C++ functions build_Omega_m(), build_bias_m(), build_MSE_m(), etc.
# are called directly from R. This file provides:
#   - A unified MSE computation pipeline for the bias-corrected estimator
#   - Scalar MSE (trace) extraction
#   - Comparison helpers between standard and compact formulas
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# MSE of the bias-corrected estimator C*_m

#' Compute the MSE matrix of the bias-corrected estimator C*_m.
#'
#' The MSE of the standard estimator is
#'   MSE(Chat_m) = Omega_m + b_m b_m'
#' where Omega_m = Cov(vec(Gamma_hat_m)) and b_m is the bias vector.
#'
#' The bias-corrected estimator transforms this as:
#'   MSE(C*_m) = (W*_m)^{-1} Omega_m (W*_m)^{-T} + b*_m (b*_m)'
#' where b*_m = (W*_m)^{-1} b_m.
#'
#' @param n1, n2    Grid dimensions.
#' @param m         Truncation parameter.
#' @param beta, lambdas, alphas, sigma  Model parameters.
#' @param mu        True mean (only needed for 16-term formula; ignored by compact).
#' @param compact   If TRUE, use the mu-free compact formula for Omega.
#' @param Nh_normalize  Normalisation for counting matrices.
#' @return  A list with:
#'   \item{MSE_standard}{MSE matrix of the standard estimator Chat_m.}
#'   \item{MSE_corrected}{MSE matrix of the bias-corrected estimator C*_m.}
#'   \item{Omega_m}{Covariance matrix of vec(Gamma_hat_m).}
#'   \item{bias_m}{Bias vector of the standard estimator.}
#'   \item{bias_corrected}{Bias vector of the corrected estimator.}
#'   \item{W_star_m}{Weight matrix used for correction.}
#'   \item{scalar_MSE_standard}{tr(MSE_standard).}
#'   \item{scalar_MSE_corrected}{tr(MSE_corrected).}
compute_mse_comparison <- function(n1, n2, m,
                                    beta, lambdas, alphas, sigma,
                                    mu = 0,
                                    compact = TRUE,
                                    Nh_normalize = TRUE) {

  # Build Omega_m (covariance of the sample autocovariance vector)
  if (compact) {
    Omega_m <- build_Omega_m_compact(n1, n2, m, beta, lambdas, alphas, sigma)
  } else {
    Omega_m <- build_Omega_m(n1, n2, m, mu, beta, lambdas, alphas, sigma)
  }

  # Bias vector of the standard estimator
  bias_m <- build_bias_m(n1, n2, m, beta, lambdas, alphas, sigma)

  # MSE of the standard estimator: Omega + b b'
  MSE_std <- Omega_m + outer(bias_m, bias_m)

  # Build weight matrix W*_m
  mat1 <- build_count_matrices(n1, Nh_normalize)
  mat2 <- build_count_matrices(n2, Nh_normalize)
  Wm   <- build_W_m(mat1, mat2, m)
  W_star_m <- Wm$W_star_m

  # Invert W*_m
  W_inv <- solve(W_star_m)

  # MSE of the corrected estimator
  Omega_corrected <- W_inv %*% Omega_m %*% t(W_inv)
  bias_corrected  <- as.vector(W_inv %*% bias_m)
  MSE_corr <- Omega_corrected + outer(bias_corrected, bias_corrected)

  list(
    MSE_standard       = MSE_std,
    MSE_corrected      = MSE_corr,
    Omega_m            = Omega_m,
    bias_m             = bias_m,
    bias_corrected     = bias_corrected,
    W_star_m           = W_star_m,
    scalar_MSE_standard  = sum(diag(MSE_std)),
    scalar_MSE_corrected = sum(diag(MSE_corr))
  )
}


#------------------------------------------------------------------------------
# Scalar MSE extraction

#' Extract the scalar (trace) MSE from an MSE matrix.
#' Equivalent to sum of diagonal = sum of individual MSE(Chat(h)) over all h.
scalar_mse <- function(MSE_matrix) {
  sum(diag(MSE_matrix))
}


#------------------------------------------------------------------------------
# MSE ratio for a single (n, m) combination

#' Compute the ratio tr(MSE_corrected) / tr(MSE_standard).
#' Values < 1 indicate the bias-corrected estimator is better.
mse_ratio <- function(n1, n2, m, beta, lambdas, alphas, sigma,
                      mu = 0, compact = TRUE, Nh_normalize = TRUE) {

  res <- compute_mse_comparison(n1, n2, m, beta, lambdas, alphas, sigma,
                                 mu = mu, compact = compact,
                                 Nh_normalize = Nh_normalize)
  res$scalar_MSE_corrected / res$scalar_MSE_standard
}


#------------------------------------------------------------------------------
# Verify compact vs 16-term formula

#' Check that the compact (mu-free) and 16-term formulas produce the same
#' Omega matrix, for a given set of parameters. Returns the max absolute
#' difference.
verify_compact_vs_16term <- function(n1, n2, m, mu, beta, lambdas, alphas, sigma) {
  Omega_16  <- build_Omega_m(n1, n2, m, mu, beta, lambdas, alphas, sigma)
  Omega_cpt <- build_Omega_m_compact(n1, n2, m, beta, lambdas, alphas, sigma)
  max(abs(Omega_16 - Omega_cpt))
}
