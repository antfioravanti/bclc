#------------------------------------------------------------------------------
# Author: Antonio Fioravanti
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# GENERATE A REGULAR LATTICE GRID
#' Create a regular rectangular lattice grid.
#'
#' @param n1  Number of rows (vertical dimension).
#' @param n2  Number of columns (horizontal dimension). Defaults to n1.
#' @return  A data.frame with columns t1, t2 (integer coordinates).
generate_grid <- function(n1, n2 = n1) {
  expand.grid(t1 = seq_len(n1), t2 = seq_len(n2))
}

#------------------------------------------------------------------------------
# SIMULATION OF GAUSSIAN PROCESSES ON THE LATTICE

#' Simulate a single realisation of a zero-mean Gaussian process on an
#' n1 x n2 lattice.
#'
#' @param n1,n2   Grid dimensions.
#' @param model    Covariance model: "modexp" (default), "gencauchy",
#'   "modexp2", "gencauchy2" or "gneiting" (see src/covariance.cpp). In the
#'   last three the margins do not depend on beta.
#' @param sigma    Variance parameter (C(0,0) = sigma^2).
#' @param alpha1,alpha2  Smoothness parameters, in (0, 2]. alpha2 defaults to
#'   alpha1.
#' @param lambda1,lambda2  Range parameters (positive). 
#'   lambda2 defaults to lambda1.
#' @param beta     Separability parameter: beta = 0 gives the separable
#'   model, beta > 0 the non-separable one. Must be in [0, 1] for "modexp2",
#'   "gencauchy2" and "gneiting".
#' @param delta    Tail parameter (> 0), used by "gencauchy", "gencauchy2"
#'   (requires alpha1 == alpha2; long memory for every beta iff delta <= 1)
#'   and "gneiting" (tau; requires delta >= beta/2).
#' @param seed     Random seed for reproducibility (NULL = no seed).
#' @param tol      Tolerance passed to MASS::mvrnorm.
#' @return  A list with components:
#'   \item{X}{n1 x n2 matrix of simulated values.}
#'   \item{grid}{Data frame of grid coordinates.}
#'   \item{Sigma}{N x N true covariance matrix used for simulation.}
#'   \item{params}{Named list of model parameters.}
simulate_lattice_process <- function(n1, n2 = n1, model   = "modexp",
                                     sigma   = 1,
                                     alpha1  = 1, alpha2  = alpha1,
                                     lambda1 = 1, lambda2 = lambda1,
                                     beta    = 0,
                                     delta   = 1,
                                     seed    = NULL,
                                     tol     = 1e-6) {

  grid  <- generate_grid(n1, n2)
  Sigma <- build_true_Sigma(c(n1, n2),
                            beta    = beta,
                            lambdas = c(lambda1, lambda2),
                            alphas  = c(alpha1, alpha2),
                            sigma   = sigma,
                            model   = model,
                            delta   = delta)

  if (!is.null(seed)) set.seed(seed)
  z <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = Sigma, tol = tol)
  X <- matrix(z, nrow = n1, ncol = n2)

  params <- list(n1 = n1, n2 = n2, sigma = sigma,
                 alpha1 = alpha1, alpha2 = alpha2,
                 lambda1 = lambda1, lambda2 = lambda2,
                 beta = beta,
                 model = model,
                 delta = delta)

  list(X = X, grid = grid, Sigma = Sigma, params = params)
}
