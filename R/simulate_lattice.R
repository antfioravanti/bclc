# =============================================================================
# Author: Antonio Fioravanti
# =============================================================================


# --------------------------------------------------------------
# GENERATE A REGULAR LATTICE GRID
#' Create a regular rectangular lattice grid.
#'
#' @param n1  Number of rows (vertical dimension).
#' @param n2  Number of columns (horizontal dimension). Defaults to n1.
#' @return  A data.frame with columns t1, t2 (integer coordinates).
generate_grid <- function(n1, n2 = n1) {
  expand.grid(t1 = seq_len(n1), t2 = seq_len(n2))
}

# -----------------------------------------------------------------------------
# SIMULATION OF GAUSSIAN PROCESSES ON THE LATTICE

#' Simulate a single realisation of a zero-mean Gaussian process on an
#' n1 x n2 lattice using the modified exponential covariance model.
#'
#' @param n1, n2   Grid dimensions.
#' @param sigma    Variance parameter (C(0,0) = sigma^2).
#' @param alpha1, alpha2  Smoothness parameters (positive).
#' @param lambda1, lambda2  Scale parameters (positive).
#' @param beta     Separability parameter (0 = separable).
#' @param seed     Random seed for reproducibility (NULL = no seed).
#' @return  A list with components:
#'   \item{X}{n1 x n2 matrix of simulated values.}
#'   \item{grid}{Data frame of grid coordinates.}
#'   \item{Sigma}{N x N true covariance matrix used for simulation.}
#'   \item{params}{Named list of model parameters.}
simulate_lattice_process <- function(n1, n2 = n1,
                                     sigma   = 1,
                                     alpha1  = 1, alpha2  = 1,
                                     lambda1 = 1, lambda2 = 1,
                                     beta    = 0,
                                     seed    = NULL) {

  grid  <- generate_grid(n1, n2)
  Sigma <- ModifiedExponentialCovariance(grid, sigma = sigma,
                                         alpha1 = alpha1, alpha2 = alpha2,
                                         lambda1 = lambda1, lambda2 = lambda2,
                                         beta = beta)

  if (!is.null(seed)) set.seed(seed)
  z <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = Sigma)
  X <- matrix(z, nrow = n1, ncol = n2)

  params <- list(n1 = n1, n2 = n2, sigma = sigma,
                 alpha1 = alpha1, alpha2 = alpha2,
                 lambda1 = lambda1, lambda2 = lambda2,
                 beta = beta)

  list(X = X, grid = grid, Sigma = Sigma, params = params)
}
