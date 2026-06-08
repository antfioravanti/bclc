#------------------------------------------------------------------------------
# Author: Antonio Fioravanti
# Basic utility functions for working with 2D lattice data
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# REVERSAL MATRIX
#' Build the \( n \times n \) reversal matrix \( \mathbf{J} \) (anti-identity).
#'
#' The matrix has ones on the anti-diagonal and zeros elsewhere.
#' Multiplication \( \mathbf{J} \%*\% \underline{x} \) reverses the entries of a vector \( \underline{x} \).
#'
#' @param n Integer. Dimension of the matrix.
#' @return An \( n \times n \) numeric matrix representing the reversal operator.
Jmat <- function(n) {
  diag(n)[, n:1, drop = FALSE]
}


#------------------------------------------------------------------------------
# POSITIVE SEMI DEFINITE CHECKER FOR A MATRIX
#' Check Whether a Matrix Is Positive Semi-Definite
#'
#' Returns `TRUE` when `M` is symmetric up to the specified tolerance and all
#' eigenvalues are greater than or equal to `-tol`.
#'
#' @param M A numeric matrix.
#' @param tol Numeric tolerance used for the symmetry check and eigenvalue test.
#'
#' @return A logical scalar.
#' @export
#------------------------------------------------------------------------------
# NEAREST PSD PROJECTION
#' Project a symmetric matrix to the nearest positive semi-definite matrix
#' (Frobenius-optimal). Floors all negative eigenvalues to zero.
#'
#' @param M A numeric symmetric matrix.
#' @return The nearest PSD matrix.
psd_project <- function(M) {
  M_sym <- (M + t(M)) / 2                  # symmetrise numerically
  ed    <- eigen(M_sym, symmetric = TRUE)
  lam   <- pmax(ed$values, 0)
  ed$vectors %*% (lam * t(ed$vectors))     # avoids diag(lam) allocation
}

#------------------------------------------------------------------------------
# POSITIVE SEMI DEFINITE CHECKER FOR A MATRIX
is_psd <- function(M, tol = 1e-8) {
  if (!isSymmetric(M, tol = tol)) return(FALSE)
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  all(ev >= -tol)
}

#------------------------------------------------------------------------------
# MATRIX TO DATA FRAME CONVERSION
#' Convert a matrix with dimnames "h1=..." / "h2=..." to a long data frame.
#' @param mat  Matrix with named rows and columns.
#' @param name_val  Name for the value column (default "val").
matrix_to_df <- function(mat, name_val = "val") {
  h1_vals <- as.numeric(gsub("h1=", "", rownames(mat)))
  h2_vals <- as.numeric(gsub("h2=", "", colnames(mat)))
  coords  <- expand.grid(h1 = h1_vals, h2 = h2_vals)
  df <- data.frame(h1 = coords$h1, h2 = coords$h2, val = as.vector(mat))
  names(df)[3] <- name_val
  df
}

#------------------------------------------------------------------------------
# GET LAG INDEX IN VEC(GAMMA)
#' Map a 2D lag vector (h1, h2) to its position in vec(Gamma).
#'
#' The full lag matrix Gamma has L1 = 2*n1-1 rows and L2 = 2*n2-1 columns.
#' vec() stacks columns (column-major order).
#'
#' @param nvec  Grid size c(n1, n2).
#' @param hvec  Lag vector c(h1, h2).
#' @return A 1-based integer index matching R's column-major vectorization,
#'   for example the position in `as.vector(Gamma)`.
get_lag_index <- function(nvec, hvec) {
  n1 <- as.integer(nvec[1])
  n2 <- as.integer(nvec[2])
  h1 <- as.integer(hvec[1])
  h2 <- as.integer(hvec[2])

  L1 <- 2L * n1 - 1L

  # Zero-based row and column position in the lag matrix
  row0 <- h1 + (n1 - 1L)
  col0 <- h2 + (n2 - 1L)

  if (row0 < 0 || row0 >= L1 || col0 < 0 || col0 >= (2L * n2 - 1L))
    stop(sprintf("Lag (%d,%d) is outside the range for grid (%d,%d).", h1, h2, n1, n2))

  # Column-major: index = col * nrow + row (0-based), then +1
  col0 * L1 + row0 + 1L
}


#------------------------------------------------------------------------------
# CENTRAL BLOCK INDICES FOR TRUNCATION
#' Get the vec-indices of the central (2m+1)^2 block within vec(Gamma).
#'
#' Used to extract the truncated lag window {-m,...,m}^2 from the full
#' L1*L2 vectorised lag matrix.
#'
#' @param nvec  Grid size c(n1, n2).
#' @param m     Truncation parameter.
#' @return  Integer vector of 1-based indices into vec(Gamma).
slice_indices <- function(nvec, m) {
  n1 <- nvec[1]; n2 <- nvec[2]
  L1 <- 2L * n1 - 1L
  L2 <- 2L * n2 - 1L

  # Row and column indices of the central block (1-based in the lag matrix)
  centre1 <- n1   # row index of h1 = 0
  centre2 <- n2   # col index of h2 = 0
  rows <- (centre1 - m):(centre1 + m)
  cols <- (centre2 - m):(centre2 + m)

  # Expand to column-major vec indices
  as.vector(outer(rows, cols, FUN = function(r, c) (c - 1L) * L1 + r))
}


#------------------------------------------------------------------------------
# LAG MATRIX PRINTER
#' Print a visual table of all lag pairs (h1, h2) and optionally their
#' vec-index, for a grid of size (n1, n2).
#'
#' This is a small helper for inspecting the lag layout of the full
#' \eqn{(2n_1 - 1) \times (2n_2 - 1)} lag matrix. When `show_index = TRUE`,
#' each cell also displays the corresponding 1-based index in R's
#' column-major vectorization.
#'
#' @param n1 Integer. Number of rows in the lattice.
#' @param n2 Integer. Number of columns in the lattice.
#' @param show_index Logical. If `TRUE`, print the associated 1-based
#'   vectorization index together with each lag pair.
#'
#' @return Invisibly returns a character matrix containing the printed lag
#'   labels.
visualise_lags <- function(n1, n2, show_index = FALSE) {
  h1 <- seq(-(n1 - 1), n1 - 1)
  h2 <- seq(-(n2 - 1), n2 - 1)
  L1 <- length(h1)

  cell <- outer(h1, h2, FUN = function(x, y) {
    if (show_index) {
      idx <- (y + n2 - 1L) * L1 + (x + n1 - 1L) + 1L
      sprintf("(%d,%d)[%d]", x, y, idx)
    } else {
      sprintf("(%d,%d)", x, y)
    }
  })
  dimnames(cell) <- list(h1 = h1, h2 = h2)
  print(cell, quote = FALSE)
  invisible(cell)
}
