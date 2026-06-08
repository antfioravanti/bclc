#------------------------------------------------------------------------------
# Author: Antonio Fioravanti
# Bias-corrected sample covariance estimator on a 2D rectangular lattice.
#------------------------------------------------------------------------------

# files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
# lapply(files, source)

#------------------------------------------------------------------------------
# BUILD WEIGHT MATRIX W_n
#' Build the full and normalised weight matrices W_n and W*_n from the
#' one-dimensional count matrices for each spatial dimension.
#'
#' @param mat1  Output of build_count_matrices(n1).
#' @param mat2  Output of build_count_matrices(n2).
#' @return  A list with:
#'   \item{W}{L1*L2 x L1*L2 unnormalised weight matrix (I - A(x)A - B(x)B + D(x)D).}
#'   \item{W_star}{L1*L2 x L1*L2 normalised weight matrix.}
build_W <- function(mat1, mat2) {

  p <- mat1$L * mat2$L

  # Unnormalised
  AkA <- kronecker(mat2$A, mat1$A)
  BkB <- kronecker(mat2$B, mat1$B)
  DkD <- kronecker(mat2$D, mat1$D)
  W   <- diag(p) - AkA - BkB + DkD

  # Normalised
  AkA_n <- kronecker(mat2$A_norm, mat1$A_norm)
  BkB_n <- kronecker(mat2$B_norm, mat1$B_norm)
  DkD_n <- kronecker(mat2$D_norm, mat1$D_norm)
  W_star <- diag(p) - AkA_n - BkB_n + DkD_n

  list(W = W, W_star = W_star)
}

#------------------------------------------------------------------------------
# Helper: extracts the central (2m+1) block of an L x L matrix.
# central_range <- function(L, m) {
#   centre <- (L + 1L) %/% 2L   # alternatie to (ceiling(L/2)-m):(ceiling(L/2)+m)
#   (centre - m):(centre + m)
# }

# # BUILD TRUNCATED WEIGHT MATRIX W_{n,m}
# #' Build the (2m+1)^2 x (2m+1)^2 reduced weight matrix W_{n,m} for a given
# #' cut-off value m, by extracting central submatrices from the
# #' one-dimensional count matrices before taking Kronecker products.
# #'
# #' This is more efficient than building the full W*_n and then subsetting,
# #' because the Kronecker product is taken on (2m+1) x (2m+1) blocks
# #' instead of L x L blocks.
# #'
# #' @param mat1  Output of build_count_matrices(n1).
# #' @param mat2  Output of build_count_matrices(n2).
# #' @param m     Truncation parameter (lags in {-m,...,m}^2).
# #' @return  A list with W_m and W_star_m.
# build_W_m <- function(mat1, mat2, m) {

#   # Central (2m+1) rows/cols of the L_i x L_i matrices
#   idx1 <- central_range(mat1$L, m)
#   idx2 <- central_range(mat2$L, m)

#   krondim <- (2L * m + 1L)^2

#   # Unnormalised
#   A1 <- mat1$A[idx1, idx1];  A2 <- mat2$A[idx2, idx2]
#   B1 <- mat1$B[idx1, idx1];  B2 <- mat2$B[idx2, idx2]
#   D1 <- mat1$D[idx1, idx1];  D2 <- mat2$D[idx2, idx2]
#   W_m <- diag(p) - kronecker(A2, A1) - kronecker(B2, B1) + kronecker(D2, D1)

#   # Normalised
#   A1n <- mat1$A_norm[idx1, idx1];  A2n <- mat2$A_norm[idx2, idx2]
#   B1n <- mat1$B_norm[idx1, idx1];  B2n <- mat2$B_norm[idx2, idx2]
#   D1n <- mat1$D_norm[idx1, idx1];  D2n <- mat2$D_norm[idx2, idx2]
#   W_star_m <- diag(krondim) - kronecker(A2n, A1n) - kronecker(B2n, B1n) + kronecker(D2n, D1n)

#   list(W_m = W_m, W_star_m = W_star_m)
# }

# build_Pm, build_Sm, build_W_m, bias_corrected_estimator live in src/estimators.cpp
#------------------------------------------------------------------------------
# BUILD WEIGHT MATRIX W_{n,m} from grid dimensions and m
#' Build both W_m and W*_m directly from grid dimensions and m.
#' Thin wrapper around build_count_matrices + build_W_m.
#'
#' @param n1, n2  Grid dimensions.
#' @param m  Truncation parameter.
#' @param Nh_normalize  Normalisation for the counting matrices.
#' @return  Same as build_W_m().
build_weight_matrices <- function(n1, n2, m, Nh_normalize = TRUE) {
  mat1 <- build_count_matrices(n1, Nh_normalize)
  mat2 <- build_count_matrices(n2, Nh_normalize)
  build_W_m(mat1, mat2, m)
}
#------------------------------------------------------------------------------
# CHECK INVERTIBILITY OF W*_m
#' For each m from 1 to m_max, check whether W*_{n,m} is singular.
#' @param mat1, mat2  Outputs of build_count_matrices for dimensions 1, 2.
#' @param m_max  Maximum truncation to test (default: min(n1, n2) - 2).
#' @param tol  Threshold for rcond below which the matrix is declared singular.
#' @return  A data.frame with columns m, rcond, singular.
check_invertibility <- function(mat1, mat2, m_max = NULL,
                                tol = .Machine$double.eps) {

  n1 <- (mat1$L + 1L) %/% 2L
  n2 <- (mat2$L + 1L) %/% 2L
  if (is.null(m_max)) m_max <- min(n1, n2) - 2L

  results <- data.frame(m = integer(m_max), rcond = numeric(m_max),
                         singular = logical(m_max))

  for (m in seq_len(m_max)) {
    Wm <- build_W_m(mat1, mat2, m)$W_star_m
    rc <- rcond(Wm)
    results$m[m]        <- m
    results$rcond[m]    <- rc
    results$singular[m] <- (rc < tol)
  }

  return(results)
}