#------------------------------------------------------------------------------
# Author: Antonio Fioravanti
# Diagnostic functions for checking structural, spectral, and algebraic
# properties of the counting factor matrices (A, B, D, A*, B*, D*) and
# the weight matrices W*_n and W*_{n,m}.
#
# These are verification/research tools — not called by the main estimator
# pipeline, but essential for the theoretical analysis in the paper.
#
# Depends on: Jmat() from utils.R.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# 1. STRUCTURAL / SYMMETRY CHECKS
#------------------------------------------------------------------------------

check_symmetric <- function(A, tol = 1e-10) {
  max(abs(A - t(A))) < tol
}

check_centrosymmetric <- function(A, tol = 1e-10) {
  J <- Jmat(nrow(A))
  max(abs(A - J %*% A %*% J)) < tol
}

check_persymmetric <- function(A, tol = 1e-10) {
  J <- Jmat(nrow(A))
  max(abs(A - J %*% t(A) %*% J)) < tol
}

check_bisymmetric <- function(A, tol = 1e-10) {
  check_symmetric(A, tol) && check_centrosymmetric(A, tol)
}

check_toeplitz <- function(A, tol = 1e-8) {
  nr <- nrow(A)
  nc <- ncol(A)
  
  for (i in seq_len(nr - 1)) {
    for (j in seq_len(nc - 1)) {
      if (abs(A[i, j] - A[i + 1, j + 1]) > tol) {
        return(FALSE)
      }
    }
  }
  
  TRUE
}


check_bttb <- function(A, block_size, tol = 1e-8) {
  n <- nrow(A)
  
  if (!is.matrix(A)) return(FALSE)
  if (n != ncol(A)) return(FALSE)
  if (n %% block_size != 0) return(FALSE)
  
  nb <- n / block_size
  
  get_block <- function(i, j) {
    rows <- ((i - 1) * block_size + 1):(i * block_size)
    cols <- ((j - 1) * block_size + 1):(j * block_size)
    A[rows, cols, drop = FALSE]
  }
  
  # 1. Check each block is Toeplitz
  for (i in seq_len(nb)) {
    for (j in seq_len(nb)) {
      if (!check_toeplitz(get_block(i, j), tol = tol)) {
        return(FALSE)
      }
    }
  }
  
  # 2. Check block Toeplitz structure
  for (i in seq_len(nb - 1)) {
    for (j in seq_len(nb - 1)) {
      B1 <- get_block(i, j)
      B2 <- get_block(i + 1, j + 1)
      
      if (max(abs(B1 - B2)) > tol) {
        return(FALSE)
      }
    }
  }
  
  TRUE
}
check_block_toeplitz <- function(A, block_size, tol = 1e-8) {
  n <- nrow(A)
  
  # Check square matrix
  if (n != ncol(A)) return(FALSE)
  
  # Check divisibility
  if (n %% block_size != 0) return(FALSE)
  
  nb <- n / block_size  # number of blocks per dimension
  
  # Function to extract block (i,j)
  get_block <- function(i, j) {
    rows <- ((i - 1) * block_size + 1):(i * block_size)
    cols <- ((j - 1) * block_size + 1):(j * block_size)
    A[rows, cols, drop = FALSE]
  }
  
  # Check block Toeplitz property
  for (i in 1:(nb - 1)) {
    for (j in 1:(nb - 1)) {
      B1 <- get_block(i, j)
      B2 <- get_block(i + 1, j + 1)
      
      if (max(abs(B1 - B2)) > tol) {
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}
check_hankel <- function(A, tol = 1e-10) {
  check_toeplitz(A %*% Jmat(nrow(A)), tol)
}


#------------------------------------------------------------------------------
# 2. SIGN / ENTRY PATTERN CHECKS
#------------------------------------------------------------------------------

check_nonnegative <- function(A, tol = 1e-10) all(A >= -tol)

check_zero_row_sum <- function(A, tol = 1e-10) max(abs(rowSums(A))) < tol

check_zero_col_sum <- function(A, tol = 1e-10) max(abs(colSums(A))) < tol

check_row_stochastic <- function(A, tol = 1e-10) {
  check_nonnegative(A, tol) && all(abs(rowSums(A) - 1) < tol)
}

check_doubly_stochastic <- function(A, tol = 1e-10) {
  check_nonnegative(A, tol) &&
    all(abs(rowSums(A) - 1) < tol) &&
    all(abs(colSums(A) - 1) < tol)
}

check_Z_matrix <- function(A, tol = 1e-10) {
  all(A[row(A) != col(A)] <= tol)
}

check_M_matrix <- function(A, tol = 1e-10) {
  offdiag_nonpos <- all(A[row(A) != col(A)] <= tol)
  ev <- eigen(A, only.values = TRUE)$values
  offdiag_nonpos && all(Re(ev) >= -tol)
}


#------------------------------------------------------------------------------
# 3. SPECTRAL PROPERTIES
#------------------------------------------------------------------------------

#' Check whether a target value is a simple eigenvalue of A.
#' Returns found, algebraic/geometric multiplicity, and defectiveness.
check_simple_eigenvalue <- function(A, target, tol = 1e-10) {
  n  <- nrow(A)
  ev <- eigen(A, only.values = TRUE)$values
  matches  <- which(Mod(ev - target) < tol)
  alg_mult <- length(matches)

  sv <- svd(A - target * diag(n))$d
  geo_mult <- sum(sv < tol)

  list(found     = alg_mult > 0,
       alg_mult  = alg_mult,
       is_simple = alg_mult == 1,
       geo_mult  = geo_mult,
       defective = alg_mult > geo_mult)
}

#' Cluster eigenvalues and report their multiplicities.
get_eigenvalue_multiplicities <- function(eigenvalues, tol = 1e-10) {
  n <- length(eigenvalues)
  assigned <- rep(FALSE, n)
  clusters <- list()

  for (i in seq_len(n)) {
    if (assigned[i]) next
    dists   <- Mod(eigenvalues - eigenvalues[i])
    members <- which(dists < tol & !assigned)
    assigned[members] <- TRUE
    clusters[[length(clusters) + 1]] <- list(
      representative = eigenvalues[i],
      alg_mult       = length(members)
    )
  }

  data.frame(
    eigenvalue = sapply(clusters, function(cl) {
      v <- cl$representative
      if (abs(Im(v)) < tol) sprintf("%.8f", Re(v))
      else sprintf("%.6f%+.6fi", Re(v), Im(v))
    }),
    re       = sapply(clusters, function(cl) Re(cl$representative)),
    im       = sapply(clusters, function(cl) Im(cl$representative)),
    modulus  = sapply(clusters, function(cl) Mod(cl$representative)),
    alg_mult = sapply(clusters, function(cl) cl$alg_mult),
    is_simple = sapply(clusters, function(cl) cl$alg_mult == 1),
    stringsAsFactors = FALSE
  )
}

#' Comprehensive spectral summary of a square matrix.
get_spectral_summary <- function(A, tol = 1e-10) {
  n <- nrow(A)
  is_sym <- check_symmetric(A, tol)
  ev <- eigen(A, symmetric = is_sym)
  eigenvalues <- ev$values
  mods <- Mod(eigenvalues)
  sv   <- svd(A)$d

  mult_table <- get_eigenvalue_multiplicities(eigenvalues, tol)
  zero_info  <- check_simple_eigenvalue(A, 0, tol)

  list(
    eigenvalues            = eigenvalues,
    spectral_radius        = max(mods),
    spectral_radius_simple = sum(abs(mods - max(mods)) < tol) == 1,
    rank                   = sum(sv > tol),
    nullity                = sum(sv <= tol),
    condition_number       = if (min(sv) > tol) max(sv) / min(sv) else Inf,
    all_real_eig           = all(abs(Im(eigenvalues)) < tol),
    all_nonneg_eig         = is_sym && all(Re(eigenvalues) >= -tol),
    multiplicity_table     = mult_table,
    n_distinct_eigenvalues = nrow(mult_table),
    all_simple             = all(mult_table$is_simple),
    zero_eigenvalue_info   = zero_info
  )
}


#------------------------------------------------------------------------------
# 4. DEFINITENESS
#------------------------------------------------------------------------------

check_psd <- function(A, tol = 1e-10) {
  if (!check_symmetric(A, tol)) return(FALSE)
  all(eigen(A, symmetric = TRUE, only.values = TRUE)$values >= -tol)
}

check_pd <- function(A, tol = 1e-10) {
  if (!check_symmetric(A, tol)) return(FALSE)
  all(eigen(A, symmetric = TRUE, only.values = TRUE)$values > tol)
}

check_indefinite <- function(A, tol = 1e-10) {
  if (!check_symmetric(A, tol)) return(FALSE)
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  any(ev > tol) && any(ev < -tol)
}


#------------------------------------------------------------------------------
# 5. OTHER PROPERTIES
#------------------------------------------------------------------------------

check_normal <- function(A, tol = 1e-10) {
  max(abs(A %*% t(A) - t(A) %*% A)) < tol
}

check_idempotent <- function(A, tol = 1e-10) {
  max(abs(A %*% A - A)) < tol
}

check_involutory <- function(A, tol = 1e-10) {
  max(abs(A %*% A - diag(nrow(A)))) < tol
}

check_diag_dominant <- function(A, strict = FALSE) {
  d <- abs(diag(A))
  o <- rowSums(abs(A)) - d
  if (strict) all(d > o) else all(d >= o)
}


#------------------------------------------------------------------------------
# 6. PERRON-FROBENIUS PROPERTIES
#------------------------------------------------------------------------------

#' Check irreducibility of a nonnegative matrix.
#' A is irreducible iff (I + |A|)^{n-1} > 0.
check_irreducible <- function(A, tol = 1e-10) {
  n <- nrow(A)
  M <- diag(n) + abs(A)
  for (k in 2:(n - 1)) M <- M %*% (diag(n) + abs(A))
  all(M > tol)
}

#' Check primitivity: irreducible + aperiodic (exists k with A^k > 0).
check_primitive <- function(A, tol = 1e-10) {
  n <- nrow(A)
  if (!check_nonnegative(A, tol) || !check_irreducible(A, tol)) return(FALSE)
  M <- A
  bound <- min((n - 1)^2 + 1, n^2)
  for (k in seq_len(bound)) {
    if (all(M > tol)) return(TRUE)
    M <- M %*% A
  }
  FALSE
}

#' Commutativity check: [A, B] = 0.
check_commute <- function(A, B, tol = 1e-10) {
  max(abs(A %*% B - B %*% A)) < tol
}

#' Check the bisymmetry relation B = J A J.
check_B_equals_JAJ <- function(A, B, tol = 1e-10) {
  J <- Jmat(nrow(A))
  max(abs(B - J %*% A %*% J)) < tol
}


#------------------------------------------------------------------------------
# 7. FULL REPORT FOR A SINGLE MATRIX
#------------------------------------------------------------------------------

#' Run all relevant checks on a single matrix and return a named list of results.
full_matrix_report <- function(A, label = "A", block_size=2, tol = 1e-10) {
  n <- nrow(A)
  list(
    label              = label,
    dim                = n,
    symmetric          = check_symmetric(A, tol),
    centrosymmetric    = check_centrosymmetric(A, tol),
    bisymmetric        = check_bisymmetric(A, tol),
    toeplitz           = check_toeplitz(A, tol),
    block_toeplitz     = check_block_toeplitz(A, block_size = block_size, tol = tol),
    bttb                = check_bttb(A, block_size = block_size, tol = tol),
    nonnegative        = check_nonnegative(A, tol),
    zero_row_sum       = check_zero_row_sum(A, tol),
    row_stochastic     = check_row_stochastic(A, tol),
    Z_matrix           = check_Z_matrix(A, tol),
    diag_dominant      = check_diag_dominant(A),
    psd                = check_psd(A, tol),
    pd                 = check_pd(A, tol),
    normal             = check_normal(A, tol),
    spectral           = get_spectral_summary(A, tol)
  )
}


#------------------------------------------------------------------------------
# 8. WEIGHT MATRIX SPECIFIC DIAGNOSTICS
#------------------------------------------------------------------------------

#' Run diagnostics on the full W*_n and its 1D components.
#'
#' @param n1, n2  Grid dimensions.
#' @param Nh_normalize  Normalisation for counting matrices.
#' @return  A list of reports for A*_1, B*_1, D*_1, A*_2, B*_2, D*_2,
#'          and the full W*_n.
diagnose_weight_matrix <- function(n1, n2, Nh_normalize = TRUE) {

  mat1 <- build_count_matrices(n1, Nh_normalize)
  mat2 <- build_count_matrices(n2, Nh_normalize)
  W    <- build_W(mat1, mat2)

  list(
    A1_star = full_matrix_report(mat1$A_norm, "A*_1"),
    B1_star = full_matrix_report(mat1$B_norm, "B*_1"),
    D1_star = full_matrix_report(mat1$D_norm, "D*_1"),
    A2_star = full_matrix_report(mat2$A_norm, "A*_2"),
    B2_star = full_matrix_report(mat2$B_norm, "B*_2"),
    D2_star = full_matrix_report(mat2$D_norm, "D*_2"),
    W_star  = full_matrix_report(W$W_star, "W*_n"),
    # Specific checks for the paper
    B1_eq_JA1J = check_B_equals_JAJ(mat1$A_norm, mat1$B_norm),
    B2_eq_JA2J = check_B_equals_JAJ(mat2$A_norm, mat2$B_norm),
    A1_B1_commute = check_commute(mat1$A_norm, mat1$B_norm),
    A2_B2_commute = check_commute(mat2$A_norm, mat2$B_norm)
  )
}
