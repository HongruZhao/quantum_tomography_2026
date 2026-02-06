## ============================================================================
## 06_fisher_and_metrics.R
## Fisher information and loss metric matrices (Frobenius, Bures, Observable)
## Based on: 06_fisher_information_and_metrics.md
## ============================================================================

source("05_mle_cvxr_solver.R")

## ==========================================================================
## FISHER INFORMATION
## ==========================================================================

#' Compute per-shot Fisher information at measurement setting a
#' [I_a(θ)]_{jk} = (1/4) Σ_b s_j(a,b) s_k(a,b) / p_{a,b}(θ)
#' @param theta Bloch vector
#' @param a Setting index
#' @param S_ab Coefficient matrix (M x d)
#' @param c_ab Constant vector (length M)
#' @param ab_df Data frame with (a, b, row) mapping
#' @param eps Small constant for numerical stability
#' @return d x d Fisher information matrix for setting a
fisher_info_setting <- function(theta, a, S_ab, c_ab, ab_df, eps = 1e-12) {
  # Get rows for setting a
  idx <- which(ab_df$a == a)
  if (length(idx) == 0) stop(sprintf("No rows found for setting a = %d", a))

  S_a <- S_ab[idx, , drop = FALSE]  # r_a x d matrix
  c_a <- c_ab[idx]                   # length r_a

  # Compute probabilities at this setting
  p_a <- c_a + 0.5 * as.numeric(S_a %*% theta)
  p_a <- pmax(p_a, eps)  # Numerical safety

  # I_a = (1/4) * t(S_a) %*% diag(1/p_a) %*% S_a
  # Efficient computation: row-wise scaling
  Ia <- 0.25 * crossprod(S_a, S_a * (1 / p_a))

  Ia
}

#' Compute Fisher information for all settings
#' @param theta Bloch vector
#' @param S_ab Coefficient matrix
#' @param c_ab Constant vector
#' @param ab_df Data frame with (a, b, row) mapping
#' @param eps Small constant for numerical stability
#' @return List of Fisher information matrices, one per setting
fisher_info_all_settings <- function(theta, S_ab, c_ab, ab_df, eps = 1e-12) {
  A <- sort(unique(ab_df$a))
  I_list <- setNames(
    lapply(A, function(a) fisher_info_setting(theta, a, S_ab, c_ab, ab_df, eps)),
    A
  )
  I_list
}

#' Compute total plug-in Fisher information
#' J_n = Σ_a N_a I_a(θ̂_n)
#' @param theta Bloch vector (current estimate)
#' @param counts Vector of counts per setting (N_a)
#' @param S_ab Coefficient matrix
#' @param c_ab Constant vector
#' @param ab_df Data frame with (a, b, row) mapping
#' @param eps Small constant for numerical stability
#' @return List with J_total (d x d), I_list (per-setting FI)
total_fisher_info <- function(theta, counts, S_ab, c_ab, ab_df, eps = 1e-12) {
  A <- sort(unique(ab_df$a))
  d <- ncol(S_ab)

  I_list <- fisher_info_all_settings(theta, S_ab, c_ab, ab_df, eps)

  J_total <- matrix(0, d, d)
  for (i in seq_along(A)) {
    a <- A[i]
    n_a <- counts[i]
    if (n_a > 0) {
      J_total <- J_total + n_a * I_list[[as.character(a)]]
    }
  }

  list(J_total = J_total, I_list = I_list, counts = counts)
}

## ==========================================================================
## LOSS METRIC MATRICES
## ==========================================================================

## --------------------------------------------------------------------------
## Frobenius (Hilbert-Schmidt) Metric
## --------------------------------------------------------------------------

#' Compute Frobenius metric matrix in Bloch coordinates
#' G_F = (1/4) [tr(σ_j σ_k)]_{jk}
#' For orthonormal basis with tr(σ_j σ_k) = 2δ_jk, G_F = (1/2) I_d
#' @param sigmas List of basis matrices
#' @return d x d metric matrix (constant, does not depend on theta)
metric_frobenius <- function(sigmas) {
  d <- length(sigmas)
  G <- matrix(0, d, d)
  for (j in 1:d) {
    for (k in 1:d) {
      G[j, k] <- Re(traceC(sigmas[[j]] %*% sigmas[[k]])) / 4
    }
  }
  G
}

## --------------------------------------------------------------------------
## Bures (SLD) Metric
## --------------------------------------------------------------------------

#' Compute Lyapunov inverse: solve Ω_ρ(X) = Y for X
#' Ω_ρ(X) = (1/2)(ρX + Xρ)
#' In eigenbasis: [Ω^{-1}_ρ(Y)]_{ij} = 2/(λ_i + λ_j) Y_{ij}
#' @param rho Density matrix (must be full rank)
#' @param Y Hermitian matrix
#' @param tol Tolerance for eigenvalue floor
#' @return X such that Ω_ρ(X) ≈ Y
lyapunov_inv <- function(rho, Y, tol = 1e-10) {
  # Ensure rho is properly Hermitian
  rho_h <- hermitianize(rho)
  N <- nrow(rho_h)

  # For complex Hermitian matrices, use eigen WITHOUT symmetric=TRUE
  # (symmetric=TRUE only works for real symmetric matrices in R)
  eg <- eigen(rho_h, symmetric = FALSE)

  U <- eg$vectors
  # Eigenvalues of Hermitian matrix are real, but may have tiny imaginary parts
  lam <- pmax(Re(eg$values), tol)

  # Transform Y to eigenbasis
  Ytilde <- Conj(t(U)) %*% Y %*% U

  # Element-wise division: X_tilde[i,j] = 2/(λ_i + λ_j) * Y_tilde[i,j]
  denom <- outer(lam, lam, "+")
  denom <- pmax(denom, 2 * tol)  # Prevent division by very small numbers
  Xtilde <- (2 / denom) * Ytilde

  # Transform back
  X <- U %*% Xtilde %*% Conj(t(U))
  hermitianize(X)
}

#' Compute Bures metric matrix in Bloch coordinates
#' [G_B(θ)]_{jk} = (1/16) tr(σ_j Ω^{-1}_{ρ(θ)}(σ_k))
#' @param theta Bloch vector
#' @param sigmas List of basis matrices
#' @param N Hilbert space dimension
#' @param tol Tolerance for Lyapunov inverse
#' @return d x d metric matrix (state-dependent)
metric_bures <- function(theta, sigmas, N, tol = 1e-10) {
  rho <- rho_of_theta(theta, sigmas, N)
  d <- length(sigmas)

  G <- matrix(0, d, d)

  # Precompute Ω^{-1}(σ_k) for each k
  Omega_inv_sig <- vector("list", d)
  for (k in 1:d) {
    Omega_inv_sig[[k]] <- lyapunov_inv(rho, sigmas[[k]], tol = tol)
  }

  # Compute metric entries
  for (j in 1:d) {
    for (k in 1:d) {
      G[j, k] <- Re(traceC(sigmas[[j]] %*% Omega_inv_sig[[k]])) / 16
    }
  }

  # Ensure symmetry
  (G + t(G)) / 2
}

## --------------------------------------------------------------------------
## Observable-Expectation Metric
## --------------------------------------------------------------------------

#' Compute observable-expectation metric matrix in Bloch coordinates
#' G_O = (1/4) Σ_j o_j o_j^T, where [o_j]_k = tr(σ_k O_j)
#' @param sigmas List of basis matrices
#' @param Obs List of Hermitian observables
#' @return d x d metric matrix (constant, does not depend on theta)
metric_observable <- function(sigmas, Obs) {
  d <- length(sigmas)
  m <- length(Obs)

  G <- matrix(0, d, d)

  for (j in 1:m) {
    # Compute o_j vector
    oj <- numeric(d)
    for (k in 1:d) {
      oj[k] <- Re(traceC(sigmas[[k]] %*% Obs[[j]]))
    }
    # Add outer product
    G <- G + tcrossprod(oj)
  }

  G / 4
}

#' Default observable family for one qubit: {X, Y, Z}
#' @param pauli Output from build_pauli_basis_1q()
#' @return List of observables
default_observables_1q <- function(pauli = NULL) {
  if (is.null(pauli)) pauli <- build_pauli_basis_1q()
  list(X = pauli$sigmas$X, Y = pauli$sigmas$Y, Z = pauli$sigmas$Z)
}

#' Default observable family for two qubits: 15 Pauli products
#' @param basis Output from build_pauli_product_basis_2q()
#' @return List of observables (unscaled Pauli products)
default_observables_2q <- function(basis = NULL) {
  if (is.null(basis)) basis <- build_pauli_product_basis_2q()

  # Use unscaled Pauli products for observables
  pauli <- basis$pauli_1q
  single <- list(I = pauli$I, X = pauli$sigmas$X, Y = pauli$sigmas$Y, Z = pauli$sigmas$Z)

  order <- basis$labels
  Obs <- list()
  for (lab in order) {
    a <- substr(lab, 1, 1)
    b <- substr(lab, 2, 2)
    Obs[[lab]] <- kron(single[[a]], single[[b]])  # Unscaled
  }
  Obs
}

## ==========================================================================
## PROXY RISK COMPUTATION
## ==========================================================================

#' Compute A-optimal proxy risk R_n(G) = tr(G J_n^{-1})
#' @param G Metric matrix (d x d)
#' @param J Total Fisher information (d x d)
#' @param ridge Ridge parameter for numerical stability
#' @return Scalar proxy risk value
proxy_risk <- function(G, J, ridge = 0) {
  d <- nrow(J)
  J_use <- J + ridge * diag(d)
  J_inv <- tryCatch(
    solve(J_use),
    error = function(e) {
      # Use pseudo-inverse if singular
      svd_J <- svd(J_use)
      tol <- max(dim(J_use)) * max(svd_J$d) * .Machine$double.eps
      pos <- svd_J$d > tol
      if (sum(pos) == 0) return(matrix(Inf, d, d))
      svd_J$v[, pos, drop = FALSE] %*%
        diag(1/svd_J$d[pos], nrow = sum(pos)) %*%
        t(svd_J$u[, pos, drop = FALSE])
    }
  )
  sum(diag(G %*% J_inv))
}

## ==========================================================================
## METRIC FUNCTION FACTORIES
## ==========================================================================

#' Create metric function for Frobenius loss
#' Returns a function theta -> G (constant)
#' @param sigmas List of basis matrices
#' @return Function(theta) returning metric matrix
make_metric_frobenius <- function(sigmas) {
  G <- metric_frobenius(sigmas)
  function(theta) G
}

#' Create metric function for Bures loss
#' Returns a function theta -> G(theta) (state-dependent)
#' @param sigmas List of basis matrices
#' @param N Hilbert space dimension
#' @return Function(theta) returning metric matrix
make_metric_bures <- function(sigmas, N) {
  function(theta) metric_bures(theta, sigmas, N)
}

#' Create metric function for observable-expectation loss
#' Returns a function theta -> G (constant)
#' @param sigmas List of basis matrices
#' @param Obs List of observables
#' @return Function(theta) returning metric matrix
make_metric_observable <- function(sigmas, Obs) {
  G <- metric_observable(sigmas, Obs)
  function(theta) G
}

## ==========================================================================
## VALIDATION
## ==========================================================================

#' Validate Fisher information is PSD
#' @param Ia Fisher information matrix
#' @param tol Tolerance
#' @return List with valid (logical) and min_eigenvalue
validate_fisher_psd <- function(Ia, tol = 1e-10) {
  eigs <- eigen(Ia, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(Re(eigs))
  list(valid = min_eig >= -tol, min_eigenvalue = min_eig)
}

#' Validate Bures metric is PSD and symmetric
#' @param G Bures metric matrix
#' @param tol Tolerance
#' @return List with validation results
validate_bures_metric <- function(G, tol = 1e-10) {
  # Check symmetry
  sym_err <- max(abs(G - t(G)))

  # Check PSD
  eigs <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(Re(eigs))

  list(
    valid = (sym_err < tol) && (min_eig >= -tol),
    symmetric = sym_err < tol,
    psd = min_eig >= -tol,
    symmetry_error = sym_err,
    min_eigenvalue = min_eig
  )
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Example usage (uncomment to test)
# cat("\n=== Fisher Information Test ===\n")
#
# # Setup
# pauli <- build_pauli_basis_1q()
# sigmas <- pauli$sigmas
# L1 <- build_library_1q_L1()
# sc <- build_Sab_cab(sigmas, L1$Q_list, 2, L1$ab_df)
#
# # Random state
# theta_test <- c(0.3, -0.2, 0.5)
# cat("Test theta:", theta_test, "\n")
#
# # Fisher info for each setting
# for (a in 1:L1$k) {
#   Ia <- fisher_info_setting(theta_test, a, sc$S_ab, sc$c_ab, L1$ab_df)
#   cat(sprintf("I_%d eigenvalues: %s\n", a,
#               paste(round(eigen(Ia)$values, 4), collapse=", ")))
# }
#
# cat("\n=== Metric Matrices Test ===\n")
#
# # Frobenius
# G_F <- metric_frobenius(sigmas)
# cat("G_Frobenius:\n")
# print(round(G_F, 4))
#
# # Bures
# G_B <- metric_bures(theta_test, sigmas, 2)
# cat("\nG_Bures:\n")
# print(round(G_B, 4))
# val_B <- validate_bures_metric(G_B)
# cat("Bures metric valid:", val_B$valid, "\n")
#
# # Observable
# Obs <- default_observables_1q(pauli)
# G_O <- metric_observable(sigmas, Obs)
# cat("\nG_Observable:\n")
# print(round(G_O, 4))

cat("06_fisher_and_metrics.R loaded successfully\n")
