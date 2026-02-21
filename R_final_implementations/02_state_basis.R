## ============================================================================
## 02_state_basis.R
## Bloch basis construction for one-qubit and two-qubit systems
## Based on: 02_state_basis_one_two_qubit.md
## ============================================================================

source("01_utilities.R")

## ==========================================================================
## ONE QUBIT (N = 2, d = 3)
## ==========================================================================

## --------------------------------------------------------------------------
## Pauli matrices
## --------------------------------------------------------------------------
#' Build Pauli basis for one qubit
#' Returns standard Pauli matrices X, Y, Z with tr(σ_j σ_k) = 2δ_jk
#' @return List with N=2, sigmas (list of X,Y,Z), I (identity)
build_pauli_basis_1q <- function() {
  # Identity
  I2 <- diag(2) + 0i

  # Pauli X: [[0,1],[1,0]]
  X <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE) + 0i

  # Pauli Y: [[0,-i],[i,0]]
  Y <- matrix(c(0, -1i, 1i, 0), 2, 2, byrow = TRUE)

  # Pauli Z: [[1,0],[0,-1]]
  Z <- matrix(c(1, 0, 0, -1), 2, 2, byrow = TRUE) + 0i

  list(
    N = 2,
    d = 3,
    sigmas = list(X = X, Y = Y, Z = Z),
    I = I2,
    labels = c("X", "Y", "Z")
  )
}

## --------------------------------------------------------------------------
## One-qubit Bloch state
## --------------------------------------------------------------------------
#' Construct density matrix from Bloch vector (one qubit)
#' ρ(θ) = (1/2)(I + θ_1 X + θ_2 Y + θ_3 Z)
#' @param theta Bloch vector (length 3)
#' @param pauli Output from build_pauli_basis_1q() (optional)
#' @return 2x2 density matrix
rho_from_bloch_1q <- function(theta, pauli = NULL) {
  if (is.null(pauli)) pauli <- build_pauli_basis_1q()
  theta <- as.numeric(theta)
  stopifnot(length(theta) == 3)

  rho <- pauli$I / 2
  for (j in 1:3) {
    rho <- rho + 0.5 * theta[j] * pauli$sigmas[[j]]
  }
  hermitianize(rho)
}

#' Extract Bloch vector from one-qubit density matrix
#' θ_j = tr(ρ σ_j)
#' @param rho 2x2 density matrix
#' @param pauli Output from build_pauli_basis_1q() (optional)
#' @return Bloch vector (length 3)
bloch_from_rho_1q <- function(rho, pauli = NULL) {
  if (is.null(pauli)) pauli <- build_pauli_basis_1q()
  theta <- numeric(3)
  for (j in 1:3) {
    theta[j] <- Re(traceC(rho %*% pauli$sigmas[[j]]))
  }
  theta
}

## ==========================================================================
## TWO QUBITS (N = 4, d = 15)
## ==========================================================================

## --------------------------------------------------------------------------
## Pauli-product basis for two qubits
## --------------------------------------------------------------------------
#' Build Pauli-product basis for two qubits
#' Uses scaled basis: σ_{(α,β)} = (1/√2)(σ_α ⊗ σ_β) for (α,β) ≠ (I,I)
#' This scaling makes tr(σ_j σ_k) = 2δ_jk, matching the one-qubit basis
#' @param order Character vector specifying ordering of 15 basis elements
#' @return List with N=4, d=15, sigmas, labels
build_pauli_product_basis_2q <- function(
    order = c("IX","IY","IZ","XI","XX","XY","XZ","YI","YX","YY","YZ","ZI","ZX","ZY","ZZ")
) {
  # Get one-qubit Paulis
  pauli <- build_pauli_basis_1q()
  I <- pauli$I
  X <- pauli$sigmas$X
  Y <- pauli$sigmas$Y
  Z <- pauli$sigmas$Z

  # Map labels to matrices
  single <- list(I = I, X = X, Y = Y, Z = Z)

  # Build the 15 traceless products with scaling 1/sqrt(2)
  sigmas <- list()
  for (lab in order) {
    a <- substr(lab, 1, 1)
    b <- substr(lab, 2, 2)
    # Scale by 1/sqrt(2) so tr(σ_j σ_k) = 2δ_jk
    sigmas[[lab]] <- (1 / sqrt(2)) * kron(single[[a]], single[[b]])
  }

  list(
    N = 4,
    d = 15,
    sigmas = sigmas,
    I = diag(4) + 0i,
    labels = order,
    pauli_1q = pauli
  )
}

## --------------------------------------------------------------------------
## Two-qubit Bloch state
## --------------------------------------------------------------------------
#' Construct density matrix from Bloch vector (two qubits)
#' With our scaled basis (tr(σ_j σ_k)=2δ): ρ(θ) = I/N + (1/2)Σ_j θ_j σ_j
#' If P = σ_α⊗σ_β (unscaled Pauli products), then σ_j = P/√2 and
#' θ_j = tr(ρ σ_j) = (1/√2) tr(ρ P).
#' @param theta Bloch vector (length 15)
#' @param basis Output from build_pauli_product_basis_2q() (optional)
#' @return 4x4 density matrix
rho_from_bloch_2q <- function(theta, basis = NULL) {
  if (is.null(basis)) basis <- build_pauli_product_basis_2q()
  theta <- as.numeric(theta)
  stopifnot(length(theta) == 15)

  N <- 4
  rho <- basis$I / N
  for (j in 1:15) {
    rho <- rho + 0.5 * theta[j] * basis$sigmas[[j]]
  }
  hermitianize(rho)
}

#' Extract Bloch vector from two-qubit density matrix
#' @param rho 4x4 density matrix
#' @param basis Output from build_pauli_product_basis_2q() (optional)
#' @return Bloch vector (length 15)
bloch_from_rho_2q <- function(rho, basis = NULL) {
  if (is.null(basis)) basis <- build_pauli_product_basis_2q()
  theta_from_rho(rho, basis$sigmas)
}

## ==========================================================================
## GENERIC BLOCH MODEL FUNCTIONS
## ==========================================================================

#' Generic function to construct rho from theta given basis
#' ρ(θ) = I/N + (1/2)Σ_j θ_j σ_j
#' @param theta Bloch vector
#' @param sigmas List of basis matrices
#' @param N Hilbert space dimension
#' @return Density matrix
rho_of_theta <- function(theta, sigmas, N) {
  d <- length(sigmas)
  theta <- as.numeric(theta)
  if (length(theta) != d) {
    stop(sprintf("theta has length %d but sigmas has %d elements", length(theta), d))
  }

  rho <- diag(N) / N + 0i
  for (j in 1:d) {
    rho <- rho + 0.5 * theta[j] * sigmas[[j]]
  }
  hermitianize(rho)
}

#' Generic function to extract theta from rho given basis
#' θ_j = tr(ρ σ_j)
#' @param rho Density matrix
#' @param sigmas List of basis matrices
#' @return Bloch vector
theta_from_rho <- function(rho, sigmas) {
  d <- length(sigmas)
  theta <- numeric(d)
  for (j in 1:d) {
    denom <- Re(traceC(sigmas[[j]] %*% sigmas[[j]]))
    if (denom == 0) stop("Encountered zero-norm basis element in theta_from_rho().")
    # For basis with tr(σ_j σ_j) ≠ 2, rescale so that rho = I/N + 0.5 Σ θ_j σ_j
    theta[j] <- 2 * Re(traceC(rho %*% sigmas[[j]])) / denom
  }
  theta
}

## ==========================================================================
## VALIDATION FUNCTIONS
## ==========================================================================

#' Validate Pauli basis orthonormality
#' Check tr(σ_j σ_k) = 2δ_jk
#' @param sigmas List of basis matrices
#' @param tol Numerical tolerance
#' @return List with valid (logical) and Gram matrix
validate_basis_orthonormality <- function(sigmas, tol = 1e-10) {
  d <- length(sigmas)
  G <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      G[i, j] <- Re(traceC(sigmas[[i]] %*% sigmas[[j]]))
    }
  }

  # Check G = 2*I
  expected <- 2 * diag(d)
  max_diff <- max(abs(G - expected))

  list(
    valid = max_diff < tol,
    Gram_matrix = G,
    max_deviation = max_diff
  )
}

#' Validate basis tracelessness
#' Check tr(σ_j) = 0 for all j
#' @param sigmas List of basis matrices
#' @param tol Numerical tolerance
#' @return List with valid (logical) and traces
validate_basis_traceless <- function(sigmas, tol = 1e-10) {
  traces <- sapply(sigmas, function(s) Re(traceC(s)))
  max_trace <- max(abs(traces))

  list(
    valid = max_trace < tol,
    traces = traces,
    max_trace = max_trace
  )
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Example usage (uncomment to test)
# cat("\n=== One Qubit Basis ===\n")
# pauli_1q <- build_pauli_basis_1q()
# cat("Dimension N =", pauli_1q$N, ", d =", pauli_1q$d, "\n")
# cat("Basis labels:", paste(pauli_1q$labels, collapse=", "), "\n")
#
# # Test round-trip: rho -> theta -> rho
# rho_test <- random_density_fullrank(2, seed=123)
# theta_test <- bloch_from_rho_1q(rho_test, pauli_1q)
# rho_recovered <- rho_from_bloch_1q(theta_test, pauli_1q)
# cat("Round-trip error:", max(Mod(rho_test - rho_recovered)), "\n")
#
# cat("\n=== Two Qubit Basis ===\n")
# pauli_2q <- build_pauli_product_basis_2q()
# cat("Dimension N =", pauli_2q$N, ", d =", pauli_2q$d, "\n")
# cat("Basis labels:", paste(pauli_2q$labels, collapse=", "), "\n")

cat("02_state_basis.R loaded successfully\n")
