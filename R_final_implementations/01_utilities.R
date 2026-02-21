## ============================================================================
## 01_utilities.R
## Core utility functions for quantum tomography simulation
## Based on: 01_existing_r_functions.md
## ============================================================================

CVXR_AVAILABLE <- requireNamespace("CVXR", quietly = TRUE)
suppressPackageStartupMessages({
  if (CVXR_AVAILABLE) library(CVXR)
  library(Matrix)
})
if (!CVXR_AVAILABLE) {
  message("CVXR package not found. CVXR-based solvers are unavailable; use solver = 'PGD'.")
}

## --------------------------------------------------------------------------
## Complex-safe trace
## --------------------------------------------------------------------------
#' Compute trace of a complex matrix
#' @param A Complex matrix
#' @return Real part of trace
traceC <- function(A) {
  Re(sum(diag(A)))
}

## --------------------------------------------------------------------------
## Hermitianize: enforce Hermiticity
## --------------------------------------------------------------------------
#' Enforce Hermiticity on a matrix
#' @param A Matrix (possibly complex)
#' @return Hermitian matrix (A + A^dagger)/2
hermitianize <- function(A) {
  (A + Conj(t(A))) / 2
}

## --------------------------------------------------------------------------
## Real embedding for complex matrices
## --------------------------------------------------------------------------
#' Real embedding: C^{N x N} --> R^{2N x 2N}
#' M ↦ [[Re M, -Im M], [Im M, Re M]]
#' Used to express complex PSD constraints in CVXR
#' @param M Complex N x N matrix
#' @return Real 2N x 2N matrix
real_embed <- function(M) {
  R <- Re(M)
  I <- Im(M)
  rbind(cbind(R, -I),
        cbind(I,  R))
}

## --------------------------------------------------------------------------
## Elementary matrix |j><k|
## --------------------------------------------------------------------------
#' Create elementary matrix with 1 at position (j,k)
#' @param j Row index
#' @param k Column index
#' @param N Dimension
#' @return N x N complex matrix with 1 at (j,k)
E_jk <- function(j, k, N) {
  M <- matrix(0 + 0i, N, N)
  M[j, k] <- 1 + 0i
  M
}

## --------------------------------------------------------------------------
## Rank-1 projector from state vector
## --------------------------------------------------------------------------
#' Create rank-1 projector from a state vector
#' @param psi Column vector (state)
#' @return Density matrix |psi><psi|
proj_rank1 <- function(psi) {
  psi <- matrix(psi, ncol = 1)
  psi %*% Conj(t(psi))
}

## --------------------------------------------------------------------------
## Random full-rank density matrix generator
## --------------------------------------------------------------------------
#' Generate a random full-rank density matrix
#' Uses Ginibre ensemble with identity mixing for guaranteed full rank
#' @param N Hilbert space dimension
#' @param eps_mix Mixing parameter with identity (default 0.05)
#' @param seed Random seed (optional)
#' @return N x N density matrix with all eigenvalues > 0
random_density_fullrank <- function(N, eps_mix = 0.05, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  # Draw complex Ginibre matrix
  G <- matrix(rnorm(N * N) + 1i * rnorm(N * N), N, N)
  # Form Wishart-type PSD matrix
  W <- G %*% Conj(t(G))
  # Normalize to trace 1
  rho0 <- W / Re(sum(diag(W)))
  # Mix with identity to ensure full rank
  rho <- (1 - eps_mix) * rho0 + eps_mix * diag(N) / N
  hermitianize(rho)
}

## --------------------------------------------------------------------------
## Random pure state + mixture with maximally mixed
## --------------------------------------------------------------------------
#' Generate a random pure-state density matrix
#' @param N Hilbert space dimension
#' @param seed Random seed (optional)
#' @return N x N rank-1 density matrix
random_pure_density <- function(N, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  psi <- rnorm(N) + 1i * rnorm(N)
  psi <- psi / sqrt(sum(Mod(psi)^2))
  hermitianize(proj_rank1(psi))
}

#' Mix a pure state with the maximally mixed state
#' rho = alpha * rho_pure + (1 - alpha) * I/N
#' @param rho_pure Rank-1 density matrix
#' @param alpha Mixing weight in [0,1]
#' @return Mixed density matrix
mix_with_maximally_mixed <- function(rho_pure, alpha) {
  if (alpha < 0 || alpha > 1) stop("alpha must be in [0,1]")
  N <- nrow(rho_pure)
  rho_mixed <- diag(N) / N
  hermitianize(alpha * rho_pure + (1 - alpha) * rho_mixed)
}

#' Mix two density matrices with weight alpha
#' rho = alpha * rho_a + (1 - alpha) * rho_b
#' @param rho_a Density matrix
#' @param rho_b Density matrix
#' @param alpha Mixing weight in [0,1]
#' @return Mixed density matrix
mix_with_state <- function(rho_a, rho_b, alpha) {
  if (alpha < 0 || alpha > 1) stop("alpha must be in [0,1]")
  hermitianize(alpha * rho_a + (1 - alpha) * rho_b)
}

## --------------------------------------------------------------------------
## Fixed pure states (non-random)
## --------------------------------------------------------------------------
#' Fixed 1-qubit pure state |+><+| (Bloch vector (1,0,0))
#' @return 2x2 density matrix
fixed_pure_state_1q_plus <- function() {
  rho_from_bloch_1q(c(1, 0, 0))
}

#' Fixed 1-qubit pure state with phase: cos(phi)|0> + i sin(phi)|1>
#' Default phi=0.2 rad
#' @param phi Angle in radians
#' @return 2x2 density matrix
fixed_pure_state_1q_yphase <- function(phi = 0.2) {
  psi <- c(cos(phi), 1i * sin(phi))
  hermitianize(proj_rank1(psi))
}

#' Fixed 1-qubit pure state |0><0|
#' @return 2x2 density matrix
fixed_pure_state_1q_0 <- function() {
  diag(c(1, 0))
}

#' Fixed 2-qubit pure state |00><00|
#' @return 4x4 density matrix
fixed_pure_state_2q_00 <- function() {
  psi <- c(1, 0, 0, 0)
  hermitianize(proj_rank1(psi))
}

#' Fixed 2-qubit Bell state |Phi+><Phi+| = (|00>+|11>)/sqrt(2)
#' @return 4x4 density matrix
fixed_pure_state_2q_phi_plus <- function() {
  psi <- c(1, 0, 0, 1) / sqrt(2)
  hermitianize(proj_rank1(psi))
}

#' Fixed 2-qubit phase-entangled state cos(phi)|00> + i sin(phi)|11>
#' @param phi Angle in radians
#' @return 4x4 density matrix
fixed_pure_state_2q_phi_phase <- function(phi = 0.2) {
  psi <- c(cos(phi), 0, 0, 1i * sin(phi))
  hermitianize(proj_rank1(psi))
}

## --------------------------------------------------------------------------
## Fixed non-symmetric states (non-random)
## --------------------------------------------------------------------------
#' Fixed 1-qubit non-symmetric mixed state (Bloch vector not axis-aligned)
#' @return 2x2 density matrix
fixed_nonsym_state_1q <- function() {
  rho_from_bloch_1q(c(0.6, 0.2, 0.1))
}

#' Fixed 2-qubit non-symmetric mixed state (not diagonal)
#' Constructed as a convex mix of a diagonal state and a pure off-diagonal state.
#' @return 4x4 density matrix
fixed_nonsym_state_2q <- function() {
  D <- diag(c(0.55, 0.25, 0.15, 0.05))
  psi <- c(1, 0.3, 0.2i, 0.1)
  psi <- psi / sqrt(sum(Mod(psi)^2))
  rho_pure <- proj_rank1(psi)
  hermitianize(0.7 * D + 0.3 * rho_pure)
}

## --------------------------------------------------------------------------
## Custom mixture rule
## --------------------------------------------------------------------------
#' Mix a pure state with a non-symmetric state using weight alpha/2
#' rho = (alpha/2) * rho_pure + (1 - alpha/2) * rho_nonsym
#' @param rho_pure Rank-1 density matrix
#' @param rho_nonsym Density matrix (non-symmetric)
#' @param alpha Mixing parameter in [0,1]
#' @return Mixed density matrix
mix_with_nonsym <- function(rho_pure, rho_nonsym, alpha) {
  if (alpha < 0 || alpha > 1) stop("alpha must be in [0,1]")
  hermitianize((alpha / 2) * rho_pure + (1 - alpha / 2) * rho_nonsym)
}

## --------------------------------------------------------------------------
## Stable seed from a string key (order-independent, reproducible)
## --------------------------------------------------------------------------
#' Deterministic seed derived from a string key
#' @param key String identifier (e.g., case_id, policy, replicate)
#' @param base Optional integer offset
#' @return Integer seed in [1, .Machine$integer.max]
stable_seed <- function(key, base = 0L) {
  key <- paste(as.character(key), collapse = "|")
  ints <- utf8ToInt(key)
  h <- 0.0
  for (v in ints) {
    h <- (h * 31 + as.numeric(v)) %% 2147483647
  }
  seed <- (h + as.numeric(base)) %% 2147483647
  seed <- as.integer(seed)
  if (is.na(seed) || seed == 0L) seed <- 1L
  seed
}

## --------------------------------------------------------------------------
## Sample random unit vector on d-sphere
## --------------------------------------------------------------------------
#' Sample a random unit vector uniformly on the d-sphere
#' @param d Dimension of vector
#' @return Unit vector of length d
sample_random_unit_vector <- function(d) {
  z <- rnorm(d)
  z / sqrt(sum(z^2))
}

## --------------------------------------------------------------------------
## Check if matrix is positive semidefinite
## --------------------------------------------------------------------------
#' Check if a Hermitian matrix is PSD
#' @param M Hermitian matrix
#' @param tol Tolerance for negative eigenvalues
#' @return TRUE if all eigenvalues >= -tol
is_psd <- function(M, tol = 1e-10) {
  eigs <- eigen(hermitianize(M), symmetric = TRUE, only.values = TRUE)$values
  all(Re(eigs) >= -tol)
}

## --------------------------------------------------------------------------
## Safe matrix inverse with ridge regularization
## --------------------------------------------------------------------------
#' Compute matrix inverse with optional ridge regularization
#' @param M Matrix to invert
#' @param ridge Ridge parameter (default 0)
#' @return Inverse of (M + ridge * I)
safe_solve <- function(M, ridge = 0) {
  d <- nrow(M)
  M_reg <- M + ridge * diag(d)
  tryCatch(
    solve(M_reg),
    error = function(e) {
      # If still singular, use pseudo-inverse
      svd_M <- svd(M_reg)
      tol <- max(dim(M_reg)) * max(svd_M$d) * .Machine$double.eps
      pos <- svd_M$d > tol
      svd_M$v[, pos, drop = FALSE] %*%
        diag(1/svd_M$d[pos], nrow = sum(pos)) %*%
        t(svd_M$u[, pos, drop = FALSE])
    }
  )
}

## --------------------------------------------------------------------------
## Kronecker product helper
## --------------------------------------------------------------------------
#' Kronecker (tensor) product of two matrices
#' @param A First matrix
#' @param B Second matrix
#' @return Kronecker product A ⊗ B
kron <- function(A, B) {
  kronecker(A, B)
}

## --------------------------------------------------------------------------
## Validate density matrix
## --------------------------------------------------------------------------
#' Validate that a matrix is a valid density matrix
#' @param rho Matrix to validate
#' @param tol Numerical tolerance
#' @return List with valid (logical) and messages (character vector)
validate_density <- function(rho, tol = 1e-10) {
  messages <- character()
  valid <- TRUE

  # Check Hermiticity
  if (max(Mod(rho - Conj(t(rho)))) > tol) {
    valid <- FALSE
    messages <- c(messages, "Not Hermitian")
  }

  # Check trace = 1
  tr_val <- Re(sum(diag(rho)))
  if (abs(tr_val - 1) > tol) {
    valid <- FALSE
    messages <- c(messages, sprintf("Trace = %.6f (not 1)", tr_val))
  }

  # Check PSD
  eigs <- eigen(hermitianize(rho), symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(Re(eigs))
  if (min_eig < -tol) {
    valid <- FALSE
    messages <- c(messages, sprintf("Min eigenvalue = %.2e (negative)", min_eig))
  }

  list(valid = valid, messages = messages, min_eigenvalue = min_eig)
}

## --------------------------------------------------------------------------
## Print utility for complex matrices
## --------------------------------------------------------------------------
#' Pretty print a complex matrix
#' @param M Complex matrix
#' @param digits Number of digits to display
print_complex_matrix <- function(M, digits = 4) {
  R <- round(Re(M), digits)
  I <- round(Im(M), digits)

  for (i in 1:nrow(M)) {
    row_str <- sapply(1:ncol(M), function(j) {
      if (abs(I[i,j]) < 10^(-digits)) {
        sprintf("%.*f", digits, R[i,j])
      } else if (I[i,j] >= 0) {
        sprintf("%.*f+%.*fi", digits, R[i,j], digits, I[i,j])
      } else {
        sprintf("%.*f%.*fi", digits, R[i,j], digits, I[i,j])
      }
    })
    cat(paste(row_str, collapse = "  "), "\n")
  }
}

cat("01_utilities.R loaded successfully\n")
