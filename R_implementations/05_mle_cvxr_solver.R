## ============================================================================
## 05_mle_cvxr_solver.R
## Stabilized MLE solver using CVXR with eigenvalue floor constraint
## Based on: 05_mle_cvxr_solver.md
## ============================================================================

source("04_measurement_library_2q.R")

## ==========================================================================
## AFFINE PROBABILITY MODEL: p_{a,b}(θ) = c_{a,b} + (1/2) s(a,b)^T θ
## ==========================================================================

#' Build the affine probability model matrices S_ab and c_ab
#' s_j(a,b) = tr(σ_j Q_{a,b}), c_{a,b} = tr(Q_{a,b})/N
#' @param sigmas List of basis matrices
#' @param Q_list List of POVMs
#' @param N Hilbert space dimension
#' @param ab_df Data frame with (a, b, row) mapping
#' @return List with S_ab (M x d matrix), c_ab (length M vector)
build_Sab_cab <- function(sigmas, Q_list, N, ab_df) {
  d <- length(sigmas)
  M <- nrow(ab_df)

  S_ab <- matrix(0, M, d)
  c_ab <- numeric(M)

  for (m in 1:M) {
    a <- ab_df$a[m]
    b <- ab_df$b[m]
    Q <- Q_list[[a]][[b]]

    # c_{a,b} = tr(Q_{a,b})/N
    c_ab[m] <- Re(traceC(Q)) / N

    # s_j(a,b) = tr(σ_j Q_{a,b})
    for (j in 1:d) {
      S_ab[m, j] <- Re(traceC(sigmas[[j]] %*% Q))
    }
  }

  list(S_ab = S_ab, c_ab = c_ab)
}

#' Compute probabilities from theta using affine model
#' p_{a,b}(θ) = c_{a,b} + (1/2) s(a,b)^T θ
#' @param theta Bloch vector
#' @param S_ab Coefficient matrix (M x d)
#' @param c_ab Constant vector (length M)
#' @return Probability vector (length M)
probs_from_theta <- function(theta, S_ab, c_ab) {
  as.numeric(c_ab + 0.5 * (S_ab %*% theta))
}

## ==========================================================================
## COUNT MANAGEMENT
## ==========================================================================

#' Build count vector Nab from sequence of (a,b) observations
#' @param a_seq Vector of setting indices
#' @param b_seq Vector of outcome indices
#' @param ab_row List mapping (a,b) to row index
#' @param M Total number of (a,b) cells
#' @return Count vector Nab
counts_from_ab <- function(a_seq, b_seq, ab_row, M) {
  Nab <- numeric(M)
  for (m in seq_along(a_seq)) {
    a <- a_seq[m]
    b <- b_seq[m]
    r <- ab_row[[a]][b]
    Nab[r] <- Nab[r] + 1
  }
  Nab
}

#' Get counts by setting (N_a = Σ_b N_{a,b})
#' @param Nab Count vector
#' @param ab_df Data frame with (a, b, row) mapping
#' @return Vector of counts per setting
counts_by_setting <- function(Nab, ab_df) {
  A <- sort(unique(ab_df$a))
  vapply(A, function(a) sum(Nab[ab_df$a == a]), numeric(1))
}

#' Increment count for a single (a,b) observation
#' @param Nab Current count vector
#' @param a Setting index
#' @param b Outcome index
#' @param ab_row List mapping (a,b) to row index
#' @return Updated count vector
increment_count <- function(Nab, a, b, ab_row) {
  idx <- ab_row[[a]][b]
  Nab[idx] <- Nab[idx] + 1
  Nab
}

## ==========================================================================
## STABILIZED MLE WITH EIGENVALUE FLOOR (CVXR)
## ==========================================================================

#' Fit theta via MLE with stabilized eigenvalue floor constraint
#' Solves: min_θ -Σ_{a,b} N_{a,b} log(p_{a,b}(θ))
#'         s.t. ρ(θ) ⪰ η I
#' @param N Hilbert space dimension
#' @param sigmas List of basis matrices
#' @param S_ab Coefficient matrix (M x d)
#' @param c_ab Constant vector (length M)
#' @param Nab Count vector (length M)
#' @param eta Eigenvalue floor (0 < eta < 1/N, default 0)
#' @param solver CVXR solver ("MOSEK", "ECOS", "SCS")
#' @param eps_log Small constant for log stability
#' @param verbose Print solver output
#' @return List with theta_hat, rho_hat, status, value
fit_theta_cvxr <- function(N, sigmas, S_ab, c_ab, Nab,
                           eta = 0,
                           solver = c("SCS", "ECOS", "MOSEK"),
                           eps_log = 1e-12,
                           verbose = FALSE) {

  solver <- match.arg(solver)
  d <- length(sigmas)
  M <- length(c_ab)

  # Check eta validity
  if (eta >= 1/N) {
    warning(sprintf("eta = %.4f must be < 1/N = %.4f; setting eta = 0", eta, 1/N))
    eta <- 0
  }

  # CVXR variable: theta

theta_var <- Variable(d)

  # Affine probabilities: p = c_ab + 0.5 * S_ab %*% theta
  S_const <- as.matrix(S_ab)
  p_expr <- c_ab + 0.5 * (S_const %*% theta_var)

  # Real-embedded PSD constraint on rho(theta)
  # ρ(θ) = I/N + (1/2) Σ_j θ_j σ_j
  SI <- diag(2 * N)  # Real-embedded identity
  Slist <- lapply(sigmas, real_embed)

  # Build affine expression for real-embedded rho
  A_affine <- SI / N
  for (j in 1:d) {
    A_affine <- A_affine + 0.5 * theta_var[j] * Slist[[j]]
  }

  # PSD variable representing real-embedded rho
  rho_var <- Variable(2 * N, 2 * N, PSD = TRUE)

  # Negative log-likelihood objective
  # Only include terms where Nab > 0 to avoid log(0) issues
  active_idx <- which(Nab > 0)
  if (length(active_idx) == 0) {
    # No data: return zero theta
    theta_hat <- rep(0, d)
    rho_hat <- rho_of_theta(theta_hat, sigmas, N)
    return(list(
      theta_hat = theta_hat,
      rho_hat = rho_hat,
      status = "trivial",
      value = 0
    ))
  }

  obj <- -sum_entries(Nab[active_idx] * log(p_expr[active_idx] + eps_log))

  # Constraints
  constraints <- list(
    rho_var == A_affine,  # rho PSD via variable
    p_expr >= eps_log     # Probabilities positive
  )

  # Add eigenvalue floor constraint if eta > 0
  # ρ(θ) ⪰ η I  ⟺  rho_embed - η * SI ⪰ 0
  if (eta > 0) {
    rho_floor <- Variable(2 * N, 2 * N, PSD = TRUE)
    constraints <- c(constraints, list(
      rho_floor == A_affine - eta * SI
    ))
  }

  # Define and solve problem
  prob <- Problem(Minimize(obj), constraints)

  # Solver selection with fallback
  if (solver == "MOSEK" && !requireNamespace("Rmosek", quietly = TRUE)) {
    solver <- "SCS"
  }
  if (solver == "ECOS" && !requireNamespace("ECOSolveR", quietly = TRUE)) {
    solver <- "SCS"
  }

  # Solver parameters
  res <- tryCatch({
    solve(prob, solver = solver, verbose = verbose,
          feastol = 1e-8, reltol = 1e-8, abstol = 1e-8,
          max_iters = 5000)
  }, error = function(e) {
    list(status = "error", value = NA, error_message = e$message)
  })

  # Extract solution
  if (res$status %in% c("optimal", "optimal_inaccurate")) {
    theta_hat <- drop(res$getValue(theta_var))
    rho_hat <- rho_of_theta(theta_hat, sigmas, N)
  } else {
    theta_hat <- rep(NA, d)
    rho_hat <- matrix(NA, N, N)
  }

  list(
    theta_hat = theta_hat,
    rho_hat = rho_hat,
    status = res$status,
    value = res$value
  )
}

## ==========================================================================
## CONVENIENCE WRAPPERS
## ==========================================================================

#' One-shot MLE fit given measurement library and data
#' @param lib Measurement library object
#' @param sigmas List of basis matrices
#' @param a_seq Vector of setting indices
#' @param b_seq Vector of outcome indices
#' @param eta Eigenvalue floor
#' @param solver CVXR solver
#' @param verbose Print solver output
#' @return MLE fit result
fit_mle_from_data <- function(lib, sigmas, a_seq, b_seq,
                              eta = 0, solver = "SCS", verbose = FALSE) {

  N <- lib$N

  # Build affine model if not cached
  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)

  # Build counts
  Nab <- counts_from_ab(a_seq, b_seq, lib$ab_row, lib$M)

  # Fit MLE
  fit_theta_cvxr(N, sigmas, sc$S_ab, sc$c_ab, Nab,
                 eta = eta, solver = solver, verbose = verbose)
}

#' Compute negative log-likelihood at a given theta
#' @param theta Bloch vector
#' @param S_ab Coefficient matrix
#' @param c_ab Constant vector
#' @param Nab Count vector
#' @param eps Small constant for log stability
#' @return NLL value
nll_theta <- function(theta, S_ab, c_ab, Nab, eps = 1e-12) {
  p <- probs_from_theta(theta, S_ab, c_ab)
  if (any(p < -1e-8)) return(Inf)
  p <- pmax(p, eps)
  -sum(Nab * log(p))
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Example usage (uncomment to test)
# cat("\n=== MLE Test (One Qubit, Library 1) ===\n")
#
# # Setup
# pauli <- build_pauli_basis_1q()
# sigmas <- pauli$sigmas
# L1 <- build_library_1q_L1()
#
# # Generate true state
# rho_true <- random_density_fullrank(2, eps_mix = 0.1, seed = 123)
# theta_true <- theta_from_rho(rho_true, sigmas)
# cat("True theta:", round(theta_true, 4), "\n")
#
# # Simulate data
# set.seed(456)
# n_shots <- 100
# a_seq <- sample(1:L1$k, n_shots, replace = TRUE)  # Uniform sampling
# b_seq <- sample_outcome_sequence(a_seq, rho_true, L1$Q_list)
#
# # Fit MLE
# fit <- fit_mle_from_data(L1, sigmas, a_seq, b_seq, eta = 1e-4, solver = "SCS")
# cat("Solver status:", fit$status, "\n")
# cat("Estimated theta:", round(fit$theta_hat, 4), "\n")
# cat("Theta error:", round(sqrt(sum((fit$theta_hat - theta_true)^2)), 4), "\n")

cat("05_mle_cvxr_solver.R loaded successfully\n")
