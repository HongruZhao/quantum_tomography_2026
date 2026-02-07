## ============================================================================
## 09_validation.R
## Validation and sanity checks for quantum tomography simulation
## Based on: 09_validation_and_sanity_checks.md
## ============================================================================

source("08_simulation_controller.R")

## ==========================================================================
## STATE PARAMETRIZATION CHECKS
## ==========================================================================

#' Validate rho_of_theta produces valid density matrices
#' @param sigmas Basis matrices
#' @param N Hilbert space dimension
#' @param n_tests Number of random tests
#' @param tol Numerical tolerance
#' @return List with validation results
validate_rho_of_theta <- function(sigmas, N, n_tests = 100, tol = 1e-10) {
  d <- length(sigmas)
  failures <- list()

  for (i in 1:n_tests) {
    # Generate random theta (within Bloch ball)
    theta <- rnorm(d)
    theta <- theta / (sqrt(sum(theta^2)) + 1) * 0.5  # Scale to be valid

    rho <- rho_of_theta(theta, sigmas, N)

    # Check Hermiticity
    herm_err <- max(Mod(rho - Conj(t(rho))))
    if (herm_err > tol) {
      failures <- c(failures, list(list(
        test = i, issue = "Hermiticity", error = herm_err
      )))
    }

    # Check trace = 1
    tr_err <- abs(Re(traceC(rho)) - 1)
    if (tr_err > tol) {
      failures <- c(failures, list(list(
        test = i, issue = "Trace", error = tr_err
      )))
    }
  }

  list(
    n_tests = n_tests,
    n_failures = length(failures),
    valid = length(failures) == 0,
    failures = failures
  )
}

#' Validate eigenvalue floor is respected by MLE
#' @param lib Measurement library
#' @param sigmas Basis matrices
#' @param eta Target eigenvalue floor
#' @param n_tests Number of tests
#' @param tol Numerical tolerance
#' @return Validation results
validate_mle_eigenvalue_floor <- function(lib, sigmas, eta = 1e-3,
                                          n_tests = 10, tol = 1e-8) {
  N <- lib$N
  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)

  failures <- list()

  for (i in 1:n_tests) {
    # Generate true state and data
    rho_true <- random_density_fullrank(N, seed = i * 100)

    # Simulate some data
    a_seq <- sample(1:lib$k, 50, replace = TRUE)
    prob_list <- born_probs_list(rho_true, lib$Q_list)
    b_seq <- sapply(a_seq, function(a) {
      sample.int(length(prob_list[[a]]), 1, prob = prob_list[[a]])
    })

    Nab <- counts_from_ab(a_seq, b_seq, lib$ab_row, lib$M)

    # Fit MLE with eigenvalue floor
    fit <- fit_theta_cvxr(N, sigmas, sc$S_ab, sc$c_ab, Nab,
                          eta = eta, solver = "SCS")

    if (fit$status %in% c("optimal", "optimal_inaccurate")) {
      rho_hat <- fit$rho_hat
      eigs <- eigen(hermitianize(rho_hat), symmetric = FALSE, only.values = TRUE)$values
      min_eig <- min(Re(eigs))

      # Check min eigenvalue >= eta - tolerance
      if (min_eig < eta - tol) {
        failures <- c(failures, list(list(
          test = i, min_eig = min_eig, eta = eta
        )))
      }
    }
  }

  list(
    n_tests = n_tests,
    eta = eta,
    n_failures = length(failures),
    valid = length(failures) == 0,
    failures = failures
  )
}

## ==========================================================================
## MEASUREMENT LIBRARY CHECKS
## ==========================================================================

#' Validate POVM completeness: Σ_b Q_{a,b} = I
#' @param lib Measurement library
#' @param tol Numerical tolerance
#' @return Validation results
validate_povm_completeness <- function(lib, tol = 1e-10) {
  N <- lib$N
  I_N <- diag(N) + 0i

  results <- list()
  all_valid <- TRUE

  for (a in 1:lib$k) {
    Q_sum <- Reduce(`+`, lib$Q_list[[a]])
    err <- max(Mod(Q_sum - I_N))

    results[[lib$setting_labels[a]]] <- list(
      valid = err < tol,
      error = err
    )

    if (err >= tol) all_valid <- FALSE
  }

  list(valid = all_valid, setting_results = results)
}

#' Validate POVM effects are PSD
#' @param lib Measurement library
#' @param tol Numerical tolerance
#' @return Validation results
validate_povm_psd <- function(lib, tol = 1e-10) {
  results <- list()
  all_valid <- TRUE

  for (a in 1:lib$k) {
    for (b in seq_along(lib$Q_list[[a]])) {
      Q <- lib$Q_list[[a]][[b]]
      eigs <- eigen(hermitianize(Q), symmetric = FALSE, only.values = TRUE)$values
      min_eig <- min(Re(eigs))

      key <- sprintf("%s_b%d", lib$setting_labels[a], b)
      results[[key]] <- list(
        valid = min_eig >= -tol,
        min_eigenvalue = min_eig
      )

      if (min_eig < -tol) all_valid <- FALSE
    }
  }

  list(valid = all_valid, effect_results = results)
}

## ==========================================================================
## BORN PROBABILITY CHECKS
## ==========================================================================

#' Validate Born probabilities sum to 1 and are non-negative
#' @param lib Measurement library
#' @param n_tests Number of random state tests
#' @param tol Numerical tolerance
#' @return Validation results
validate_born_probabilities <- function(lib, n_tests = 50, tol = 1e-10) {
  N <- lib$N
  failures <- list()

  for (i in 1:n_tests) {
    rho <- random_density_fullrank(N, seed = i * 17)
    prob_list <- born_probs_list(rho, lib$Q_list)

    for (a in 1:lib$k) {
      p <- prob_list[[a]]

      # Check non-negativity
      if (any(p < -tol)) {
        failures <- c(failures, list(list(
          test = i, setting = a, issue = "negative", min_p = min(p)
        )))
      }

      # Check sum to 1
      sum_err <- abs(sum(p) - 1)
      if (sum_err > tol) {
        failures <- c(failures, list(list(
          test = i, setting = a, issue = "sum", sum_error = sum_err
        )))
      }
    }
  }

  list(
    n_tests = n_tests,
    n_failures = length(failures),
    valid = length(failures) == 0,
    failures = failures
  )
}

## ==========================================================================
## FISHER INFORMATION CHECKS
## ==========================================================================

#' Validate Fisher information is PSD
#' @param lib Measurement library
#' @param sigmas Basis matrices
#' @param n_tests Number of tests
#' @param tol Tolerance
#' @return Validation results
validate_fisher_psd_all <- function(lib, sigmas, n_tests = 20, tol = 1e-10) {
  N <- lib$N
  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)

  failures <- list()

  for (i in 1:n_tests) {
    # Random theta (within valid region)
    theta <- rnorm(length(sigmas)) * 0.3

    for (a in 1:lib$k) {
      Ia <- fisher_info_setting(theta, a, sc$S_ab, sc$c_ab, lib$ab_df)
      eigs <- eigen(Ia, symmetric = TRUE, only.values = TRUE)$values
      min_eig <- min(Re(eigs))

      if (min_eig < -tol) {
        failures <- c(failures, list(list(
          test = i, setting = a, min_eigenvalue = min_eig
        )))
      }
    }
  }

  list(
    n_tests = n_tests,
    n_failures = length(failures),
    valid = length(failures) == 0,
    failures = failures
  )
}

#' Validate total FI monotonicity (adding samples increases FI)
#' @param lib Measurement library
#' @param sigmas Basis matrices
#' @param tol Tolerance
#' @return Validation results
validate_fi_monotonicity <- function(lib, sigmas, tol = 1e-10) {
  N <- lib$N
  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)

  theta <- rnorm(length(sigmas)) * 0.3

  # Compute all per-setting FI
  I_list <- fisher_info_all_settings(theta, sc$S_ab, sc$c_ab, lib$ab_df)

  # Check each Ia is PSD
  all_psd <- TRUE
  for (a in names(I_list)) {
    Ia <- I_list[[a]]
    min_eig <- min(eigen(Ia, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig < -tol) all_psd <- FALSE
  }

  # J_n + I_a should be >= J_n (in PSD sense)
  # This is automatic since Ia >= 0

  list(
    valid = all_psd,
    message = if(all_psd) "All Ia are PSD, monotonicity guaranteed"
              else "Some Ia have negative eigenvalues"
  )
}

## ==========================================================================
## BURES METRIC CHECKS
## ==========================================================================

#' Validate Bures metric computation
#' @param sigmas Basis matrices
#' @param N Hilbert space dimension
#' @param n_tests Number of tests
#' @param tol Tolerance
#' @return Validation results
validate_bures_metric_all <- function(sigmas, N, n_tests = 20, tol = 1e-10) {
  failures <- list()

  for (i in 1:n_tests) {
    theta <- rnorm(length(sigmas)) * 0.3
    G_B <- metric_bures(theta, sigmas, N)

    # Check symmetry
    sym_err <- max(abs(G_B - t(G_B)))
    if (sym_err > tol) {
      failures <- c(failures, list(list(
        test = i, issue = "symmetry", error = sym_err
      )))
    }

    # Check PSD
    eigs <- eigen(G_B, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(Re(eigs))
    if (min_eig < -tol) {
      failures <- c(failures, list(list(
        test = i, issue = "PSD", min_eigenvalue = min_eig
      )))
    }
  }

  list(
    n_tests = n_tests,
    n_failures = length(failures),
    valid = length(failures) == 0,
    failures = failures
  )
}

#' Validate Lyapunov inverse: Ω_ρ(Ω_ρ^{-1}(Y)) ≈ Y
#' @param N Hilbert space dimension
#' @param n_tests Number of tests
#' @param tol Tolerance
#' @return Validation results
validate_lyapunov_inverse <- function(N, n_tests = 20, tol = 1e-10) {
  failures <- list()

  for (i in 1:n_tests) {
    # Random full-rank state
    rho <- random_density_fullrank(N, seed = i * 23)

    # Random Hermitian Y
    Y_raw <- matrix(rnorm(N*N) + 1i*rnorm(N*N), N, N)
    Y <- hermitianize(Y_raw)

    # Compute Ω^{-1}(Y)
    X <- lyapunov_inv(rho, Y)

    # Check Ω(X) ≈ Y
    Omega_X <- 0.5 * (rho %*% X + X %*% rho)
    err <- max(Mod(Omega_X - Y))

    if (err > tol) {
      failures <- c(failures, list(list(
        test = i, error = err
      )))
    }
  }

  list(
    n_tests = n_tests,
    n_failures = length(failures),
    valid = length(failures) == 0,
    failures = failures
  )
}

## ==========================================================================
## MUB CHECKS (Library B)
## ==========================================================================

#' Full MUB validation for Library B
#' @param tol Tolerance
#' @return Validation results
validate_library_B_full <- function(tol = 1e-10) {
  LB <- build_library_2q_B()

  results <- list()

  # 1. POVM completeness
  results$completeness <- validate_povm_completeness(LB, tol)

  # 2. Effects PSD
  results$psd <- validate_povm_psd(LB, tol)

  # 3. MUB unbiasedness
  results$mub <- validate_mub_unbiasedness(LB, tol)

  # 4. Basis orthonormality
  results$ortho <- validate_basis_orthonormality_2q(LB, tol)

  # Overall
  results$valid <- all(c(
    results$completeness$valid,
    results$psd$valid,
    results$mub$valid,
    results$ortho$valid
  ))

  results
}

## ==========================================================================
## SELECTOR SANITY CHECKS
## ==========================================================================
#' Validate selector consistency: exact and GI1 should give finite scores
#' @param lib Measurement library
#' @param sigmas Basis matrices
#' @param n_tests Number of tests
#' @return Validation results
validate_selector_sanity <- function(lib, sigmas, n_tests = 10) {
  N <- lib$N
  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)

  failures <- list()

  for (i in 1:n_tests) {
    # Generate some data
    rho_true <- random_density_fullrank(N, seed = i * 31)
    a_seq <- rep(1:lib$k, each = 10)  # Ensure coverage
    prob_list <- born_probs_list(rho_true, lib$Q_list)
    b_seq <- sapply(a_seq, function(a) {
      sample.int(length(prob_list[[a]]), 1, prob = prob_list[[a]])
    })

    Nab <- counts_from_ab(a_seq, b_seq, lib$ab_row, lib$M)
    counts <- counts_by_setting(Nab, lib$ab_df)

    # Fit MLE
    fit <- fit_theta_cvxr(N, sigmas, sc$S_ab, sc$c_ab, Nab,
                          eta = 1e-3, solver = "SCS")

    if (!fit$status %in% c("optimal", "optimal_inaccurate")) {
      next
    }

    theta_hat <- fit$theta_hat

    # Frobenius metric
    G <- metric_frobenius(sigmas)

    # Test exact selector
    sel_exact <- select_next_setting_metric(
      theta_hat, counts, sc$S_ab, sc$c_ab, lib$ab_df,
      G, method = "exact", ridge = 1e-8
    )

    if (!all(is.finite(sel_exact$scores))) {
      failures <- c(failures, list(list(
        test = i, method = "exact", issue = "non-finite scores"
      )))
    }

    # Test GI1 selector
    sel_gi1 <- select_next_setting_metric(
      theta_hat, counts, sc$S_ab, sc$c_ab, lib$ab_df,
      G, method = "GI1", ridge = 1e-8
    )

    if (!all(is.finite(sel_gi1$scores))) {
      failures <- c(failures, list(list(
        test = i, method = "GI1", issue = "non-finite scores"
      )))
    }
  }

  list(
    n_tests = n_tests,
    n_failures = length(failures),
    valid = length(failures) == 0,
    failures = failures
  )
}

## ==========================================================================
## COMPREHENSIVE VALIDATION RUNNER
## ==========================================================================

#' Run all validation checks
#' @param verbose Print results
#' @return Summary of all validations
run_all_validations <- function(verbose = TRUE) {
  results <- list()

  if (verbose) cat("=== Running Validation Suite ===\n\n")

  # ===== One-qubit validations =====
  if (verbose) cat("--- One Qubit ---\n")

  pauli <- build_pauli_basis_1q()
  sigmas_1q <- pauli$sigmas
  L1 <- build_library_1q_L1()

  # State parametrization
  if (verbose) cat("rho_of_theta: ")
  results$rho_1q <- validate_rho_of_theta(sigmas_1q, N = 2)
  if (verbose) cat(if(results$rho_1q$valid) "PASS\n" else "FAIL\n")

  # POVM checks
  if (verbose) cat("POVM completeness: ")
  results$povm_1q <- validate_povm_completeness(L1)
  if (verbose) cat(if(results$povm_1q$valid) "PASS\n" else "FAIL\n")

  if (verbose) cat("POVM PSD: ")
  results$povm_psd_1q <- validate_povm_psd(L1)
  if (verbose) cat(if(results$povm_psd_1q$valid) "PASS\n" else "FAIL\n")

  # Born probabilities
  if (verbose) cat("Born probabilities: ")
  results$born_1q <- validate_born_probabilities(L1)
  if (verbose) cat(if(results$born_1q$valid) "PASS\n" else "FAIL\n")

  # Fisher info
  if (verbose) cat("Fisher info PSD: ")
  results$fisher_1q <- validate_fisher_psd_all(L1, sigmas_1q)
  if (verbose) cat(if(results$fisher_1q$valid) "PASS\n" else "FAIL\n")

  # Bures metric
  if (verbose) cat("Bures metric: ")
  results$bures_1q <- validate_bures_metric_all(sigmas_1q, N = 2)
  if (verbose) cat(if(results$bures_1q$valid) "PASS\n" else "FAIL\n")

  # Lyapunov inverse
  if (verbose) cat("Lyapunov inverse: ")
  results$lyap_1q <- validate_lyapunov_inverse(N = 2)
  if (verbose) cat(if(results$lyap_1q$valid) "PASS\n" else "FAIL\n")

  # Selector sanity
  if (verbose) cat("Selector sanity: ")
  results$selector_1q <- validate_selector_sanity(L1, sigmas_1q)
  if (verbose) cat(if(results$selector_1q$valid) "PASS\n" else "FAIL\n")

  # ===== Two-qubit validations =====
  if (verbose) cat("\n--- Two Qubits ---\n")

  basis_2q <- build_pauli_product_basis_2q()
  sigmas_2q <- basis_2q$sigmas
  LA <- build_library_2q_A()
  LB <- build_library_2q_B()

  # State parametrization
  if (verbose) cat("rho_of_theta: ")
  results$rho_2q <- validate_rho_of_theta(sigmas_2q, N = 4)
  if (verbose) cat(if(results$rho_2q$valid) "PASS\n" else "FAIL\n")

  # Library A POVM
  if (verbose) cat("Library A POVM: ")
  results$povm_A <- validate_povm_completeness(LA)
  if (verbose) cat(if(results$povm_A$valid) "PASS\n" else "FAIL\n")

  # Library B full
  if (verbose) cat("Library B (MUB): ")
  results$lib_B <- validate_library_B_full()
  if (verbose) cat(if(results$lib_B$valid) "PASS\n" else "FAIL\n")

  # Fisher info (Library A)
  if (verbose) cat("Fisher info PSD (A): ")
  results$fisher_A <- validate_fisher_psd_all(LA, sigmas_2q, n_tests = 10)
  if (verbose) cat(if(results$fisher_A$valid) "PASS\n" else "FAIL\n")

  # Bures metric
  if (verbose) cat("Bures metric: ")
  results$bures_2q <- validate_bures_metric_all(sigmas_2q, N = 4, n_tests = 10)
  if (verbose) cat(if(results$bures_2q$valid) "PASS\n" else "FAIL\n")

  # Summary
  all_valid <- all(sapply(results, function(r) {
    if (is.list(r) && "valid" %in% names(r)) r$valid else TRUE
  }))

  if (verbose) {
    cat("\n=== Summary ===\n")
    cat(sprintf("All validations: %s\n", if(all_valid) "PASS" else "FAIL"))
  }

  results$all_valid <- all_valid
  results
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Run all validations
# val_results <- run_all_validations(verbose = TRUE)

cat("09_validation.R loaded successfully\n")
