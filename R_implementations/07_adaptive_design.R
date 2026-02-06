## ============================================================================
## 07_adaptive_design.R
## Adaptive design policies: Uniform, Exact one-step, GI1 (first-order proxy)
## Based on: 07_adaptive_design_policies.md
## ============================================================================

source("06_fisher_and_metrics.R")

## ==========================================================================
## METRIC-WEIGHTED A-OPTIMAL SELECTION
## ==========================================================================

#' Select next measurement setting using metric-weighted A-optimal criterion
#' Exact: argmin_a tr(G (J_n + I_a)^{-1})
#' GI1:   argmax_a tr(G J_n^{-1} I_a J_n^{-1})
#' @param theta_hat Current Bloch estimate
#' @param counts Vector of counts per setting
#' @param S_ab Coefficient matrix
#' @param c_ab Constant vector
#' @param ab_df Data frame with (a, b, row) mapping
#' @param G Metric matrix (d x d)
#' @param method Selection method ("exact" or "GI1")
#' @param ridge Ridge parameter for numerical stability
#' @param eps Tolerance for Fisher info computation
#' @return List with a_next (selected setting), scores, method_used
select_next_setting_metric <- function(theta_hat, counts, S_ab, c_ab, ab_df,
                                       G, method = c("exact", "GI1"),
                                       ridge = 1e-8, eps = 1e-12) {
  method <- match.arg(method)
  d <- ncol(S_ab)

  # Compute total Fisher info and per-setting FI
  tf <- total_fisher_info(theta_hat, counts, S_ab, c_ab, ab_df, eps)
  J_n <- tf$J_total
  I_list <- tf$I_list
  k <- length(I_list)

  # Regularize J_n if needed
  J_reg <- J_n + ridge * diag(d)

  scores <- numeric(k)
  names(scores) <- names(I_list)

  if (method == "exact") {
    # Exact: score[a] = tr(G (J_n + I_a)^{-1})  -> minimize
    for (i in seq_along(I_list)) {
      Ia <- I_list[[i]]
      J_cand <- J_reg + Ia
      J_cand_inv <- tryCatch(
        solve(J_cand),
        # J_cand already includes ridge via J_reg, avoid double-adding.
        error = function(e) safe_solve(J_cand, ridge = 0)
      )
      scores[i] <- sum(diag(G %*% J_cand_inv))
    }
    a_next <- as.integer(names(scores)[which.min(scores)])
    method_used <- "exact"

  } else {  # GI1
    # GI1: score[a] = tr(G J_n^{-1} I_a J_n^{-1})  -> maximize
    J_inv <- tryCatch(
      solve(J_reg),
      # J_reg already includes ridge, avoid double-adding.
      error = function(e) safe_solve(J_reg, ridge = 0)
    )

    for (i in seq_along(I_list)) {
      Ia <- I_list[[i]]
      # tr(G J^{-1} I_a J^{-1}) = tr(J^{-1} G J^{-1} I_a) for symmetric G
      GJinv <- G %*% J_inv
      scores[i] <- sum(diag(GJinv %*% Ia %*% J_inv))
    }
    a_next <- as.integer(names(scores)[which.max(scores)])
    method_used <- "GI1"
  }

  list(
    a_next = a_next,
    scores = scores,
    method_used = method_used,
    J_n = J_n,
    I_list = I_list
  )
}

## ==========================================================================
## FULL ADAPTIVE DESIGN SEQUENCE
## ==========================================================================

#' Get deterministic initialization settings for each library
#' @param lib_name Library name (L1, L2, L3, A, B)
#' @param k Total number of settings in library
#' @return Vector of setting indices to use for initialization
get_init_settings <- function(lib_name, k) {
  switch(lib_name,
    "L1" = 1:3,              # X, Y, Z (3 settings)
    "L2" = 1:3,              # X, Y, Z first 3 of 9 (3 settings)
    "L3" = 1:4,              # All 4 random axes (4 settings)
    "A"  = 1:15,             # All 15 Pauli-parity PVMs (15 settings)
    "B"  = 1:5,              # All 5 MUB bases (5 settings)
    1:min(k, 3)              # Default: first 3 settings
  )
}

#' Run full adaptive design sequence with metric-weighted selection
#' Uses DETERMINISTIC initialization: measure each init setting once
#' @param lib Measurement library object
#' @param sigmas List of basis matrices
#' @param rho_true True density matrix (for sampling)
#' @param n_total Total number of samples
#' @param metric_fun Function(theta) returning metric matrix G
#' @param policy Design policy ("uniform", "exact", "GI1")
#' @param lib_name Library name for initialization (L1, L2, L3, A, B)
#' @param eta_mle Eigenvalue floor for MLE
#' @param ridge Ridge parameter for selection
#' @param solver CVXR solver
#' @param seed Random seed
#' @param verbose Print progress
#' @param check_every How often to recompute MLE (default 1 = every step)
#' @return List with full trajectory data
adaptive_design_sequence_metric <- function(
    lib, sigmas, rho_true, n_total,
    metric_fun,
    policy = c("uniform", "exact", "GI1"),
    lib_name = NULL,
    eta_mle = 1e-4,
    ridge = 1e-8,
    solver = "SCS",
    seed = NULL,
    verbose = FALSE,
    check_every = 1
) {
  policy <- match.arg(policy)
  if (!is.null(seed)) set.seed(seed)

  N <- lib$N
  d <- length(sigmas)
  k <- lib$k

  # Build affine model
  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)
  S_ab <- sc$S_ab
  c_ab <- sc$c_ab

  # Pre-compute Born probabilities
  prob_list <- born_probs_list(rho_true, lib$Q_list)

  # Compute TRUE theta from rho_true for oracle risk evaluation (plots)
  theta_true <- theta_from_rho(rho_true, sigmas)
  G_true <- metric_fun(theta_true)

  # Get deterministic initialization settings
  init_settings <- get_init_settings(lib_name, k)
  n_init <- length(init_settings)

  # Initialize storage
  a_seq <- integer(n_total)
  b_seq <- integer(n_total)
  Nab <- numeric(lib$M)
  risk_history <- numeric(n_total)        # n * R_n(G) using TRUE theta
  theta_history <- matrix(NA, n_total, d)

  theta_hat <- rep(0, d)

  for (t in 1:n_total) {
    # DETERMINISTIC INITIALIZATION: use fixed settings for first n_init samples
    if (t <= n_init) {
      # Initialization phase: cycle through init_settings
      a_t <- init_settings[t]
    } else if (policy == "uniform") {
      # Uniform policy: random selection after initialization
      a_t <- sample.int(k, 1)
    } else {
      # Adaptive selection using current MLE
      counts <- counts_by_setting(Nab, lib$ab_df)
      G_hat <- metric_fun(theta_hat)
      ridge_sel <- if (t <= n_init) ridge else 0
      sel <- select_next_setting_metric(
        theta_hat, counts, S_ab, c_ab, lib$ab_df,
        G_hat, method = policy, ridge = ridge_sel
      )
      a_t <- sel$a_next
    }

    # Sample outcome from true state
    probs_a <- prob_list[[a_t]]
    b_t <- sample.int(length(probs_a), size = 1, prob = probs_a)

    # Update
    a_seq[t] <- a_t
    b_seq[t] <- b_t
    Nab <- increment_count(Nab, a_t, b_t, lib$ab_row)

    # Recompute MLE
    if (t %% check_every == 0 || t == n_total || t == n_init) {
      fit <- fit_theta_cvxr(N, sigmas, S_ab, c_ab, Nab,
                            eta = eta_mle, solver = solver, verbose = FALSE)

      if (fit$status %in% c("optimal", "optimal_inaccurate")) {
        theta_hat <- fit$theta_hat
      }
    }

    theta_history[t, ] <- theta_hat

    # Compute SCALED oracle risk n * R_n(G) at TRUE theta (for plotting)
    # Adaptive selection itself is based on MLE (theta_hat).
    counts_new <- counts_by_setting(Nab, lib$ab_df)
    tf_true <- total_fisher_info(theta_true, counts_new, S_ab, c_ab, lib$ab_df)
    ridge_eval <- if (t <= n_init) ridge else 0
    raw_risk <- proxy_risk(G_true, tf_true$J_total, ridge = ridge_eval)
    risk_history[t] <- t * raw_risk  # Scaled risk: n * R_n(G)

    if (verbose && t %% 50 == 0) {
      cat(sprintf("t=%d: scaled_risk=%.4f\n", t, risk_history[t]))
    }
  }

  list(
    policy = policy,
    n_total = n_total,
    n_init = n_init,
    init_settings = init_settings,
    a_seq = a_seq,
    b_seq = b_seq,
    Nab = Nab,
    theta_hat = theta_hat,
    theta_true = theta_true,
    theta_history = theta_history,
    risk_history = risk_history,
    counts_final = counts_by_setting(Nab, lib$ab_df)
  )
}

## ==========================================================================
## BATCH RUNNER FOR MONTE CARLO
## ==========================================================================

#' Run single replicate of adaptive design simulation
#' @param lib Measurement library object
#' @param sigmas List of basis matrices
#' @param n_total Total samples per replicate
#' @param metric_fun Function(theta) -> G
#' @param policy Design policy
#' @param lib_name Library name (L1, L2, L3, A, B) for initialization
#' @param sim_cfg Simulation configuration list
#' @param rep_seed Seed for this replicate
#' @return Risk curve (length n_total)
run_one_replicate <- function(lib, sigmas, n_total, metric_fun, policy,
                              lib_name, sim_cfg, rep_seed, rho_true = NULL) {
  N <- lib$N

  # Generate true state (or reuse fixed state if provided)
  if (is.null(rho_true)) {
    rho_true <- random_density_fullrank(N, eps_mix = 0.05, seed = rep_seed)
  }

  # Run adaptive sequence with deterministic initialization
  result <- adaptive_design_sequence_metric(
    lib = lib,
    sigmas = sigmas,
    rho_true = rho_true,
    n_total = n_total,
    metric_fun = metric_fun,
    policy = policy,
    lib_name = lib_name,
    eta_mle = sim_cfg$eta_mle,
    ridge = sim_cfg$ridge_init,
    solver = sim_cfg$solver,
    seed = rep_seed,
    verbose = FALSE,
    check_every = sim_cfg$check_every
  )

  result$risk_history
}

#' Run Monte Carlo simulation for one (library, loss, policy) combination
#' @param lib Measurement library
#' @param sigmas Basis matrices
#' @param n_total Samples per trajectory
#' @param metric_fun Metric function
#' @param policy Policy name
#' @param lib_name Library name (L1, L2, L3, A, B) for initialization
#' @param n_rep Number of replicates
#' @param sim_cfg Configuration
#' @param seed_base Base seed (used only if seed_key is NULL)
#' @param seed_key Stable string key for replicate seeding (order-independent)
#' @param verbose Print progress
#' @return Matrix of risk curves (n_rep x n_total)
run_monte_carlo <- function(lib, sigmas, n_total, metric_fun, policy,
                            lib_name, n_rep, sim_cfg, seed_base = 1,
                            rho_true = NULL, seed_key = NULL, verbose = FALSE) {
  risk_matrix <- matrix(NA, n_rep, n_total)

  for (r in 1:n_rep) {
    if (verbose && r %% 10 == 0) {
      cat(sprintf("  Replicate %d/%d\n", r, n_rep))
    }

    if (!is.null(seed_key)) {
      rep_seed <- stable_seed(paste0(seed_key, "_rep_", r), base = seed_base)
    } else {
      rep_seed <- seed_base + r * 100
    }
    risk_matrix[r, ] <- run_one_replicate(
      lib, sigmas, n_total, metric_fun, policy, lib_name, sim_cfg, rep_seed,
      rho_true = rho_true
    )
  }

  risk_matrix
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Example usage (uncomment to test)
# cat("\n=== Adaptive Design Test ===\n")
#
# # Setup
# pauli <- build_pauli_basis_1q()
# sigmas <- pauli$sigmas
# L1 <- build_library_1q_L1()
#
# # Generate true state
# rho_true <- random_density_fullrank(2, eps_mix = 0.1, seed = 123)
#
# # Metric functions
# metric_frob <- make_metric_frobenius(sigmas)
#
# # Run short adaptive sequence
# result <- adaptive_design_sequence_metric(
#   lib = L1,
#   sigmas = sigmas,
#   rho_true = rho_true,
#   n_total = 50,
#   metric_fun = metric_frob,
#   policy = "GI1",
#   eta_mle = 1e-4,
#   solver = "SCS",
#   seed = 456,
#   verbose = TRUE
# )
#
# cat("\nFinal risk:", tail(result$risk_history, 1), "\n")
# cat("Counts per setting:", result$counts_final, "\n")

cat("07_adaptive_design.R loaded successfully\n")
