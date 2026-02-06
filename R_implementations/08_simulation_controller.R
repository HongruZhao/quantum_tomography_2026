## ============================================================================
## 08_simulation_controller.R
## Simulation controller: orchestrates experiments and generates 15 plots
## Based on: 08_simulation_controller.md
## ============================================================================

source("07_adaptive_design.R")

## ==========================================================================
## SIMULATION CONFIGURATION
## ==========================================================================

#' Default simulation configuration
#' @return List of simulation parameters
default_sim_config <- function() {
  list(
    # Monte Carlo parameters
    n_rep = 100,           # Number of replicates
    seed_base = 1,         # Base RNG seed
    case_seeds = NULL,     # Optional named list mapping case_id -> seed (replicate seeding)
    system_seeds = NULL,   # Optional named list mapping system ("1q","2q") -> seed (true state)
    true_state_alpha = 0.99,# Mixing weight for pure-state vs maximally mixed

    # Sample size
    n_total = 400,         # Total shots per replicate
    check_every = 1,       # Recompute MLE every N steps

    # Stabilization parameters
    eta_mle = 1e-3,        # Eigenvalue floor for MLE
    ridge_init = 1e-6,     # Ridge for early selection

    # CVXR options
    solver = "SCS",
    eps_log = 1e-12,
    verbose_solver = FALSE
  )
}

## ==========================================================================
## EXPERIMENT REGISTRY
## ==========================================================================

#' Get all one-qubit libraries
#' @return List of library constructors and metadata
get_1q_libraries <- function() {
  list(
    L1 = list(
      name = "L1",
      full_name = "Pauli PVMs (X,Y,Z)",
      builder = build_library_1q_L1
    ),
    L2 = list(
      name = "L2",
      full_name = "Nine-axis PVMs",
      builder = build_library_1q_L2
    ),
    L3 = list(
      name = "L3",
      full_name = "Random 4-axis PVMs",
      builder = function() build_library_1q_L3(seed = 42)
    )
  )
}

#' Get all two-qubit libraries
#' @return List of library constructors and metadata
get_2q_libraries <- function() {
  list(
    A = list(
      name = "A",
      full_name = "15 Pauli-parity PVMs",
      builder = build_library_2q_A
    ),
    B = list(
      name = "B",
      full_name = "5 MUB bases",
      builder = build_library_2q_B
    )
  )
}

#' Get loss metrics for a given basis
#' @param sigmas List of basis matrices
#' @param N Hilbert space dimension
#' @param Obs List of observables (optional)
#' @return List of metric function factories and metadata
get_loss_metrics <- function(sigmas, N) {
  list(
    frobenius = list(
      name = "frobenius",
      full_name = "Frobenius (Hilbert-Schmidt)",
      make_metric = function(lib) make_metric_frobenius(sigmas)
    ),
    bures = list(
      name = "bures",
      full_name = "Bures (SLD)",
      make_metric = function(lib) make_metric_bures(sigmas, N)
    ),
    observable = list(
      name = "observable",
      full_name = "Observable-expectation",
      make_metric = function(lib) {
        if (is.null(lib$Obs)) stop("Library is missing Obs list for observable metric.")
        make_metric_observable(sigmas, lib$Obs)
      }
    )
  )
}

## ==========================================================================
## SINGLE EXPERIMENT RUNNER
## ==========================================================================

#' Run single experiment: one (system, library, loss) combination
#' @param system "1q" or "2q"
#' @param lib_info Library info from registry
#' @param loss_info Loss info from registry
#' @param sim_cfg Simulation configuration
#' @param verbose Print progress
#' @return List with results for all three policies
run_single_experiment <- function(system, lib_info, loss_info, sim_cfg,
                                  case_seed = NULL, verbose = TRUE) {

  if (verbose) {
    cat(sprintf("\n=== Running: %s / %s / %s ===\n",
                system, lib_info$full_name, loss_info$full_name))
  }
  case_id <- sprintf("%s_%s_%s", system, lib_info$name, loss_info$name)

  # Build library
  lib <- lib_info$builder()
  N <- lib$N

  # Build basis
  if (N == 2) {
    basis <- build_pauli_basis_1q()
  } else {
    basis <- build_pauli_product_basis_2q()
  }
  sigmas <- basis$sigmas

  # Build metric function
  metric_fun <- loss_info$make_metric(lib)
  # True state seed is shared across all cases within the same system ("1q"/"2q").
  rho_base_seed <- sim_cfg$seed_base
  if (!is.null(sim_cfg$system_seeds) &&
      !is.null(sim_cfg$system_seeds[[system]])) {
    rho_base_seed <- as.integer(sim_cfg$system_seeds[[system]])
  }
  if (system == "1q") {
    rho_mix <- fixed_pure_state_1q_0()
    rho_pure <- fixed_pure_state_1q_yphase(phi = 0.2)
    rho_true_fixed <- mix_with_state(rho_pure, rho_mix, sim_cfg$true_state_alpha)
  } else {
    rho_pure <- fixed_pure_state_2q_phi_phase(phi = 0.2)
    rho_true_fixed <- mix_with_maximally_mixed(rho_pure, sim_cfg$true_state_alpha)
  }

  # Run all three policies
  policies <- c("uniform", "exact", "GI1")
  results <- list()

  for (policy in policies) {
    if (verbose) cat(sprintf("  Policy: %s\n", policy))

    # Monte Carlo
    policy_seed_key <- paste0(case_id, "_", policy)
    rep_base_seed <- sim_cfg$seed_base
    if (!is.null(sim_cfg$case_seeds) &&
        !is.null(sim_cfg$case_seeds[[case_id]])) {
      rep_base_seed <- as.integer(sim_cfg$case_seeds[[case_id]])
    } else if (!is.null(case_seed)) {
      rep_base_seed <- as.integer(case_seed)
    }
    risk_mat <- run_monte_carlo(
      lib = lib,
      sigmas = sigmas,
      n_total = sim_cfg$n_total,
      metric_fun = metric_fun,
      policy = policy,
      lib_name = lib_info$name,
      n_rep = sim_cfg$n_rep,
      sim_cfg = sim_cfg,
      rho_true = rho_true_fixed,
      seed_base = rep_base_seed,
      seed_key = policy_seed_key,
      verbose = FALSE
    )

    # Compute mean and SE
    mean_risk <- colMeans(risk_mat, na.rm = TRUE)
    sd_risk <- apply(risk_mat, 2, sd, na.rm = TRUE)
    se_risk <- sd_risk / sqrt(sim_cfg$n_rep)

    results[[policy]] <- list(
      policy = policy,
      risk_matrix = risk_mat,
      mean_risk = mean_risk,
      sd_risk = sd_risk,
      se_risk = se_risk
    )
  }

  list(
    case_id = case_id,
    system = system,
    library = lib_info$name,
    library_full = lib_info$full_name,
    loss = loss_info$name,
    loss_full = loss_info$full_name,
    n_total = sim_cfg$n_total,
    n_rep = sim_cfg$n_rep,
    results = results
  )
}

## ==========================================================================
## PLOTTING FUNCTIONS
## ==========================================================================

#' Plot risk curves for a single experiment
#' @param exp_result Result from run_single_experiment
#' @param filename Output filename (optional)
#' @param ... Additional plot parameters
plot_risk_curves <- function(exp_result, filename = NULL,
                             x_start = NULL, x_end = NULL, ...) {

  n_total <- exp_result$n_total
  if (is.null(x_start)) x_start <- min(20, n_total)
  if (is.null(x_end)) x_end <- min(100, n_total)
  n_seq <- x_start:x_end

  # Extract mean curves
  uniform_risk <- exp_result$results$uniform$mean_risk[x_start:x_end]
  exact_risk <- exp_result$results$exact$mean_risk[x_start:x_end]
  gi1_risk <- exp_result$results$GI1$mean_risk[x_start:x_end]

  # Determine y range
  all_risk <- c(uniform_risk, exact_risk, gi1_risk)
  y_range <- range(all_risk[is.finite(all_risk)], na.rm = TRUE)

  # Open file if specified
  if (!is.null(filename)) {
    png(filename, width = 800, height = 600)
  }

  # Plot
  plot(n_seq, uniform_risk, type = "l", col = "gray50", lwd = 2,
       xlim = c(x_start, x_end), ylim = y_range,
       xlab = "Number of samples (n)",
       ylab = expression("Scaled Proxy Risk " * n %.% R[n](G)),
       main = sprintf("%s / %s / %s",
                      ifelse(exp_result$system == "1q", "One qubit", "Two qubits"),
                      exp_result$library_full,
                      exp_result$loss_full),
       ...)
  lines(n_seq, exact_risk, col = "blue", lwd = 2)
  lines(n_seq, gi1_risk, col = "red", lwd = 2, lty = 2)

  legend("topright",
         legend = c("Uniform", "Exact", "GI1"),
         col = c("gray50", "blue", "red"),
         lwd = 2, lty = c(1, 1, 2),
         bty = "n")

  grid()

  if (!is.null(filename)) {
    dev.off()
    cat(sprintf("Saved: %s\n", filename))
  }
}

#' Create tidy data frame from experiment result
#' @param exp_result Result from run_single_experiment
#' @return Data frame in tidy format
experiment_to_df <- function(exp_result) {
  n_total <- exp_result$n_total

  df_list <- lapply(names(exp_result$results), function(policy) {
    r <- exp_result$results[[policy]]
    data.frame(
      system = exp_result$system,
      library = exp_result$library,
      loss = exp_result$loss,
      policy = policy,
      n = 1:n_total,
      risk_mean = r$mean_risk,
      risk_sd = r$sd_risk,
      risk_se = r$se_risk
    )
  })

  do.call(rbind, df_list)
}

## ==========================================================================
## FULL SIMULATION RUNNER (15 EXPERIMENTS)
## ==========================================================================

#' Run all 15 experiments: (3+2 libraries) × 3 losses
#' @param sim_cfg Simulation configuration
#' @param output_dir Directory for output files
#' @param verbose Print progress
#' @return List with all results
run_all_experiments <- function(sim_cfg = NULL, output_dir = "results",
                                verbose = TRUE) {

  if (is.null(sim_cfg)) sim_cfg <- default_sim_config()

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }

  all_results <- list()
  all_dfs <- list()
  exp_count <- 0

  # ===== One-qubit experiments =====
  libs_1q <- get_1q_libraries()
  pauli <- build_pauli_basis_1q()
  losses_1q <- get_loss_metrics(pauli$sigmas, N = 2)

  for (lib_name in names(libs_1q)) {
    lib_info <- libs_1q[[lib_name]]

    for (loss_name in names(losses_1q)) {
      loss_info <- losses_1q[[loss_name]]
      exp_count <- exp_count + 1

      if (verbose) {
        cat(sprintf("\n[%d/15] 1q / %s / %s\n", exp_count, lib_name, loss_name))
      }

      result <- run_single_experiment(
        system = "1q",
        lib_info = lib_info,
        loss_info = loss_info,
        sim_cfg = sim_cfg,
        verbose = verbose
      )

      # Store
      key <- sprintf("1q_%s_%s", lib_name, loss_name)
      all_results[[key]] <- result
      all_dfs[[key]] <- experiment_to_df(result)

      # Plot
      plot_file <- file.path(plots_dir, sprintf("%s.png", key))
      plot_risk_curves(result, filename = plot_file)
    }
  }

  # ===== Two-qubit experiments =====
  libs_2q <- get_2q_libraries()
  basis_2q <- build_pauli_product_basis_2q()
  losses_2q <- get_loss_metrics(basis_2q$sigmas, N = 4)

  for (lib_name in names(libs_2q)) {
    lib_info <- libs_2q[[lib_name]]

    for (loss_name in names(losses_2q)) {
      loss_info <- losses_2q[[loss_name]]
      exp_count <- exp_count + 1

      if (verbose) {
        cat(sprintf("\n[%d/15] 2q / %s / %s\n", exp_count, lib_name, loss_name))
      }

      result <- run_single_experiment(
        system = "2q",
        lib_info = lib_info,
        loss_info = loss_info,
        sim_cfg = sim_cfg,
        verbose = verbose
      )

      # Store
      key <- sprintf("2q_%s_%s", lib_name, loss_name)
      all_results[[key]] <- result
      all_dfs[[key]] <- experiment_to_df(result)

      # Plot
      plot_file <- file.path(plots_dir, sprintf("%s.png", key))
      plot_risk_curves(result, filename = plot_file)
    }
  }

  # Combine all data frames
  combined_df <- do.call(rbind, all_dfs)

  # Save results
  results_file <- file.path(output_dir, "results.csv")
  write.csv(combined_df, results_file, row.names = FALSE)
  cat(sprintf("\nSaved results to: %s\n", results_file))

  # Save config
  config_file <- file.path(output_dir, "config.rds")
  saveRDS(sim_cfg, config_file)

  list(
    results = all_results,
    combined_df = combined_df,
    config = sim_cfg,
    output_dir = output_dir
  )
}

## ==========================================================================
## QUICK TEST RUNNER
## ==========================================================================

#' Run a quick test with reduced parameters
#' @param n_rep Number of replicates
#' @param n_total Samples per replicate
#' @return Test results
run_quick_test <- function(n_rep = 10, n_total = 100) {
  sim_cfg <- default_sim_config()
  sim_cfg$n_rep <- n_rep
  sim_cfg$n_total <- n_total
  sim_cfg$check_every <- 5  # Less frequent MLE updates

  cat("Running quick test with reduced parameters...\n")
  cat(sprintf("  n_rep = %d, n_total = %d\n", n_rep, n_total))

  # Just run one experiment as test
  pauli <- build_pauli_basis_1q()
  lib_info <- list(name = "L1", full_name = "Pauli PVMs", builder = build_library_1q_L1)
  loss_info <- list(name = "frobenius", full_name = "Frobenius",
                    make_metric = function() make_metric_frobenius(pauli$sigmas))

  result <- run_single_experiment(
    system = "1q",
    lib_info = lib_info,
    loss_info = loss_info,
    sim_cfg = sim_cfg,
    verbose = TRUE
  )

  # Plot
  plot_risk_curves(result)

  result
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Example: run quick test
# test_result <- run_quick_test(n_rep = 5, n_total = 50)

# Example: run full simulation (takes a long time!)
# full_results <- run_all_experiments(verbose = TRUE)

cat("08_simulation_controller.R loaded successfully\n")
