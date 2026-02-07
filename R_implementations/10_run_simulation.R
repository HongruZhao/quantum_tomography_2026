## ============================================================================
## 10_run_simulation.R
## Run all 15 adaptive-design simulations and generate plots
## Follows quantum_tomography_tutorial.html and code_logic_overview.html
## ============================================================================

cat("==========================================================\n")
cat(" QUANTUM TOMOGRAPHY ADAPTIVE DESIGN SIMULATION\n")
cat("==========================================================\n\n")

## ==========================================================================
## 0. SETUP: Load all dependencies
## ==========================================================================

cat("Loading R implementation files...\n")

# Set working directory to R_implementations folder
script_dir <- tryCatch({
  # Prefer --file= argument when running via Rscript
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  # Fallback: if sourced interactively
  this_file <- sys.frame(1)$ofile
  if (!is.null(this_file) && this_file != "") {
    return(dirname(normalizePath(this_file)))
  }
  NULL
}, error = function(e) NULL)

if (!is.null(script_dir) && script_dir != "") {
  setwd(script_dir)
} else if (dir.exists("R_implementations")) {
  setwd("R_implementations")
}

# Source all implementation files in order
source("01_utilities.R")
source("02_state_basis.R")
source("03_measurement_library_1q.R")
source("04_measurement_library_2q.R")
source("05_mle_cvxr_solver.R")
source("06_fisher_and_metrics.R")
source("07_adaptive_design.R")
source("08_simulation_controller.R")

cat("All dependencies loaded successfully.\n\n")

## ==========================================================================
## 1. SIMULATION CONFIGURATION
## ==========================================================================

detect_parallel_workers <- function(default = 1L) {
  nc <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(nc) || nc < 1L) nc <- NA_integer_

  if (is.na(nc) || nc <= 1L) {
    out <- tryCatch(system("getconf _NPROCESSORS_ONLN", intern = TRUE), error = function(e) character(0))
    if (length(out) > 0) {
      val <- suppressWarnings(as.integer(out[1]))
      if (!is.na(val) && val > 1L) nc <- val
    }
  }

  if ((is.na(nc) || nc <= 1L) && identical(Sys.info()[["sysname"]], "Linux")) {
    out <- tryCatch(system("nproc", intern = TRUE), error = function(e) character(0))
    if (length(out) > 0) {
      val <- suppressWarnings(as.integer(out[1]))
      if (!is.na(val) && val > 1L) nc <- val
    }
  }

  if (is.na(nc) || nc < 1L) nc <- as.integer(default)
  as.integer(nc)
}

# Use default config and override key parameters as needed
simulation_config <- default_sim_config()

# Monte Carlo parameters
simulation_config$n_rep <- 1
simulation_config$seed_base <- 2026
simulation_config$case_seeds <- NULL  # optional: list("1q_L1_frobenius" = 2026, ...)
simulation_config$system_seeds <- list(
  "1q" = simulation_config$seed_base,
  "2q" = simulation_config$seed_base + 1
)
simulation_config$true_state_alpha <- 0.99

# Sample size and recompute cadence
simulation_config$n_total <- 400
simulation_config$check_every <- 1

# Stabilization parameters
simulation_config$eta_mle <- 1e-3
simulation_config$ridge_init <- 1e-8

# Solver
simulation_config$solver <- "SCS"
simulation_config$n_workers <- {
  detect_parallel_workers(default = 1L)
}

# Plot window (common across all cases)
plot_start <- 20
plot_end <- 400

cat("SIMULATION CONFIGURATION:\n")
cat(sprintf("  - Monte Carlo replicates: %d\n", simulation_config$n_rep))
cat(sprintf("  - Samples per trajectory: %d\n", simulation_config$n_total))
cat(sprintf("  - MLE recompute every: %d steps\n", simulation_config$check_every))
cat(sprintf("  - Eigenvalue floor (eta): %.2e\n", simulation_config$eta_mle))
cat(sprintf("  - CVXR solver: %s\n", simulation_config$solver))
cat(sprintf("  - Parallel workers: %d\n", simulation_config$n_workers))
cat(sprintf("  - True-state mix alpha: %.2f\n", simulation_config$true_state_alpha))
cat(sprintf("  - Plot window: n = %d..%d\n", plot_start, plot_end))
cat("\n")

## ==========================================================================
## 2. EXPERIMENT GRID (15 CASES)
## ==========================================================================

experiment_grid <- list(
  # One-qubit (N=2)
  list(system = "1q", lib = "L1", loss = "frobenius",  case_num = 1),
  list(system = "1q", lib = "L1", loss = "bures",      case_num = 2),
  list(system = "1q", lib = "L1", loss = "observable", case_num = 3),
  list(system = "1q", lib = "L2", loss = "frobenius",  case_num = 4),
  list(system = "1q", lib = "L2", loss = "bures",      case_num = 5),
  list(system = "1q", lib = "L2", loss = "observable", case_num = 6),
  list(system = "1q", lib = "L3", loss = "frobenius",  case_num = 7),
  list(system = "1q", lib = "L3", loss = "bures",      case_num = 8),
  list(system = "1q", lib = "L3", loss = "observable", case_num = 9),

  # Two-qubit (N=4)
  list(system = "2q", lib = "A",  loss = "frobenius",  case_num = 10),
  list(system = "2q", lib = "A",  loss = "bures",      case_num = 11),
  list(system = "2q", lib = "A",  loss = "observable", case_num = 12),
  list(system = "2q", lib = "B",  loss = "frobenius",  case_num = 13),
  list(system = "2q", lib = "B",  loss = "bures",      case_num = 14),
  list(system = "2q", lib = "B",  loss = "observable", case_num = 15)
)

## ==========================================================================
## 3. HELPERS
## ==========================================================================

get_library <- function(lib_name) {
  switch(lib_name,
    "L1" = build_library_1q_L1(),
    "L2" = build_library_1q_L2(),
    "L3" = build_library_1q_L3(seed = 42),
    "A"  = build_library_2q_A(),
    "B"  = build_library_2q_B(),
    stop(paste("Unknown library:", lib_name))
  )
}

get_basis <- function(system) {
  if (system == "1q") build_pauli_basis_1q() else build_pauli_product_basis_2q()
}

get_metric_function <- function(loss_name, sigmas, N, lib) {
  switch(loss_name,
    "frobenius"  = make_metric_frobenius(sigmas),
    "bures"      = make_metric_bures(sigmas, N),
    "observable" = {
      if (is.null(lib$Obs)) stop("Library is missing Obs list for observable metric.")
      make_metric_observable(sigmas, lib$Obs)
    },
    stop(paste("Unknown loss:", loss_name))
  )
}

# Build the fixed true state by system, matching the manuscript table setup:
# rho_star(alpha) = alpha * |psi><psi| + (1-alpha) * I/N.
build_true_state_for_system <- function(system, alpha) {
  if (system == "1q") {
    rho_pure <- fixed_pure_state_1q_yphase(phi = 0.2)
  } else if (system == "2q") {
    rho_pure <- fixed_pure_state_2q_phi_phase(phi = 0.2)
  } else {
    stop(paste("Unknown system:", system))
  }
  mix_with_maximally_mixed(rho_pure, alpha)
}

solve_optimal_design <- function(I_list, G, solver = "SCS") {
  k <- length(I_list)
  d <- nrow(G)

  pi <- CVXR::Variable(k)
  W <- CVXR::Variable(d, d, symmetric = TRUE)

  I_expr <- 0
  for (i in seq_len(k)) {
    I_expr <- I_expr + pi[i] * I_list[[i]]
  }
  I_expr <- (I_expr + t(I_expr)) / 2

  I_block <- CVXR::bmat(list(
    list(I_expr, diag(d)),
    list(diag(d), W)
  ))

  obj <- CVXR::Minimize(CVXR::sum_entries(G * W))
  constraints <- list(sum(pi) == 1, pi >= 0, I_block %>>% 0)

  prob <- CVXR::Problem(obj, constraints)
  result <- CVXR::solve(prob, solver = solver, verbose = FALSE, max_iters = 50000, eps = 1e-6)

  list(
    pi = result$getValue(pi),
    status = result$status
  )
}

solve_optimal_design_mirror <- function(I_list, G, ridge = 1e-8, max_iter = 5000,
                                        tol = 1e-10, eta = 1.0) {
  k <- length(I_list)
  d <- nrow(G)
  pi <- rep(1 / k, k)

  f_val <- Inf
  for (iter in 1:max_iter) {
    I_pi <- Reduce("+", lapply(seq_len(k), function(i) pi[i] * I_list[[i]]))
    I_pi <- I_pi + ridge * diag(d)
    J_inv <- solve(I_pi)
    f_val <- sum(diag(G %*% J_inv))

    grad <- numeric(k)
    for (i in seq_len(k)) {
      grad[i] <- -sum(diag(G %*% J_inv %*% I_list[[i]] %*% J_inv))
    }

    eta_i <- eta
    improved <- FALSE
    for (bt in 1:25) {
      pi_new <- pi * exp(-eta_i * grad)
      pi_new <- pi_new / sum(pi_new)
      I_pi_new <- Reduce("+", lapply(seq_len(k), function(i) pi_new[i] * I_list[[i]]))
      I_pi_new <- I_pi_new + ridge * diag(d)
      f_new <- sum(diag(G %*% solve(I_pi_new)))

      if (is.finite(f_new) && f_new <= f_val) {
        improved <- TRUE
        break
      }
      eta_i <- eta_i * 0.5
    }

    if (!improved) break

    if (max(abs(pi_new - pi)) < tol || abs(f_val - f_new) / max(1, abs(f_val)) < tol) {
      pi <- pi_new
      f_val <- f_new
      break
    }

    pi <- pi_new
    f_val <- f_new
  }

  list(pi = pi, value = f_val, iter = iter)
}

compute_plot_ylim <- function(risk_values, oracle_limit = NULL, floor_zero = TRUE) {
  valid_risks <- risk_values[is.finite(risk_values)]

  if (length(valid_risks) == 0) {
    y_min <- 0
    y_max <- 1
  } else {
    y_range <- range(valid_risks, na.rm = TRUE)
    span <- y_range[2] - y_range[1]
    y_pad <- if (is.finite(span) && span > 0) {
      0.05 * span
    } else {
      max(1e-6, abs(y_range[1]) * 0.05 + 1e-6)
    }
    y_min <- y_range[1] - y_pad
    y_max <- y_range[2] + y_pad
  }

  if (!is.null(oracle_limit) && is.finite(oracle_limit)) {
    oracle_eps <- max(1e-6, abs(oracle_limit) * 0.02)
    y_min <- min(y_min, oracle_limit - oracle_eps)
    y_max <- max(y_max, oracle_limit + oracle_eps)
  }

  if (floor_zero) y_min <- max(0, y_min)
  if (!is.finite(y_min) || !is.finite(y_max) || y_max <= y_min) {
    y_min <- if (is.finite(y_min)) y_min else 0
    y_max <- y_min + max(1e-6, abs(y_min) * 0.05 + 1e-6)
  }

  c(y_min, y_max)
}

plot_experiment <- function(exp_result, output_dir, plot_start, plot_end) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  n_total <- exp_result$n_total
  x_start <- min(plot_start, n_total)
  x_end <- min(plot_end, n_total)
  n_seq <- x_start:x_end

  uniform_risk <- exp_result$results$uniform$mean_risk[x_start:x_end]
  exact_risk <- exp_result$results$exact$mean_risk[x_start:x_end]
  gi1_risk <- exp_result$results$GI1$mean_risk[x_start:x_end]

  all_risks <- c(uniform_risk, exact_risk, gi1_risk)
  y_lim <- compute_plot_ylim(all_risks, oracle_limit = exp_result$oracle_limit, floor_zero = TRUE)

  system_name <- ifelse(exp_result$system == "1q", "One-Qubit", "Two-Qubit")
  lib_full <- switch(exp_result$library,
    "L1" = "L1: Pauli PVMs",
    "L2" = "L2: Nine-axis",
    "L3" = "L3: Random 4-axis",
    "A"  = "A: Pauli-parity",
    "B"  = "B: 5 MUBs"
  )
  loss_full <- switch(exp_result$loss,
    "frobenius"  = "Frobenius",
    "bures"      = "Bures",
    "observable" = "Observable"
  )

  filename <- file.path(output_dir, sprintf("case%02d_%s_%s_%s.png",
                                            exp_result$case_num,
                                            exp_result$system,
                                            exp_result$library,
                                            exp_result$loss))
  png(filename, width = 900, height = 650, res = 100)
  par(mar = c(5, 5, 4, 2))

  plot(n_seq, uniform_risk, type = "l", col = "gray40", lwd = 2.5,
       xlim = c(x_start, x_end), ylim = y_lim,
       xlab = "Number of Samples (n)",
       ylab = expression("Scaled Oracle Risk " * n %.% R[n](G)),
       main = sprintf("Case %d: %s / %s / %s", exp_result$case_num, system_name, lib_full, loss_full),
       cex.lab = 1.2, cex.main = 1.2)
  lines(n_seq, exact_risk, col = "blue", lwd = 2.5)
  lines(n_seq, gi1_risk, col = "red", lwd = 2.5, lty = 2)

  if (!is.null(exp_result$oracle_limit)) {
    abline(h = exp_result$oracle_limit, col = "black", lwd = 2, lty = 3)
    legend("topright",
           legend = c("Uniform", "Exact", "GI1", "Oracle limit"),
           col = c("gray40", "blue", "red", "black"),
           lwd = 2.5, lty = c(1, 1, 2, 3),
           bg = "white", box.lwd = 1)
  } else {
    legend("topright",
           legend = c("Uniform", "Exact", "GI1"),
           col = c("gray40", "blue", "red"),
           lwd = 2.5, lty = c(1, 1, 2),
           bg = "white", box.lwd = 1)
  }

  grid(lty = "dotted", col = "gray80")

  mtext(sprintf("Oracle risk at true theta | n_rep=%d", exp_result$n_rep),
        side = 1, line = 4, cex = 0.9, col = "gray50")

  dev.off()
}

## ==========================================================================
## 4. RUN SINGLE EXPERIMENT
## ==========================================================================

run_experiment <- function(exp, cfg) {
  case_id <- sprintf("%s_%s_%s", exp$system, exp$lib, exp$loss)
  system_str <- ifelse(exp$system == "1q", "1-Qubit", "2-Qubit")

  cat(sprintf("\n[Case %d/15] %s / %s / %s\n", exp$case_num, system_str, exp$lib, exp$loss))
  cat("----------------------------------------\n")

  lib <- get_library(exp$lib)
  basis <- get_basis(exp$system)
  sigmas <- basis$sigmas
  N <- lib$N
  metric_fun <- get_metric_function(exp$loss, sigmas, N, lib)
  rho_true_fixed <- build_true_state_for_system(exp$system, cfg$true_state_alpha)

  # Compute oracle limit for plotting: lim_{n} n * R_n(G) = tr(G I_pi*^{-1})
  theta_true <- theta_from_rho(rho_true_fixed, sigmas)
  G_true <- metric_fun(theta_true)
  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)
  I_list <- fisher_info_all_settings(theta_true, sc$S_ab, sc$c_ab, lib$ab_df)
  k <- length(I_list)
  I_uniform <- Reduce("+", I_list) / k
  uniform_limit <- proxy_risk(G_true, I_uniform)
  opt <- solve_optimal_design(I_list, G_true, solver = cfg$solver)
  pi_opt <- as.numeric(opt$pi)
  I_opt <- Reduce("+", lapply(seq_len(k), function(i) pi_opt[i] * I_list[[i]]))
  oracle_limit <- proxy_risk(G_true, I_opt)
  if (!is.finite(oracle_limit) || oracle_limit > uniform_limit * (1 + 1e-8)) {
    opt_md <- solve_optimal_design_mirror(I_list, G_true)
    pi_opt <- opt_md$pi
    I_opt <- Reduce("+", lapply(seq_len(k), function(i) pi_opt[i] * I_list[[i]]))
    oracle_limit <- proxy_risk(G_true, I_opt)
  }

  results <- list()
  for (policy in c("uniform", "exact", "GI1")) {
    cat(sprintf("  Running %s policy (%d replicates)...", toupper(policy), cfg$n_rep))
    flush.console()

    policy_seed_key <- paste0(case_id, "_", policy)
    rep_base_seed <- cfg$seed_base
    if (!is.null(cfg$case_seeds) && !is.null(cfg$case_seeds[[case_id]])) {
      rep_base_seed <- as.integer(cfg$case_seeds[[case_id]])
    }
    risk_matrix <- run_monte_carlo(
      lib = lib,
      sigmas = sigmas,
      n_total = cfg$n_total,
      metric_fun = metric_fun,
      policy = policy,
      lib_name = exp$lib,
      n_rep = cfg$n_rep,
      sim_cfg = cfg,
      rho_true = rho_true_fixed,
      seed_base = rep_base_seed,
      seed_key = policy_seed_key,
      verbose = FALSE
    )

    mean_risk <- colMeans(risk_matrix, na.rm = TRUE)
    sd_risk <- apply(risk_matrix, 2, sd, na.rm = TRUE)
    se_risk <- sd_risk / sqrt(cfg$n_rep)

    cat(sprintf(" Done! Final scaled risk: %.4f\n", tail(mean_risk, 1)))

    results[[policy]] <- list(
      policy = policy,
      risk_matrix = risk_matrix,
      mean_risk = mean_risk,
      sd_risk = sd_risk,
      se_risk = se_risk
    )
  }

  list(
    case_id = case_id,
    case_num = exp$case_num,
    system = exp$system,
    library = exp$lib,
    loss = exp$loss,
    N = N,
    n_total = cfg$n_total,
    n_rep = cfg$n_rep,
    oracle_limit = oracle_limit,
    results = results
  )
}

## ==========================================================================
## 5. RUN ALL EXPERIMENTS
## ==========================================================================

run_full_simulation <- function(cfg = simulation_config,
                                output_dir = "results",
                                plot_start = 20,
                                plot_end = 100) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  all_results <- list()
  all_data <- data.frame()

  n_workers <- if (!is.null(cfg$n_workers) && is.finite(cfg$n_workers)) {
    max(1L, as.integer(cfg$n_workers))
  } else {
    detect_parallel_workers(default = 1L)
  }
  n_workers <- min(n_workers, length(experiment_grid))

  if (.Platform$OS.type == "windows" && n_workers > 1L) {
    cat(sprintf("Parallel mclapply unavailable on Windows; falling back to sequential run.\n"))
    n_workers <- 1L
  }

  if (n_workers > 1L) {
    cat(sprintf("Running %d cases in parallel using %d workers...\n",
                length(experiment_grid), n_workers))
    exp_results <- parallel::mclapply(
      experiment_grid,
      function(exp) run_experiment(exp, cfg),
      mc.cores = n_workers,
      mc.preschedule = FALSE
    )
  } else {
    cat("Running cases sequentially...\n")
    exp_results <- lapply(experiment_grid, function(exp) run_experiment(exp, cfg))
  }

  for (exp_result in exp_results) {
    all_results[[exp_result$case_id]] <- exp_result

    for (policy in names(exp_result$results)) {
      r <- exp_result$results[[policy]]
      df_policy <- data.frame(
        case_num = exp_result$case_num,
        system = exp_result$system,
        library = exp_result$library,
        loss = exp_result$loss,
        policy = policy,
        n = 1:cfg$n_total,
        risk_mean = r$mean_risk,
        risk_sd = r$sd_risk,
        risk_se = r$se_risk
      )
      all_data <- rbind(all_data, df_policy)
    }

    plot_experiment(exp_result, plots_dir, plot_start, plot_end)
  }

  # Save outputs
  write.csv(all_data, file.path(output_dir, "results.csv"), row.names = FALSE)
  saveRDS(cfg, file.path(output_dir, "simulation_config.rds"))
  saveRDS(all_results, file.path(output_dir, "full_results.rds"))

  list(results = all_results, data = all_data, config = cfg, output_dir = output_dir)
}

## ==========================================================================
## 6. COMBINED 15-PANEL PLOT
## ==========================================================================

generate_combined_plot <- function(all_results, output_dir = "results",
                                   plot_start = 20, plot_end = 100) {
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  png(file.path(plots_dir, "all_15_cases_combined.png"),
      width = 2400, height = 2000, res = 150)
  par(mfrow = c(5, 3), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1))

  case_order <- c(
    "1q_L1_frobenius", "1q_L1_bures", "1q_L1_observable",
    "1q_L2_frobenius", "1q_L2_bures", "1q_L2_observable",
    "1q_L3_frobenius", "1q_L3_bures", "1q_L3_observable",
    "2q_A_frobenius",  "2q_A_bures",  "2q_A_observable",
    "2q_B_frobenius",  "2q_B_bures",  "2q_B_observable"
  )

  for (case_id in case_order) {
    r <- all_results[[case_id]]
    if (is.null(r)) next

    x_start <- min(plot_start, r$n_total)
    x_end <- min(plot_end, r$n_total)
    n_seq <- x_start:x_end

    uniform <- r$results$uniform$mean_risk[x_start:x_end]
    exact <- r$results$exact$mean_risk[x_start:x_end]
    gi1 <- r$results$GI1$mean_risk[x_start:x_end]

    all_r <- c(uniform, exact, gi1)
    y_lim <- compute_plot_ylim(all_r, oracle_limit = r$oracle_limit, floor_zero = TRUE)

    plot(n_seq, uniform, type = "l", col = "gray40", lwd = 1.5,
         xlim = c(x_start, x_end), ylim = y_lim,
         xlab = "n", ylab = expression(n %.% R[n](G)),
         main = sprintf("%s / %s / %s", toupper(r$system), r$library, r$loss),
         cex.main = 1.0)
    lines(n_seq, exact, col = "blue", lwd = 1.5)
    lines(n_seq, gi1, col = "red", lwd = 1.5, lty = 2)
    if (!is.null(r$oracle_limit)) {
      abline(h = r$oracle_limit, col = "black", lwd = 1.2, lty = 3)
    }
    grid(lty = "dotted", col = "gray80")
  }

  mtext("Adaptive Quantum Tomography: 15 Experiment Cases (Oracle Risk)",
        outer = TRUE, cex = 1.5, line = 0.5)
  par(fig = c(0, 1, 0, 0.03), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", horiz = TRUE, legend = c("Uniform", "Exact", "GI1"),
         col = c("gray40", "blue", "red"), lwd = 2, lty = c(1, 1, 2),
         cex = 1.2, bty = "n")
  dev.off()

  cat("Saved combined plot: all_15_cases_combined.png\n")
}

## ==========================================================================
## 7. MAIN ENTRYPOINT
## ==========================================================================

cat("\nRunning full 15-case simulation...\n")
full_results <- run_full_simulation(
  cfg = simulation_config,
  output_dir = "results",
  plot_start = plot_start,
  plot_end = plot_end
)

generate_combined_plot(full_results$results,
                       output_dir = full_results$output_dir,
                       plot_start = plot_start,
                       plot_end = plot_end)

cat("\nSimulation complete. Outputs:\n")
cat(sprintf("  - results/results.csv\n"))
cat(sprintf("  - results/full_results.rds\n"))
cat(sprintf("  - results/simulation_config.rds\n"))
cat(sprintf("  - results/plots/*.png\n"))
cat("\n")
