## ============================================================================
## 11_run_1q_mse_mineig.R
## One-qubit diagnostics: MSE, Frobenius distance, and MLE min eigenvalue
## ============================================================================

cat("==========================================================\n")
cat(" 1Q DIAGNOSTICS: MSE / DISTANCE / MIN EIGENVALUE\n")
cat("==========================================================\n\n")

## ==========================================================================
## 0) SETUP
## ==========================================================================

script_dir <- tryCatch({
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (file.exists(script_path)) {
      return(dirname(normalizePath(script_path)))
    }
  }
  this_file <- sys.frame(1)$ofile
  if (!is.null(this_file) && this_file != "" && file.exists(this_file)) {
    return(dirname(normalizePath(this_file)))
  }
  NULL
}, error = function(e) NULL)

if (!is.null(script_dir) && script_dir != "") {
  setwd(script_dir)
} else if (dir.exists("R_implementations")) {
  setwd("R_implementations")
}

source("01_utilities.R")
source("02_state_basis.R")
source("03_measurement_library_1q.R")
source("04_measurement_library_2q.R")
source("05_mle_cvxr_solver.R")
source("06_fisher_and_metrics.R")
source("07_adaptive_design.R")
source("08_simulation_controller.R")

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
  as.integer(max(1L, nc))
}

compute_plot_ylim <- function(values, ref_line = NULL, floor_zero = TRUE) {
  valid <- values[is.finite(values)]
  if (length(valid) == 0) {
    y_min <- 0
    y_max <- 1
  } else {
    y_range <- range(valid)
    span <- y_range[2] - y_range[1]
    pad <- if (span > 0) 0.05 * span else max(1e-6, abs(y_range[1]) * 0.05 + 1e-6)
    y_min <- y_range[1] - pad
    y_max <- y_range[2] + pad
  }
  if (!is.null(ref_line) && is.finite(ref_line)) {
    ref_pad <- max(1e-6, abs(ref_line) * 0.05)
    y_min <- min(y_min, ref_line - ref_pad)
    y_max <- max(y_max, ref_line + ref_pad)
  }
  if (floor_zero) y_min <- max(0, y_min)
  if (!is.finite(y_min) || !is.finite(y_max) || y_max <= y_min) {
    y_min <- if (is.finite(y_min)) y_min else 0
    y_max <- y_min + max(1e-6, abs(y_min) * 0.05 + 1e-6)
  }
  c(y_min, y_max)
}

summarize_matrix <- function(mat) {
  list(
    mean = colMeans(mat, na.rm = TRUE),
    sd = apply(mat, 2, sd, na.rm = TRUE),
    se = apply(mat, 2, sd, na.rm = TRUE) / sqrt(nrow(mat))
  )
}

get_library <- function(lib_name) {
  switch(lib_name,
    "L1" = build_library_1q_L1(),
    "L2" = build_library_1q_L2(),
    "L3" = build_library_1q_L3(seed = 42),
    stop(paste("Unknown library:", lib_name))
  )
}

get_metric_function <- function(loss_name, sigmas, N, lib) {
  switch(loss_name,
    "frobenius" = make_metric_frobenius(sigmas),
    "bures" = make_metric_bures(sigmas, N),
    "observable" = make_metric_observable(sigmas, lib$Obs),
    stop(paste("Unknown loss:", loss_name))
  )
}

build_true_state_1q <- function(alpha) {
  rho_pure <- fixed_pure_state_1q_yphase(phi = 0.2)
  mix_with_maximally_mixed(rho_pure, alpha)
}

## ==========================================================================
## 1) CONFIG (1Q ONLY)
## ==========================================================================

cfg <- default_sim_config()
cfg$n_rep <- 100
cfg$seed_base <- 2026
cfg$true_state_alpha <- 0.99
cfg$n_total <- 100
cfg$check_every <- 1
cfg$eta_mle <- 1e-3
cfg$ridge_init <- 1e-8
cfg$solver <- "SCS"
cfg$n_workers <- min(9L, detect_parallel_workers(default = 1L))

plot_start <- 20
plot_end <- 100

cat("Configuration:\n")
cat(sprintf("  n_rep = %d\n", cfg$n_rep))
cat(sprintf("  n_total = %d\n", cfg$n_total))
cat(sprintf("  alpha = %.3f\n", cfg$true_state_alpha))
cat(sprintf("  eta_mle = %.2e\n", cfg$eta_mle))
cat(sprintf("  solver = %s\n", cfg$solver))
cat(sprintf("  workers = %d\n", cfg$n_workers))
cat(sprintf("  plot window = %d..%d\n\n", plot_start, plot_end))

experiment_grid_1q <- list(
  list(system = "1q", lib = "L1", loss = "frobenius", case_num = 1),
  list(system = "1q", lib = "L1", loss = "bures", case_num = 2),
  list(system = "1q", lib = "L1", loss = "observable", case_num = 3),
  list(system = "1q", lib = "L2", loss = "frobenius", case_num = 4),
  list(system = "1q", lib = "L2", loss = "bures", case_num = 5),
  list(system = "1q", lib = "L2", loss = "observable", case_num = 6),
  list(system = "1q", lib = "L3", loss = "frobenius", case_num = 7),
  list(system = "1q", lib = "L3", loss = "bures", case_num = 8),
  list(system = "1q", lib = "L3", loss = "observable", case_num = 9)
)

## ==========================================================================
## 2) RUN ONE EXPERIMENT
## ==========================================================================

run_experiment_1q_diagnostics <- function(exp, cfg) {
  case_id <- sprintf("%s_%s_%s", exp$system, exp$lib, exp$loss)
  cat(sprintf("[Case %d/9] %s\n", exp$case_num, case_id))

  lib <- get_library(exp$lib)
  basis <- build_pauli_basis_1q()
  sigmas <- basis$sigmas
  N <- lib$N
  metric_fun <- get_metric_function(exp$loss, sigmas, N, lib)
  rho_true <- build_true_state_1q(cfg$true_state_alpha)

  policy_results <- list()
  for (policy in c("uniform", "exact", "GI1")) {
    cat(sprintf("  - policy %s ... ", policy))
    flush.console()

    risk_mat <- matrix(NA_real_, cfg$n_rep, cfg$n_total)
    mse_mat <- matrix(NA_real_, cfg$n_rep, cfg$n_total)
    dist_mat <- matrix(NA_real_, cfg$n_rep, cfg$n_total)
    mineig_mat <- matrix(NA_real_, cfg$n_rep, cfg$n_total)

    seed_key <- paste0(case_id, "_", policy)
    for (r in seq_len(cfg$n_rep)) {
      rep_seed <- stable_seed(paste0(seed_key, "_rep_", r), base = cfg$seed_base)
      res <- adaptive_design_sequence_metric(
        lib = lib,
        sigmas = sigmas,
        rho_true = rho_true,
        n_total = cfg$n_total,
        metric_fun = metric_fun,
        policy = policy,
        lib_name = exp$lib,
        eta_mle = cfg$eta_mle,
        ridge = cfg$ridge_init,
        solver = cfg$solver,
        seed = rep_seed,
        verbose = FALSE,
        check_every = cfg$check_every,
        record_density_metrics = TRUE
      )
      risk_mat[r, ] <- res$risk_history
      mse_mat[r, ] <- res$mse_history
      dist_mat[r, ] <- res$distance_history
      mineig_mat[r, ] <- res$min_eig_history
    }

    policy_results[[policy]] <- list(
      risk = summarize_matrix(risk_mat),
      mse = summarize_matrix(mse_mat),
      distance = summarize_matrix(dist_mat),
      min_eig = summarize_matrix(mineig_mat)
    )
    cat("done\n")
  }

  list(
    case_num = exp$case_num,
    case_id = case_id,
    system = exp$system,
    library = exp$lib,
    loss = exp$loss,
    n_total = cfg$n_total,
    n_rep = cfg$n_rep,
    eta_mle = cfg$eta_mle,
    results = policy_results
  )
}

## ==========================================================================
## 3) PLOTTING
## ==========================================================================

plot_case_diagnostics <- function(exp_result, output_dir, plot_start = 20, plot_end = 100) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  x_start <- min(plot_start, exp_result$n_total)
  x_end <- min(plot_end, exp_result$n_total)
  n_seq <- x_start:x_end

  pols <- c("uniform", "exact", "GI1")
  cols <- c("gray40", "blue", "red")
  ltys <- c(1, 1, 2)

  pull_metric <- function(metric_name) {
    lapply(pols, function(p) exp_result$results[[p]][[metric_name]]$mean[x_start:x_end])
  }

  mse_vals <- pull_metric("mse")
  dist_vals <- pull_metric("distance")
  mineig_vals <- pull_metric("min_eig")

  lib_full <- switch(exp_result$library,
    "L1" = "L1: Pauli PVMs",
    "L2" = "L2: Nine-axis",
    "L3" = "L3: Random 4-axis"
  )

  f_out <- file.path(output_dir, sprintf("case%02d_%s_%s_%s_diagnostics.png",
                                         exp_result$case_num,
                                         exp_result$system,
                                         exp_result$library,
                                         exp_result$loss))
  png(f_out, width = 1000, height = 1050, res = 110)
  par(mfrow = c(3, 1), mar = c(4, 5, 3, 1), oma = c(1, 1, 2, 1))

  y_mse <- compute_plot_ylim(unlist(mse_vals), floor_zero = TRUE)
  plot(n_seq, mse_vals[[1]], type = "l", col = cols[1], lwd = 2.2, lty = ltys[1],
       xlim = c(x_start, x_end), ylim = y_mse,
       xlab = "Number of Samples (n)", ylab = "MSE of density matrix")
  lines(n_seq, mse_vals[[2]], col = cols[2], lwd = 2.2, lty = ltys[2])
  lines(n_seq, mse_vals[[3]], col = cols[3], lwd = 2.2, lty = ltys[3])
  grid(lty = "dotted", col = "gray82")
  legend("topright", legend = c("Uniform", "Exact", "GI1"),
         col = cols, lwd = 2.2, lty = ltys, bg = "white", box.lwd = 1)
  title(main = "MSE: mean(|rho_hat - rho_true|^2)")

  y_dist <- compute_plot_ylim(unlist(dist_vals), floor_zero = TRUE)
  plot(n_seq, dist_vals[[1]], type = "l", col = cols[1], lwd = 2.2, lty = ltys[1],
       xlim = c(x_start, x_end), ylim = y_dist,
       xlab = "Number of Samples (n)", ylab = "Frobenius distance")
  lines(n_seq, dist_vals[[2]], col = cols[2], lwd = 2.2, lty = ltys[2])
  lines(n_seq, dist_vals[[3]], col = cols[3], lwd = 2.2, lty = ltys[3])
  grid(lty = "dotted", col = "gray82")
  legend("topright", legend = c("Uniform", "Exact", "GI1"),
         col = cols, lwd = 2.2, lty = ltys, bg = "white", box.lwd = 1)
  title(main = "Distance: Frobenius norm ||rho_hat - rho_true||_F")

  y_mineig <- compute_plot_ylim(unlist(mineig_vals), ref_line = exp_result$eta_mle, floor_zero = FALSE)
  plot(n_seq, mineig_vals[[1]], type = "l", col = cols[1], lwd = 2.2, lty = ltys[1],
       xlim = c(x_start, x_end), ylim = y_mineig,
       xlab = "Number of Samples (n)", ylab = "Minimum eigenvalue of MLE")
  lines(n_seq, mineig_vals[[2]], col = cols[2], lwd = 2.2, lty = ltys[2])
  lines(n_seq, mineig_vals[[3]], col = cols[3], lwd = 2.2, lty = ltys[3])
  abline(h = exp_result$eta_mle, col = "black", lwd = 1.8, lty = 3)
  grid(lty = "dotted", col = "gray82")
  legend("topright", legend = c("Uniform", "Exact", "GI1", "eta_mle"),
         col = c(cols, "black"), lwd = c(2.2, 2.2, 2.2, 1.8), lty = c(ltys, 3),
         bg = "white", box.lwd = 1)
  title(main = "MLE stability: minimum eigenvalue")

  mtext(sprintf("Case %d: One-Qubit / %s / %s   (n_rep=%d)",
                exp_result$case_num, lib_full, exp_result$loss, exp_result$n_rep),
        outer = TRUE, cex = 1.15, line = 0.3, font = 2)
  dev.off()

  f_out
}

plot_combined_metric <- function(all_results, metric_name, ylab, filename,
                                 plot_start = 20, plot_end = 100, eta_line = FALSE) {
  png(filename, width = 2400, height = 2000, res = 150)
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1))

  case_order <- sprintf("1q_%s_%s", rep(c("L1", "L2", "L3"), each = 3),
                        rep(c("frobenius", "bures", "observable"), times = 3))
  for (cid in case_order) {
    r <- all_results[[cid]]
    if (is.null(r)) next

    x_start <- min(plot_start, r$n_total)
    x_end <- min(plot_end, r$n_total)
    n_seq <- x_start:x_end

    y_u <- r$results$uniform[[metric_name]]$mean[x_start:x_end]
    y_e <- r$results$exact[[metric_name]]$mean[x_start:x_end]
    y_g <- r$results$GI1[[metric_name]]$mean[x_start:x_end]
    y_lim <- compute_plot_ylim(c(y_u, y_e, y_g),
                               ref_line = if (eta_line) r$eta_mle else NULL,
                               floor_zero = !eta_line)

    plot(n_seq, y_u, type = "l", col = "gray40", lwd = 1.8,
         xlim = c(x_start, x_end), ylim = y_lim,
         xlab = "n", ylab = ylab,
         main = sprintf("%s / %s", r$library, r$loss), cex.main = 1.0)
    lines(n_seq, y_e, col = "blue", lwd = 1.8)
    lines(n_seq, y_g, col = "red", lwd = 1.8, lty = 2)
    if (eta_line) abline(h = r$eta_mle, col = "black", lwd = 1.4, lty = 3)
    grid(lty = "dotted", col = "gray82")
  }

  mtext(sprintf("One-Qubit Diagnostics (%s) | n=%d..%d", metric_name, plot_start, plot_end),
        outer = TRUE, cex = 1.4, line = 0.5)
  par(fig = c(0, 1, 0, 0.03), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  if (eta_line) {
    legend("center", horiz = TRUE, legend = c("Uniform", "Exact", "GI1", "eta_mle"),
           col = c("gray40", "blue", "red", "black"), lwd = c(2, 2, 2, 1.6),
           lty = c(1, 1, 2, 3), cex = 1.15, bty = "n")
  } else {
    legend("center", horiz = TRUE, legend = c("Uniform", "Exact", "GI1"),
           col = c("gray40", "blue", "red"), lwd = 2, lty = c(1, 1, 2),
           cex = 1.15, bty = "n")
  }
  dev.off()
}

## ==========================================================================
## 4) RUN ALL 1Q CASES
## ==========================================================================

run_all_1q <- function(cfg, output_dir = "results_1q_diagnostics",
                       plot_start = 20, plot_end = 100) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  n_workers <- min(cfg$n_workers, length(experiment_grid_1q))
  if (.Platform$OS.type == "windows" && n_workers > 1L) n_workers <- 1L

  if (n_workers > 1L) {
    cat(sprintf("Running 1q cases in parallel with %d workers...\n", n_workers))
    exp_results <- parallel::mclapply(
      experiment_grid_1q,
      function(exp) run_experiment_1q_diagnostics(exp, cfg),
      mc.cores = n_workers,
      mc.preschedule = FALSE
    )
  } else {
    cat("Running 1q cases sequentially...\n")
    exp_results <- lapply(experiment_grid_1q, function(exp) run_experiment_1q_diagnostics(exp, cfg))
  }

  all_results <- list()
  all_df <- data.frame()
  for (res in exp_results) {
    all_results[[res$case_id]] <- res
    for (policy in c("uniform", "exact", "GI1")) {
      rr <- res$results[[policy]]
      dfp <- data.frame(
        case_num = res$case_num,
        case_id = res$case_id,
        system = res$system,
        library = res$library,
        loss = res$loss,
        policy = policy,
        n = seq_len(res$n_total),
        mse_mean = rr$mse$mean,
        mse_sd = rr$mse$sd,
        mse_se = rr$mse$se,
        distance_mean = rr$distance$mean,
        distance_sd = rr$distance$sd,
        distance_se = rr$distance$se,
        min_eig_mean = rr$min_eig$mean,
        min_eig_sd = rr$min_eig$sd,
        min_eig_se = rr$min_eig$se
      )
      all_df <- rbind(all_df, dfp)
    }
    plot_case_diagnostics(res, plots_dir, plot_start = plot_start, plot_end = plot_end)
  }

  write.csv(all_df, file.path(output_dir, "diagnostics_1q_curves.csv"), row.names = FALSE)
  saveRDS(all_results, file.path(output_dir, "diagnostics_1q_full_results.rds"))
  saveRDS(cfg, file.path(output_dir, "diagnostics_1q_config.rds"))

  plot_combined_metric(
    all_results = all_results,
    metric_name = "mse",
    ylab = "MSE",
    filename = file.path(plots_dir, "combined_1q_mse.png"),
    plot_start = plot_start,
    plot_end = plot_end,
    eta_line = FALSE
  )
  plot_combined_metric(
    all_results = all_results,
    metric_name = "distance",
    ylab = "Frobenius distance",
    filename = file.path(plots_dir, "combined_1q_distance.png"),
    plot_start = plot_start,
    plot_end = plot_end,
    eta_line = FALSE
  )
  plot_combined_metric(
    all_results = all_results,
    metric_name = "min_eig",
    ylab = "Min eigenvalue",
    filename = file.path(plots_dir, "combined_1q_min_eig.png"),
    plot_start = plot_start,
    plot_end = plot_end,
    eta_line = TRUE
  )

  cat("\nOutputs written to:\n")
  cat(sprintf("  - %s\n", file.path(output_dir, "diagnostics_1q_curves.csv")))
  cat(sprintf("  - %s\n", file.path(output_dir, "diagnostics_1q_full_results.rds")))
  cat(sprintf("  - %s\n", file.path(output_dir, "diagnostics_1q_config.rds")))
  cat(sprintf("  - %s\n", file.path(plots_dir, "*.png")))
  invisible(list(all_results = all_results, curves = all_df, cfg = cfg, output_dir = output_dir))
}

## ==========================================================================
## 5) MAIN
## ==========================================================================

run_all_1q(
  cfg = cfg,
  output_dir = "results_1q_diagnostics_nrep100",
  plot_start = plot_start,
  plot_end = plot_end
)
