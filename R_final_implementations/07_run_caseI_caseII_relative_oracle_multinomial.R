## ============================================================================
## 07_run_caseI_caseII_relative_oracle_multinomial.R
## Case I/II (alpha = 0.5, 0.9) Monte Carlo study with multinomial batches.
## - 15-case grid (1q L1/L2/L3 and 2q A/B x frobenius/bures/observable)
## - Policies: permutation-baseline(key='uniform') / exact / GI1 / oracle_GI1
## - Adaptive selection count n_total
## - Multinomial trials per adaptive selection n_trials
## - MLE refreshed every adaptive selection (check_every = 1)
## ============================================================================

cat("==============================================================\n")
cat(" CASE I/II MULTINOMIAL ADAPTIVE TOMOGRAPHY (15 CASES)\n")
cat("==============================================================\n\n")

script_dir <- tryCatch({
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  this_file <- sys.frame(1)$ofile
  if (!is.null(this_file) && this_file != "") {
    return(dirname(normalizePath(this_file)))
  }
  NULL
}, error = function(e) NULL)

if (!is.null(script_dir) && script_dir != "") {
  if (!file.exists(file.path(script_dir, "00_simulation_parameters.R"))) {
    candidate <- file.path(script_dir, "R_final_implementations")
    if (file.exists(file.path(candidate, "00_simulation_parameters.R"))) script_dir <- candidate
  }
  setwd(script_dir)
} else if (dir.exists("R_final_implementations")) {
  setwd("R_final_implementations")
} else if (dir.exists("R_implementations_Feb")) {
  setwd("R_implementations_Feb")
}

source("01_utilities.R")
source("02_state_basis.R")
source("03_measurement_library_1q.R")
source("04_measurement_library_2q.R")
source("05_mle_cvxr_solver.R")
source("06_fisher_and_metrics.R")
source(file.path(getwd(), "00_simulation_parameters.R"))

## ============================================================================
## Helpers
## ============================================================================

detect_parallel_workers <- function(default = 1L) {
  nc <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(nc) || nc < 1L) nc <- as.integer(default)
  as.integer(max(1L, nc))
}

parse_cli_value <- function(prefix, args = commandArgs(trailingOnly = TRUE), default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[1])
}

as_int_or <- function(x, default) {
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) default else v
}

as_num_or <- function(x, default) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) default else v
}

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_alpha_vec <- function(x, default = c(0.5, 0.9)) {
  vals <- suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(default)
  unique(vals)
}

parse_int_vec <- function(x, default = integer()) {
  vals <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals)]
  vals <- vals[!is.na(vals)]
  vals <- unique(vals)
  if (length(vals) == 0) return(default)
  vals
}

alpha_case_name <- function(alpha) {
  if (abs(alpha - 0.5) < 1e-12) return("Case I")
  if (abs(alpha - 0.9) < 1e-12) return("Case II")
  sprintf("alpha=%.3f", alpha)
}

alpha_tag <- function(alpha) {
  gsub("\\.", "p", format(alpha, nsmall = 1, trim = TRUE))
}

safe_inverse <- function(M, ridge = 0) {
  d <- nrow(M)
  M2 <- (M + t(M)) / 2 + ridge * diag(d)
  solve(M2)
}

matrix_sqrt_psd <- function(M, tol = 0) {
  eg <- eigen(hermitianize(M), symmetric = FALSE)
  vals <- Re(eg$values)
  if (any(!is.finite(vals))) stop("matrix_sqrt_psd: non-finite eigenvalues.")
  if (tol > 0 && any(vals < -tol)) stop("matrix_sqrt_psd: matrix is not PSD.")
  vals <- pmax(vals, 0)
  U <- eg$vectors
  out <- U %*% diag(sqrt(vals), nrow = length(vals)) %*% Conj(t(U))
  hermitianize(out)
}

bures_distance_density <- function(rho, sigma, tol = 0) {
  sqrt_rho <- matrix_sqrt_psd(rho, tol = tol)
  inner <- hermitianize(sqrt_rho %*% sigma %*% sqrt_rho)
  vals <- Re(eigen(inner, symmetric = FALSE, only.values = TRUE)$values)
  if (any(!is.finite(vals))) stop("bures_distance_density: non-finite eigenvalues.")
  if (tol > 0 && any(vals < -tol)) stop("bures_distance_density: negative eigenvalue in fidelity matrix.")
  vals <- pmax(vals, 0)
  root_fid <- sum(sqrt(vals))
  if (!is.finite(root_fid)) stop("bures_distance_density: non-finite fidelity.")
  dist_sq <- 2 - 2 * root_fid
  if (tol > 0 && dist_sq < -tol) stop("bures_distance_density: negative squared distance.")
  sqrt(max(0, dist_sq))
}

increment_count_batch <- function(Nab, a, batch_counts, ab_row) {
  idx <- ab_row[[a]]
  Nab[idx] <- Nab[idx] + as.numeric(batch_counts)
  Nab
}

get_init_settings <- function(lib_name, k) {
  switch(lib_name,
    "L1" = 1:3,
    "L2" = 1:3,
    "L3" = 1:4,
    "A" = 1:15,
    "B" = 1:5,
    1:min(k, 3)
  )
}

x_idx_for_system <- function(system, n_total, cfg) {
  x_start <- if (identical(system, "1q")) as.integer(cfg$plot_1q_start) else as.integer(cfg$plot_2q_start)
  x_end <- as.integer(n_total)
  if (x_start > x_end) x_start <- 1L
  seq.int(x_start, x_end)
}

get_case_grid <- function() {
  list(
    list(case_num = 1,  case_id = "1q_L1_frobenius",  system = "1q", library = "L1", loss = "frobenius",  library_full = "L1: Pauli PVMs"),
    list(case_num = 2,  case_id = "1q_L1_bures",      system = "1q", library = "L1", loss = "bures",      library_full = "L1: Pauli PVMs"),
    list(case_num = 3,  case_id = "1q_L1_observable", system = "1q", library = "L1", loss = "observable", library_full = "L1: Pauli PVMs"),
    list(case_num = 4,  case_id = "1q_L2_frobenius",  system = "1q", library = "L2", loss = "frobenius",  library_full = "L2: Nine-axis PVMs"),
    list(case_num = 5,  case_id = "1q_L2_bures",      system = "1q", library = "L2", loss = "bures",      library_full = "L2: Nine-axis PVMs"),
    list(case_num = 6,  case_id = "1q_L2_observable", system = "1q", library = "L2", loss = "observable", library_full = "L2: Nine-axis PVMs"),
    list(case_num = 7,  case_id = "1q_L3_frobenius",  system = "1q", library = "L3", loss = "frobenius",  library_full = "L3: Random 4-axis PVMs"),
    list(case_num = 8,  case_id = "1q_L3_bures",      system = "1q", library = "L3", loss = "bures",      library_full = "L3: Random 4-axis PVMs"),
    list(case_num = 9,  case_id = "1q_L3_observable", system = "1q", library = "L3", loss = "observable", library_full = "L3: Random 4-axis PVMs"),
    list(case_num = 10, case_id = "2q_A_frobenius",   system = "2q", library = "A",  loss = "frobenius",  library_full = "A: 15 Pauli-parity PVMs"),
    list(case_num = 11, case_id = "2q_A_bures",       system = "2q", library = "A",  loss = "bures",      library_full = "A: 15 Pauli-parity PVMs"),
    list(case_num = 12, case_id = "2q_A_observable",  system = "2q", library = "A",  loss = "observable", library_full = "A: 15 Pauli-parity PVMs"),
    list(case_num = 13, case_id = "2q_B_frobenius",   system = "2q", library = "B",  loss = "frobenius",  library_full = "B: 5 MUB bases"),
    list(case_num = 14, case_id = "2q_B_bures",       system = "2q", library = "B",  loss = "bures",      library_full = "B: 5 MUB bases"),
    list(case_num = 15, case_id = "2q_B_observable",  system = "2q", library = "B",  loss = "observable", library_full = "B: 5 MUB bases")
  )
}

get_library <- function(lib_name, l3_seed = 42) {
  switch(lib_name,
    "L1" = build_library_1q_L1(),
    "L2" = build_library_1q_L2(),
    "L3" = build_library_1q_L3(seed = l3_seed),
    "A" = build_library_2q_A(),
    "B" = build_library_2q_B(),
    stop(sprintf("Unknown library: %s", lib_name))
  )
}

get_basis <- function(system) {
  if (system == "1q") return(build_pauli_basis_1q())
  if (system == "2q") return(build_pauli_product_basis_2q())
  stop(sprintf("Unknown system: %s", system))
}

build_true_state <- function(system, alpha) {
  if (system == "1q") {
    rho_pure <- fixed_pure_state_1q_yphase(phi = 0.2)
    return(mix_with_maximally_mixed(rho_pure, alpha))
  }
  if (system == "2q") {
    rho_pure <- fixed_pure_state_2q_phi_phase(phi = 0.2)
    return(mix_with_maximally_mixed(rho_pure, alpha))
  }
  stop(sprintf("Unknown system: %s", system))
}

get_metric_function <- function(loss_name, sigmas, N, lib, bures_tol = 0) {
  switch(loss_name,
    "frobenius" = make_metric_frobenius(sigmas),
    "bures" = function(theta) metric_bures(theta, sigmas, N, tol = bures_tol),
    "observable" = {
      if (is.null(lib$Obs)) stop("Library missing Obs for observable metric.")
      make_metric_observable(sigmas, lib$Obs)
    },
    stop(sprintf("Unknown loss: %s", loss_name))
  )
}

select_next_setting_metric_batch <- function(theta_hat, counts, S_ab, c_ab, ab_df,
                                             G, method = c("exact", "GI1"),
                                             n_trials = 100L,
                                             ridge = 0,
                                             eps = 0) {
  method <- match.arg(method)

  tf <- total_fisher_info(theta_hat, counts, S_ab, c_ab, ab_df, eps = eps)
  J_n <- tf$J_total
  I_list <- tf$I_list
  d <- ncol(S_ab)

  J_reg <- (J_n + t(J_n)) / 2 + ridge * diag(d)
  scores <- numeric(length(I_list))
  names(scores) <- names(I_list)

  if (method == "exact") {
    for (i in seq_along(I_list)) {
      Ia <- I_list[[i]]
      J_cand <- J_reg + as.numeric(n_trials) * Ia
      scores[i] <- tryCatch(
        sum(diag(G %*% safe_inverse(J_cand, ridge = 0))),
        error = function(e) Inf
      )
    }
    a_next <- as.integer(names(scores)[which.min(scores)])
  } else {
    J_inv <- tryCatch(safe_inverse(J_reg, ridge = 0), error = function(e) NULL)
    if (is.null(J_inv)) {
      stop("GI1 selection failed: singular accumulated Fisher matrix and ridge=0.")
    }
    for (i in seq_along(I_list)) {
      Ia <- I_list[[i]]
      scores[i] <- as.numeric(n_trials) * sum(diag(G %*% J_inv %*% Ia %*% J_inv))
    }
    a_next <- as.integer(names(scores)[which.max(scores)])
  }

  list(a_next = a_next, scores = scores)
}

solve_optimal_design <- function(I_list, G, solver = "SCS") {
  if (!isTRUE(get("CVXR_AVAILABLE", inherits = TRUE))) {
    return(list(status = "cvxr_unavailable", pi = NULL))
  }

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
  res <- tryCatch(
    CVXR::solve(prob, solver = solver, verbose = FALSE, max_iters = 50000),
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(list(status = "cvxr_error", pi = NULL))
  }

  list(status = res$status, pi = as.numeric(res$getValue(pi)))
}

solve_optimal_design_mirror <- function(I_list, G, ridge = 0, max_iter = 5000,
                                        tol = 1e-10, eta = 1.0) {
  k <- length(I_list)
  d <- nrow(G)
  pi <- rep(1 / k, k)

  f_val <- Inf
  for (iter in seq_len(max_iter)) {
    I_pi <- Reduce("+", lapply(seq_len(k), function(i) pi[i] * I_list[[i]]))
    I_pi <- (I_pi + t(I_pi)) / 2 + ridge * diag(d)
    J_inv <- tryCatch(safe_inverse(I_pi, ridge = 0), error = function(e) NULL)
    if (is.null(J_inv)) break
    f_val <- sum(diag(G %*% J_inv))

    grad <- numeric(k)
    for (i in seq_len(k)) {
      grad[i] <- -sum(diag(G %*% J_inv %*% I_list[[i]] %*% J_inv))
    }

    eta_i <- eta
    improved <- FALSE
    pi_new <- pi
    f_new <- f_val

    for (bt in 1:25) {
      cand <- pi * exp(-eta_i * grad)
      cand <- cand / sum(cand)
      I_c <- Reduce("+", lapply(seq_len(k), function(i) cand[i] * I_list[[i]]))
      I_c <- (I_c + t(I_c)) / 2 + ridge * diag(d)
      f_c <- tryCatch(sum(diag(G %*% safe_inverse(I_c, ridge = 0))), error = function(e) Inf)
      if (is.finite(f_c) && f_c <= f_val) {
        pi_new <- cand
        f_new <- f_c
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

  list(pi = pi, value = f_val)
}

build_case_object <- function(exp, alpha, cfg) {
  lib <- get_library(exp$library, l3_seed = cfg$l3_seed)
  basis <- get_basis(exp$system)
  sigmas <- basis$sigmas
  N <- lib$N

  rho_true <- build_true_state(exp$system, alpha)
  theta_true <- theta_from_rho(rho_true, sigmas)

  metric_fun <- get_metric_function(exp$loss, sigmas, N, lib, bures_tol = cfg$bures_tol)
  G_true <- metric_fun(theta_true)

  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)
  I_list <- fisher_info_all_settings(theta_true, sc$S_ab, sc$c_ab, lib$ab_df)

  k <- length(I_list)
  I_uniform <- Reduce("+", I_list) / k
  uniform_limit <- proxy_risk(G_true, I_uniform, ridge = cfg$ridge_selection)

  opt <- solve_optimal_design(I_list, G_true, solver = cfg$oracle_solver)
  oracle_status <- opt$status
  pi_opt <- opt$pi

  use_mirror <- is.null(pi_opt) || any(!is.finite(pi_opt)) ||
    !(oracle_status %in% c("optimal", "optimal_inaccurate"))

  if (!use_mirror) {
    I_opt <- Reduce("+", lapply(seq_len(k), function(i) pi_opt[i] * I_list[[i]]))
    oracle_limit <- proxy_risk(G_true, I_opt, ridge = cfg$ridge_selection)
    if (!is.finite(oracle_limit) || oracle_limit <= 0 || oracle_limit > uniform_limit * (1 + 1e-8)) {
      use_mirror <- TRUE
    }
  }

  if (use_mirror) {
    md <- solve_optimal_design_mirror(I_list, G_true, ridge = cfg$ridge_selection)
    pi_opt <- md$pi
    I_opt <- Reduce("+", lapply(seq_len(k), function(i) pi_opt[i] * I_list[[i]]))
    oracle_limit <- proxy_risk(G_true, I_opt, ridge = cfg$ridge_selection)
    oracle_status <- paste0(oracle_status, " -> mirror")
  }

  init_settings <- get_init_settings(exp$library, lib$k)
  if (length(init_settings) > cfg$n_total) {
    stop(sprintf("n_total=%d is too small for initialization length %d in case %s",
                 cfg$n_total, length(init_settings), exp$case_id))
  }

  list(
    case_num = exp$case_num,
    case_id = exp$case_id,
    system = exp$system,
    library = exp$library,
    loss = exp$loss,
    library_full = exp$library_full,
    alpha = alpha,
    N = N,
    lib = lib,
    sigmas = sigmas,
    rho_true = rho_true,
    theta_true = theta_true,
    metric_fun = metric_fun,
    G_true = G_true,
    S_ab = sc$S_ab,
    c_ab = sc$c_ab,
    oracle_limit = oracle_limit,
    uniform_limit = uniform_limit,
    ratio_uniform_over_oracle = uniform_limit / oracle_limit,
    oracle_status = oracle_status,
    init_settings = init_settings,
    n_init = length(init_settings),
    prob_list_true = born_probs_list(rho_true, lib$Q_list),
    oracle_mse_limit = if (exp$loss == "frobenius") oracle_limit / (N * N) else NA_real_
  )
}

mle_fit_success <- function(fit) {
  fit$status %in% c("optimal", "optimal_inaccurate", "trivial")
}

solver_map_val <- function(fit, solver_name, field, fallback = NA_real_) {
  sm <- fit$solver_status_map
  if (is.null(sm) || !is.data.frame(sm) || !(field %in% names(sm))) return(fallback)
  idx <- which(sm$solver == solver_name)
  if (length(idx) == 0L) return(fallback)
  v <- suppressWarnings(as.numeric(sm[[field]][idx[1]]))
  if (length(v) == 0L || is.na(v)) fallback else v
}

solver_map_status <- function(fit, solver_name, fallback = NA_character_) {
  sm <- fit$solver_status_map
  if (is.null(sm) || !is.data.frame(sm) || !("status" %in% names(sm))) return(fallback)
  idx <- which(sm$solver == solver_name)
  if (length(idx) == 0L) return(fallback)
  as.character(sm$status[idx[1]])
}

extract_fit_diagnostics <- function(fit, cfg, cvx_runtime_sec = NA_real_, pgd_runtime_sec = NA_real_) {
  cvx_backend <- cfg$cvx_solver_preference
  if (cvx_backend %in% c("PGD", "AUTO", "AUTO_PGD")) cvx_backend <- "SCS"

  cvx_status <- solver_map_status(fit, cvx_backend, fallback = if (identical(fit$solver_used, cvx_backend)) fit$status else NA_character_)
  pgd_status <- solver_map_status(fit, "PGD", fallback = if (identical(fit$solver_used, "PGD")) fit$status else NA_character_)

  cvx_loss <- solver_map_val(fit, cvx_backend, "train_loss", fallback = NA_real_)
  pgd_loss <- solver_map_val(fit, "PGD", "train_loss", fallback = NA_real_)

  if (!is.finite(cvx_loss) && identical(fit$solver_used, cvx_backend)) cvx_loss <- as.numeric(fit$value)
  if (!is.finite(pgd_loss) && identical(fit$solver_used, "PGD")) pgd_loss <- as.numeric(fit$value)

  selected_loss <- if (identical(fit$solver_used, "PGD")) pgd_loss else cvx_loss
  if (!is.finite(selected_loss)) selected_loss <- as.numeric(fit$value)

  list(
    solver_used = as.character(fit$solver_used %||% NA_character_),
    cvx_status = as.character(cvx_status),
    pgd_status = as.character(pgd_status),
    cvx_train_loss = as.numeric(cvx_loss),
    pgd_train_loss = as.numeric(pgd_loss),
    selected_train_loss = as.numeric(selected_loss),
    cvx_runtime_sec = as.numeric(cvx_runtime_sec),
    pgd_runtime_sec = as.numeric(pgd_runtime_sec)
  )
}

fit_mle_step <- function(case_obj, Nab, theta_hat, cfg) {
  N <- case_obj$N
  sigmas <- case_obj$sigmas

  mode <- cfg$solver_selection_mode_resolved
  cvx_backend <- cfg$cvx_solver_preference
  if (cvx_backend %in% c("PGD", "AUTO", "AUTO_PGD")) cvx_backend <- "SCS"

  if (identical(mode, "always_dual")) {
    fit <- fit_theta_dual_compare(
      N = N,
      sigmas = sigmas,
      S_ab = case_obj$S_ab,
      c_ab = case_obj$c_ab,
      Nab = Nab,
      eta = cfg$eta_mle,
      cvx_solver = cvx_backend,
      eps_log = cfg$nll_eps_log,
      verbose = FALSE,
      theta_start = theta_hat,
      truncate_density = TRUE,
      allow_eta_zero = (cfg$eta_mle <= 0)
    )
    diag <- list(
      solver_used = as.character(fit$solver_used %||% fit$solver_dual_selected %||% NA_character_),
      cvx_status = as.character(fit$cvx_status %||% NA_character_),
      pgd_status = as.character(fit$pgd_status %||% NA_character_),
      cvx_train_loss = as.numeric(fit$cvx_train_loss %||% NA_real_),
      pgd_train_loss = as.numeric(fit$pgd_train_loss %||% NA_real_),
      selected_train_loss = as.numeric(fit$selected_train_loss %||% fit$value %||% NA_real_),
      cvx_runtime_sec = as.numeric(fit$cvx_runtime_sec %||% NA_real_),
      pgd_runtime_sec = as.numeric(fit$pgd_runtime_sec %||% NA_real_)
    )
    return(list(fit = fit, diag = diag))
  }

  if (identical(mode, "single_pgd") || identical(cfg$production_solver, "PGD")) {
    fit <- fit_theta_cvxr(
      N = N,
      sigmas = sigmas,
      S_ab = case_obj$S_ab,
      c_ab = case_obj$c_ab,
      Nab = Nab,
      eta = cfg$eta_mle,
      eps_log = cfg$nll_eps_log,
      solver = "PGD",
      verbose = FALSE,
      theta_start = theta_hat,
      truncate_density = TRUE,
      allow_eta_zero = (cfg$eta_mle <= 0)
    )
    return(list(fit = fit, diag = extract_fit_diagnostics(fit, cfg)))
  }

  t_cvx0 <- proc.time()[["elapsed"]]
  fit_cvx <- tryCatch(
    fit_theta_cvxr(
      N = N,
      sigmas = sigmas,
      S_ab = case_obj$S_ab,
      c_ab = case_obj$c_ab,
      Nab = Nab,
      eta = cfg$eta_mle,
      eps_log = cfg$nll_eps_log,
      solver = cvx_backend,
      verbose = FALSE,
      theta_start = theta_hat,
      truncate_density = TRUE,
      allow_eta_zero = (cfg$eta_mle <= 0)
    ),
    error = function(e) {
      list(
        theta_hat = theta_hat,
        rho_hat = rho_of_theta(theta_hat, sigmas, N),
        status = "error",
        value = NA_real_,
        solver_used = NA_character_,
        solver_status_map = NULL
      )
    }
  )
  cvx_elapsed <- proc.time()[["elapsed"]] - t_cvx0

  if (mle_fit_success(fit_cvx)) {
    return(list(fit = fit_cvx, diag = extract_fit_diagnostics(fit_cvx, cfg, cvx_runtime_sec = cvx_elapsed, pgd_runtime_sec = NA_real_)))
  }

  t_pgd0 <- proc.time()[["elapsed"]]
  fit_pgd <- fit_theta_cvxr(
    N = N,
    sigmas = sigmas,
    S_ab = case_obj$S_ab,
    c_ab = case_obj$c_ab,
    Nab = Nab,
    eta = cfg$eta_mle,
    eps_log = cfg$nll_eps_log,
    solver = "PGD",
    verbose = FALSE,
    theta_start = theta_hat,
    truncate_density = TRUE,
    allow_eta_zero = (cfg$eta_mle <= 0)
  )
  pgd_elapsed <- proc.time()[["elapsed"]] - t_pgd0

  diag <- extract_fit_diagnostics(fit_pgd, cfg, cvx_runtime_sec = cvx_elapsed, pgd_runtime_sec = pgd_elapsed)
  if (is.na(diag$cvx_status) || diag$cvx_status == "") diag$cvx_status <- as.character(fit_cvx$status %||% "error")
  diag$cvx_train_loss <- suppressWarnings(as.numeric(fit_cvx$value %||% NA_real_))
  list(fit = fit_pgd, diag = diag)
}

run_one_replicate_batch <- function(case_obj, policy, cfg, rep_seed) {
  if (!is.null(rep_seed)) set.seed(rep_seed)

  lib <- case_obj$lib
  N <- case_obj$N
  sigmas <- case_obj$sigmas
  d <- length(sigmas)

  Nab <- numeric(lib$M)
  theta_hat <- rep(0, d)
  rho_hat_cur <- diag(N) / N

  n_total <- as.integer(cfg$n_total)
  n_trials <- as.integer(cfg$n_trials)

  scaled_proxy <- numeric(n_total)
  scaled_mse <- numeric(n_total)
  scaled_fro_sq <- numeric(n_total)
  scaled_bures_sq <- numeric(n_total)

  keep_traj <- isTRUE(cfg$save_full_trajectory)
  if (keep_traj) {
    action_idx <- integer(n_total)
    action_label <- character(n_total)
    theta_path <- matrix(NA_real_, nrow = n_total, ncol = d)
    solver_used <- character(n_total)
    cvx_status <- character(n_total)
    pgd_status <- character(n_total)
    cvx_train_loss <- rep(NA_real_, n_total)
    pgd_train_loss <- rep(NA_real_, n_total)
    selected_train_loss <- rep(NA_real_, n_total)
    cvx_runtime_sec <- rep(NA_real_, n_total)
    pgd_runtime_sec <- rep(NA_real_, n_total)
  }

  for (t in seq_len(n_total)) {
    if (t <= case_obj$n_init) {
      a_t <- case_obj$init_settings[t]
    } else if (policy == "uniform") {
      if (identical(cfg$baseline_policy_mode, "permutation_cyclic")) {
        # Deterministic permutation baseline over settings after initialization.
        a_t <- ((t - case_obj$n_init - 1L) %% lib$k) + 1L
      } else {
        a_t <- sample.int(lib$k, size = 1)
      }
    } else {
      counts <- counts_by_setting(Nab, lib$ab_df)
      if (policy == "oracle_GI1") {
        sel <- select_next_setting_metric_batch(
          theta_hat = case_obj$theta_true,
          counts = counts,
          S_ab = case_obj$S_ab,
          c_ab = case_obj$c_ab,
          ab_df = lib$ab_df,
          G = case_obj$G_true,
          method = "GI1",
          n_trials = n_trials,
          ridge = cfg$ridge_selection,
          eps = cfg$fisher_eps
        )
      } else {
        G_hat <- case_obj$metric_fun(theta_hat)
        sel <- select_next_setting_metric_batch(
          theta_hat = theta_hat,
          counts = counts,
          S_ab = case_obj$S_ab,
          c_ab = case_obj$c_ab,
          ab_df = lib$ab_df,
          G = G_hat,
          method = policy,
          n_trials = n_trials,
          ridge = cfg$ridge_selection,
          eps = cfg$fisher_eps
        )
      }
      a_t <- sel$a_next
    }

    probs_a <- case_obj$prob_list_true[[a_t]]
    batch_counts <- as.integer(rmultinom(1, size = n_trials, prob = probs_a)[, 1])
    Nab <- increment_count_batch(Nab, a_t, batch_counts, lib$ab_row)

    fit_out <- fit_mle_step(case_obj, Nab, theta_hat, cfg)
    fit <- fit_out$fit
    diag_fit <- fit_out$diag

    if (mle_fit_success(fit)) {
      theta_hat <- fit$theta_hat
      if (!is.null(fit$rho_hat) && is.matrix(fit$rho_hat)) {
        rho_hat_cur <- hermitianize(fit$rho_hat)
      } else {
        rho_hat_cur <- hermitianize(rho_of_theta(theta_hat, sigmas, N))
      }
    }

    rho_hat_cur <- enforce_rho_floor(rho_hat_cur, eta = 0, sigmas = sigmas)$rho_hat
    diff <- rho_hat_cur - case_obj$rho_true
    mse_val <- mean(Mod(diff)^2)
    fro_sq <- sum(Mod(diff)^2)
    bures_sq <- bures_distance_density(rho_hat_cur, case_obj$rho_true, tol = cfg$bures_tol)^2

    counts_now <- counts_by_setting(Nab, lib$ab_df)
    tf_true <- total_fisher_info(case_obj$theta_true, counts_now,
                                 case_obj$S_ab, case_obj$c_ab, lib$ab_df,
                                 eps = cfg$fisher_eps)
    proxy_raw <- proxy_risk(case_obj$G_true, tf_true$J_total, ridge = cfg$ridge_selection)

    n_eff <- t * n_trials
    scaled_proxy[t] <- n_eff * proxy_raw
    scaled_mse[t] <- n_eff * mse_val
    scaled_fro_sq[t] <- n_eff * fro_sq
    scaled_bures_sq[t] <- n_eff * bures_sq

    if (keep_traj) {
      action_idx[t] <- a_t
      action_label[t] <- case_obj$lib$setting_labels[a_t]
      theta_path[t, ] <- theta_hat
      solver_used[t] <- as.character(diag_fit$solver_used %||% NA_character_)
      cvx_status[t] <- as.character(diag_fit$cvx_status %||% NA_character_)
      pgd_status[t] <- as.character(diag_fit$pgd_status %||% NA_character_)
      cvx_train_loss[t] <- as.numeric(diag_fit$cvx_train_loss %||% NA_real_)
      pgd_train_loss[t] <- as.numeric(diag_fit$pgd_train_loss %||% NA_real_)
      selected_train_loss[t] <- as.numeric(diag_fit$selected_train_loss %||% NA_real_)
      cvx_runtime_sec[t] <- as.numeric(diag_fit$cvx_runtime_sec %||% NA_real_)
      pgd_runtime_sec[t] <- as.numeric(diag_fit$pgd_runtime_sec %||% NA_real_)
    }
  }

  out <- list(
    scaled_proxy = scaled_proxy,
    scaled_mse = scaled_mse,
    scaled_fro_sq = scaled_fro_sq,
    scaled_bures_sq = scaled_bures_sq,
    rep_seed = rep_seed
  )

  if (keep_traj) {
    out$trajectory <- list(
      action_idx = action_idx,
      action_label = action_label,
      theta_hat = theta_path,
      scaled_proxy = scaled_proxy,
      scaled_mse = scaled_mse,
      scaled_fro_sq = scaled_fro_sq,
      scaled_bures_sq = scaled_bures_sq,
      oracle_limit_step = rep(case_obj$oracle_limit, n_total),
      oracle_gi1_proxy_step = if (policy == "oracle_GI1") scaled_proxy else rep(NA_real_, n_total),
      solver_used = solver_used,
      cvx_status = cvx_status,
      pgd_status = pgd_status,
      cvx_train_loss = cvx_train_loss,
      pgd_train_loss = pgd_train_loss,
      selected_train_loss = selected_train_loss,
      cvx_runtime_sec = cvx_runtime_sec,
      pgd_runtime_sec = pgd_runtime_sec
    )
  }

  out
}

summarize_matrix <- function(mat) {
  list(
    mean = colMeans(mat, na.rm = TRUE),
    sd = apply(mat, 2, sd, na.rm = TRUE),
    se = apply(mat, 2, sd, na.rm = TRUE) / sqrt(nrow(mat))
  )
}

make_trajectory_chunk <- function(rep_ids, reps, n_total, d) {
  list(
    replicate_id = as.integer(rep_ids),
    rep_seed = vapply(reps, function(x) as.integer(x$rep_seed), integer(1)),
    action_idx = do.call(rbind, lapply(reps, function(x) x$trajectory$action_idx)),
    action_label = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$action_label, nrow = 1))),
    theta_hat = lapply(reps, function(x) x$trajectory$theta_hat),
    scaled_proxy = do.call(rbind, lapply(reps, function(x) x$trajectory$scaled_proxy)),
    scaled_mse = do.call(rbind, lapply(reps, function(x) x$trajectory$scaled_mse)),
    scaled_fro_sq = do.call(rbind, lapply(reps, function(x) x$trajectory$scaled_fro_sq)),
    scaled_bures_sq = do.call(rbind, lapply(reps, function(x) x$trajectory$scaled_bures_sq)),
    oracle_limit_step = do.call(rbind, lapply(reps, function(x) x$trajectory$oracle_limit_step)),
    oracle_gi1_proxy_step = do.call(rbind, lapply(reps, function(x) x$trajectory$oracle_gi1_proxy_step)),
    solver_used = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$solver_used, nrow = 1))),
    cvx_status = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$cvx_status, nrow = 1))),
    pgd_status = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$pgd_status, nrow = 1))),
    cvx_train_loss = do.call(rbind, lapply(reps, function(x) x$trajectory$cvx_train_loss)),
    pgd_train_loss = do.call(rbind, lapply(reps, function(x) x$trajectory$pgd_train_loss)),
    selected_train_loss = do.call(rbind, lapply(reps, function(x) x$trajectory$selected_train_loss)),
    cvx_runtime_sec = do.call(rbind, lapply(reps, function(x) x$trajectory$cvx_runtime_sec)),
    pgd_runtime_sec = do.call(rbind, lapply(reps, function(x) x$trajectory$pgd_runtime_sec)),
    n_total = n_total,
    theta_dim = d
  )
}

run_policy_monte_carlo <- function(case_obj, policy, cfg, trajectory_file = NULL) {
  rep_idx <- seq_len(cfg$n_rep)
  alpha_id <- sprintf("%.1f", case_obj$alpha)
  d <- length(case_obj$sigmas)

  mat_proxy <- matrix(NA_real_, nrow = cfg$n_rep, ncol = cfg$n_total)
  mat_mse <- matrix(NA_real_, nrow = cfg$n_rep, ncol = cfg$n_total)
  mat_fro <- matrix(NA_real_, nrow = cfg$n_rep, ncol = cfg$n_total)
  mat_bures <- matrix(NA_real_, nrow = cfg$n_rep, ncol = cfg$n_total)

  chunk_size <- max(1L, as.integer(cfg$trajectory_chunk_size))
  part_files <- character(0)

  for (chunk_start in seq.int(1L, cfg$n_rep, by = chunk_size)) {
    chunk_end <- min(cfg$n_rep, chunk_start + chunk_size - 1L)
    ids <- rep_idx[chunk_start:chunk_end]

    rep_fun <- function(r) {
      seed_key <- sprintf("alpha_%s_%s_%s_rep_%d", alpha_id, case_obj$case_id, policy, r)
      rep_seed <- stable_seed(seed_key, base = cfg$seed_base)
      run_one_replicate_batch(case_obj, policy, cfg, rep_seed = rep_seed)
    }

    if (.Platform$OS.type != "windows" && cfg$n_workers > 1L) {
      reps <- parallel::mclapply(ids, rep_fun, mc.cores = cfg$n_workers)
    } else {
      reps <- lapply(ids, rep_fun)
    }

    for (j in seq_along(ids)) {
      rid <- ids[j]
      mat_proxy[rid, ] <- reps[[j]]$scaled_proxy
      mat_mse[rid, ] <- reps[[j]]$scaled_mse
      mat_fro[rid, ] <- reps[[j]]$scaled_fro_sq
      mat_bures[rid, ] <- reps[[j]]$scaled_bures_sq
    }

    if (isTRUE(cfg$save_full_trajectory) && !is.null(trajectory_file)) {
      part_idx <- length(part_files) + 1L
      part_file <- sub("\\.rds$", sprintf(".part%04d.rds", part_idx), trajectory_file)
      chunk_payload <- make_trajectory_chunk(ids, reps, cfg$n_total, d)
      saveRDS(chunk_payload, part_file)
      part_files <- c(part_files, part_file)
    }
  }

  out <- list(
    policy = policy,
    proxy = summarize_matrix(mat_proxy),
    mse = summarize_matrix(mat_mse),
    fro_sq = summarize_matrix(mat_fro),
    bures_sq = summarize_matrix(mat_bures)
  )

  if (isTRUE(cfg$save_full_matrices)) {
    out$matrices <- list(
      proxy = mat_proxy,
      mse = mat_mse,
      fro_sq = mat_fro,
      bures_sq = mat_bures
    )
  }

  if (isTRUE(cfg$save_full_trajectory) && !is.null(trajectory_file)) {
    traj_meta <- list(
      format = "chunked_replicate_trajectory_v1",
      case_id = case_obj$case_id,
      alpha = case_obj$alpha,
      policy = policy,
      n_rep = cfg$n_rep,
      n_total = cfg$n_total,
      n_trials = cfg$n_trials,
      theta_dim = d,
      part_files = normalizePath(part_files, winslash = "/", mustWork = FALSE)
    )
    saveRDS(traj_meta, trajectory_file)
    out$trajectory_artifact <- traj_meta
  }

  out
}

policy_style <- function(policy) {
  switch(policy,
    "uniform" = list(label = "Permutation", col = "gray40", lty = 1),
    "exact" = list(label = "Exact", col = "blue", lty = 1),
    "GI1" = list(label = "GI1", col = "red", lty = 2),
    "oracle_GI1" = list(label = "Oracle GI1", col = "darkgreen", lty = 4),
    list(label = policy, col = "black", lty = 1)
  )
}

policy_order <- function(case_result) {
  desired <- c("uniform", "exact", "GI1", "oracle_GI1")
  available <- names(case_result$results)
  desired[desired %in% available]
}

compute_case_comment <- function(case_result) {
  n_end <- case_result$config$n_total
  u <- case_result$results$uniform$proxy$mean[n_end]
  e <- case_result$results$exact$proxy$mean[n_end]
  g <- case_result$results$GI1$proxy$mean[n_end]
  o <- case_result$results$oracle_GI1$proxy$mean[n_end]
  exact_better <- is.finite(e) && is.finite(u) && (e < u)
  gi1_better <- is.finite(g) && is.finite(u) && (g < u)
  oracle_gi1_better <- is.finite(o) && is.finite(u) && (o < u)
  comment <- sprintf("Exact<Permutation: %s | GI1<Permutation: %s | OracleGI1<Permutation: %s",
                     ifelse(exact_better, "YES", "NO"),
                     ifelse(gi1_better, "YES", "NO"),
                     ifelse(oracle_gi1_better, "YES", "NO"))
  list(
    exact_better = exact_better,
    gi1_better = gi1_better,
    oracle_gi1_better = oracle_gi1_better,
    comment = comment,
    uniform_final = u,
    exact_final = e,
    gi1_final = g,
    oracle_gi1_final = o
  )
}

plot_case_scaled_proxy <- function(case_result, out_file, cfg) {
  n_idx <- x_idx_for_system(case_result$system, case_result$config$n_total, cfg)
  n_seq <- n_idx

  p_ord <- policy_order(case_result)
  vals <- unlist(lapply(p_ord, function(p) case_result$results[[p]]$proxy$mean[n_idx]), use.names = FALSE)
  vals <- c(vals, case_result$oracle_limit)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) vals <- c(0, 1)
  y_lim <- range(vals)
  span <- y_lim[2] - y_lim[1]
  pad <- if (span > 0) 0.05 * span else max(1e-6, abs(y_lim[1]) * 0.05 + 1e-6)
  y_lim <- c(max(0, y_lim[1] - pad), y_lim[2] + pad)

  png(out_file, width = 1100, height = 760, res = 120)
  par(mar = c(5, 5, 4, 2))

  first <- p_ord[1]
  st <- policy_style(first)
  plot(n_seq, case_result$results[[first]]$proxy$mean[n_idx], type = "l",
       col = st$col, lty = st$lty, lwd = 2.5,
       xlab = "Adaptive step (t)", ylab = expression(n[eff] %.% "Proxy risk"),
       main = sprintf("Case %02d | alpha=%.1f | %s | %s",
                      case_result$case_num,
                      case_result$alpha,
                      case_result$case_id,
                      case_result$loss),
       ylim = y_lim)

  for (p in p_ord[-1]) {
    stp <- policy_style(p)
    lines(n_seq, case_result$results[[p]]$proxy$mean[n_idx], col = stp$col, lty = stp$lty, lwd = 2.5)
  }

  abline(h = case_result$oracle_limit, col = "black", lty = 3, lwd = 2)
  grid(lty = "dotted", col = "gray85")

  lg <- c(vapply(p_ord, function(p) policy_style(p)$label, character(1)), "Oracle limit")
  lc <- c(vapply(p_ord, function(p) policy_style(p)$col, character(1)), "black")
  ll <- c(vapply(p_ord, function(p) policy_style(p)$lty, numeric(1)), 3)
  legend("topright", legend = lg, col = lc, lty = ll, lwd = c(rep(2.5, length(p_ord)), 2), bg = "white")

  mtext(sprintf("Window: t=%d..%d | n_eff = t x n_trials (n_trials=%d). %s",
                min(n_idx), max(n_idx), case_result$config$n_trials, case_result$comment),
        side = 1, line = 3.5, cex = 0.9)

  dev.off()
}

plot_alpha_relative_panel <- function(alpha_results, out_file, cfg) {
  ordered <- alpha_results[order(vapply(alpha_results, function(x) as.numeric(x$case_num), numeric(1)))]

  png(out_file, width = 2500, height = 2000, res = 140)
  par(mfrow = c(5, 3), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1))

  for (cr in ordered) {
    n_idx <- x_idx_for_system(cr$system, cr$config$n_total, cfg)
    n_seq <- n_idx
    oracle <- cr$oracle_limit
    p_ord <- policy_order(cr)

    rel_map <- lapply(p_ord, function(p) cr$results[[p]]$proxy$mean[n_idx] / oracle)
    names(rel_map) <- p_ord

    y_all <- unlist(rel_map, use.names = FALSE)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) y_all <- c(0.5, 1.5)
    y_hi <- max(y_all)
    span <- max(y_hi - 1, 0)
    pad <- if (span > 0) 0.08 * span else 0.08
    y_lim <- c(1, max(1.02, y_hi + pad))

    first <- p_ord[1]
    plot(n_seq, rel_map[[first]], type = "l", col = policy_style(first)$col, lty = policy_style(first)$lty, lwd = 1.8,
         xlab = "t", ylab = expression((n[eff] %.% "Proxy risk") / "Oracle limit"),
         ylim = y_lim,
         main = sprintf("%02d %s | %s\n%s",
                        cr$case_num, cr$library, cr$loss,
                        ifelse(cr$exact_better && cr$gi1_better && cr$oracle_gi1_better,
                               "Adaptive better than permutation", "Check ranking")),
         cex.main = 0.95)
    if (length(p_ord) > 1L) {
      for (p in p_ord[-1]) {
        st <- policy_style(p)
        lines(n_seq, rel_map[[p]], col = st$col, lty = st$lty, lwd = 1.8)
      }
    }
    abline(h = 1, col = "black", lty = 3, lwd = 1.4)
    grid(lty = "dotted", col = "gray85")
  }

  mtext(sprintf("Relative Oracle Panel (alpha=%.1f) | 1q:10..500, 2q:20..500", ordered[[1]]$alpha), outer = TRUE, cex = 1.5, line = 0.7)
  par(fig = c(0, 1, 0, 0.04), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", horiz = TRUE,
         legend = c("Permutation", "Exact", "GI1", "Oracle GI1", "Oracle ratio = 1"),
         col = c("gray40", "blue", "red", "darkgreen", "black"),
         lty = c(1, 1, 2, 4, 3), lwd = c(2, 2, 2, 2, 1.5), bty = "n")

  dev.off()
}

plot_alpha_metric_panel <- function(alpha_results, loss_target,
                                   metric_slot,
                                   oracle_getter,
                                   ylab,
                                   title_text,
                                   out_file,
                                   cfg) {
  keep <- alpha_results[vapply(alpha_results, function(x) identical(x$loss, loss_target), logical(1))]
  keep <- keep[order(vapply(keep, function(x) as.numeric(x$case_num), numeric(1)))]
  if (length(keep) == 0) return(invisible(NULL))

  png(out_file, width = 2100, height = 1700, res = 140)
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1))

  for (cr in keep) {
    n_idx <- x_idx_for_system(cr$system, cr$config$n_total, cfg)
    n_seq <- n_idx

    p_ord <- policy_order(cr)
    curve_map <- lapply(p_ord, function(p) cr$results[[p]][[metric_slot]]$mean[n_idx])
    names(curve_map) <- p_ord
    oracle_line <- oracle_getter(cr)

    y_all <- c(unlist(curve_map, use.names = FALSE), oracle_line)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) y_all <- c(0, 1)
    y_lim <- range(y_all)
    span <- y_lim[2] - y_lim[1]
    pad <- if (span > 0) 0.07 * span else 0.1
    y_lim <- c(max(0, y_lim[1] - pad), y_lim[2] + pad)

    first <- p_ord[1]
    plot(n_seq, curve_map[[first]], type = "l", col = policy_style(first)$col, lty = policy_style(first)$lty, lwd = 2,
         xlab = "Adaptive step (t)", ylab = ylab,
         ylim = y_lim,
         main = sprintf("%02d %s | %s\n%s", cr$case_num, cr$library, cr$loss, cr$comment),
         cex.main = 0.9)
    if (length(p_ord) > 1L) {
      for (p in p_ord[-1]) {
        st <- policy_style(p)
        lines(n_seq, curve_map[[p]], col = st$col, lty = st$lty, lwd = 2)
      }
    }
    abline(h = oracle_line, col = "black", lwd = 1.6, lty = 3)
    grid(lty = "dotted", col = "gray85")
  }

  mtext(sprintf("%s | 1q:10..500, 2q:20..500", title_text), outer = TRUE, cex = 1.35, line = 0.7)
  par(fig = c(0, 1, 0, 0.04), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", horiz = TRUE,
         legend = c("Permutation", "Exact", "GI1", "Oracle GI1", "Oracle limit"),
         col = c("gray40", "blue", "red", "darkgreen", "black"),
         lty = c(1, 1, 2, 4, 3), lwd = c(2, 2, 2, 2, 1.6), bty = "n")

  dev.off()
}

policy_for_selection <- function(policy) {
  if (policy %in% c("exact", "GI1")) return(policy)
  if (policy == "oracle_GI1") return("GI1")
  NA_character_
}

run_solver_pilot <- function(case_grid, alphas, cfg) {
  cvxr_available <- exists("CVXR_AVAILABLE", inherits = TRUE) && isTRUE(get("CVXR_AVAILABLE", inherits = TRUE))
  pilot_alpha <- cfg$pilot_alpha
  if (!(pilot_alpha %in% alphas)) pilot_alpha <- alphas[1]

  exp_by_case <- setNames(case_grid, vapply(case_grid, function(x) as.character(x$case_num), character(1)))
  pilot_cases <- as.character(cfg$pilot_cases)
  pilot_cases <- pilot_cases[pilot_cases %in% names(exp_by_case)]
  if (length(pilot_cases) == 0L) {
    pilot_cases <- names(exp_by_case)[seq_len(min(2L, length(exp_by_case)))]
  }

  rows <- list()
  row_idx <- 0L

  pilot_cfg <- cfg
  pilot_cfg$n_total <- as.integer(cfg$pilot_n_total)
  pilot_cfg$n_rep <- as.integer(cfg$pilot_n_rep)
  pilot_cfg$solver_selection_mode_resolved <- "always_dual"
  pilot_cfg$save_full_trajectory <- FALSE

  for (cid in pilot_cases) {
    exp <- exp_by_case[[cid]]
    case_obj <- build_case_object(exp, pilot_alpha, pilot_cfg)

    for (rep_id in seq_len(pilot_cfg$n_rep)) {
      rep_seed <- stable_seed(sprintf("pilot_alpha_%s_case_%s_rep_%d", pilot_alpha, cid, rep_id), base = cfg$seed_base)
      set.seed(rep_seed)

      Nab <- numeric(case_obj$lib$M)
      theta_hat <- rep(0, length(case_obj$sigmas))

      for (t in seq_len(pilot_cfg$n_total)) {
        if (t <= case_obj$n_init) {
          a_t <- case_obj$init_settings[t]
        } else {
          counts <- counts_by_setting(Nab, case_obj$lib$ab_df)
          G_hat <- case_obj$metric_fun(theta_hat)
          sel <- select_next_setting_metric_batch(
            theta_hat = theta_hat,
            counts = counts,
            S_ab = case_obj$S_ab,
            c_ab = case_obj$c_ab,
            ab_df = case_obj$lib$ab_df,
            G = G_hat,
            method = "GI1",
            n_trials = pilot_cfg$n_trials,
            ridge = pilot_cfg$ridge_selection,
            eps = pilot_cfg$fisher_eps
          )
          a_t <- sel$a_next
        }

        probs_a <- case_obj$prob_list_true[[a_t]]
        batch_counts <- as.integer(rmultinom(1, size = pilot_cfg$n_trials, prob = probs_a)[, 1])
        Nab <- increment_count_batch(Nab, a_t, batch_counts, case_obj$lib$ab_row)

        fit <- fit_theta_dual_compare(
          N = case_obj$N,
          sigmas = case_obj$sigmas,
          S_ab = case_obj$S_ab,
          c_ab = case_obj$c_ab,
          Nab = Nab,
          eta = pilot_cfg$eta_mle,
          cvx_solver = pilot_cfg$cvx_solver_preference,
          eps_log = pilot_cfg$nll_eps_log,
          verbose = FALSE,
          theta_start = theta_hat,
          truncate_density = TRUE,
          allow_eta_zero = (pilot_cfg$eta_mle <= 0)
        )

        if (mle_fit_success(fit)) {
          theta_hat <- fit$theta_hat
        }

        row_idx <- row_idx + 1L
        rows[[row_idx]] <- data.frame(
          alpha = pilot_alpha,
          case_num = as.integer(exp$case_num),
          case_id = exp$case_id,
          rep = rep_id,
          t = t,
          cvx_status = as.character(fit$cvx_status %||% NA_character_),
          pgd_status = as.character(fit$pgd_status %||% NA_character_),
          cvx_train_loss = as.numeric(fit$cvx_train_loss %||% NA_real_),
          pgd_train_loss = as.numeric(fit$pgd_train_loss %||% NA_real_),
          cvx_runtime_sec = as.numeric(fit$cvx_runtime_sec %||% NA_real_),
          pgd_runtime_sec = as.numeric(fit$pgd_runtime_sec %||% NA_real_),
          selected_solver = as.character(fit$solver_dual_selected %||% fit$solver_used %||% NA_character_),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  step_df <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()

  is_ok <- function(x) x %in% c("optimal", "optimal_inaccurate", "trivial")
  cvx_ok <- if (nrow(step_df) > 0) mean(is_ok(step_df$cvx_status), na.rm = TRUE) else NA_real_
  pgd_ok <- if (nrow(step_df) > 0) mean(is_ok(step_df$pgd_status), na.rm = TRUE) else NA_real_

  mean_cvx_loss <- if (nrow(step_df) > 0) mean(step_df$cvx_train_loss, na.rm = TRUE) else NA_real_
  mean_pgd_loss <- if (nrow(step_df) > 0) mean(step_df$pgd_train_loss, na.rm = TRUE) else NA_real_
  mean_cvx_time <- if (nrow(step_df) > 0) mean(step_df$cvx_runtime_sec, na.rm = TRUE) else NA_real_
  mean_pgd_time <- if (nrow(step_df) > 0) mean(step_df$pgd_runtime_sec, na.rm = TRUE) else NA_real_

  chosen <- "CVX"
  reason <- "default_to_cvx"

  if (!cvxr_available) {
    chosen <- "PGD"
    reason <- "cvxr_unavailable"
  } else if (is.finite(cvx_ok) && is.finite(pgd_ok) && cvx_ok < 0.5 && pgd_ok > cvx_ok + 0.1) {
    chosen <- "PGD"
    reason <- "cvx_fail_rate_high"
  } else if (is.finite(mean_pgd_loss) && is.finite(mean_cvx_loss) && mean_pgd_loss < mean_cvx_loss - 1e-6) {
    chosen <- "PGD"
    reason <- "pgd_lower_train_loss"
  } else if (is.finite(mean_pgd_time) && is.finite(mean_cvx_time) &&
             mean_pgd_time > 1.25 * mean_cvx_time &&
             (!is.finite(mean_pgd_loss) || !is.finite(mean_cvx_loss) || mean_cvx_loss <= mean_pgd_loss + 1e-6)) {
    chosen <- "CVX"
    reason <- "pgd_slower_without_loss_gain"
  } else if (is.finite(cvx_ok) && is.finite(pgd_ok) && cvx_ok >= pgd_ok) {
    chosen <- "CVX"
    reason <- "cvx_success_rate_not_worse"
  }

  summary_df <- data.frame(
    pilot_alpha = pilot_alpha,
    pilot_cases = paste(pilot_cases, collapse = ","),
    n_rep = cfg$pilot_n_rep,
    n_total = cfg$pilot_n_total,
    cvx_solver_preference = cfg$cvx_solver_preference,
    cvx_success_rate = cvx_ok,
    pgd_success_rate = pgd_ok,
    mean_cvx_train_loss = mean_cvx_loss,
    mean_pgd_train_loss = mean_pgd_loss,
    mean_cvx_runtime_sec = mean_cvx_time,
    mean_pgd_runtime_sec = mean_pgd_time,
    chosen_solver = chosen,
    decision_reason = reason,
    stringsAsFactors = FALSE
  )

  list(step_df = step_df, summary_df = summary_df, chosen_solver = chosen, reason = reason)
}

bind_rows_or_empty <- function(x) {
  if (length(x) == 0) return(data.frame())
  do.call(rbind, x)
}

## ============================================================================
## Configuration
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_profile <- parse_cli_value("--config_profile=", args = args, default = "full")
config_file <- parse_cli_value("--config_file=", args = args, default = NULL)
cfg <- get_caseI_caseII_config(profile = config_profile, config_file = config_file, args = args)

# Backward-compatible alias for legacy CLI flag used in older runs.
ridge_cli <- parse_cli_value("--ridge=", args = args, default = NA_character_)
if (!is.na(ridge_cli)) cfg$ridge_selection <- as_num_or(ridge_cli, cfg$ridge_selection)

output_dir <- parse_cli_value("--output_dir=", args = args, default = cfg$output_dir %||% "")
if (!nzchar(output_dir)) {
  alpha_tag_join <- paste(vapply(cfg$alphas, function(a) gsub("\\.", "p", format(a, nsmall = 1, trim = TRUE)), character(1)), collapse = "_")
  output_dir <- sprintf("results_caseI_caseII_alpha%s_nrep%d_n%d_trials%d", alpha_tag_join, cfg$n_rep, cfg$n_total, cfg$n_trials)
}
alphas <- cfg$alphas
alpha_names <- vapply(alphas, alpha_case_name, FUN.VALUE = character(1))

if (cfg$n_rep < 1L) stop("n_rep must be >= 1")
if (cfg$n_total < 1L) stop("n_total must be >= 1")
if (cfg$n_trials < 1L) stop("n_trials must be >= 1")
if (cfg$n_workers < 1L) cfg$n_workers <- 1L
if (cfg$trajectory_chunk_size < 1L) cfg$trajectory_chunk_size <- 1L

options(qt_pgd_max_iter = as.integer(cfg$pgd_max_iter))

cat("Configuration:\n")
cat(sprintf("  config_profile = %s\n", config_profile))
cat(sprintf("  n_rep = %d\n", cfg$n_rep))
cat(sprintf("  n_total (adaptive steps) = %d\n", cfg$n_total))
cat(sprintf("  n_trials (multinomial size) = %d\n", cfg$n_trials))
cat(sprintf("  check_every = %d\n", cfg$check_every))
cat(sprintf("  eta_mle = %.2e\n", cfg$eta_mle))
cat(sprintf("  ridge_selection = %.2e\n", cfg$ridge_selection))
cat(sprintf("  fisher_eps = %.2e\n", cfg$fisher_eps))
cat(sprintf("  nll_eps_log = %.2e\n", cfg$nll_eps_log))
cat(sprintf("  bures_tol = %.2e\n", cfg$bures_tol))
cat(sprintf("  solver_selection_mode = %s\n", cfg$solver_selection_mode))
cat(sprintf("  cvx_solver_preference = %s\n", cfg$cvx_solver_preference))
cat(sprintf("  solver (legacy) = %s\n", cfg$solver))
cat(sprintf("  oracle_solver = %s\n", cfg$oracle_solver))
cat(sprintf("  n_workers = %d\n", cfg$n_workers))
cat(sprintf("  pgd_max_iter = %d\n", cfg$pgd_max_iter))
cat(sprintf("  save_full_matrices = %s\n", ifelse(cfg$save_full_matrices, "TRUE", "FALSE")))
cat(sprintf("  save_full_trajectory = %s\n", ifelse(cfg$save_full_trajectory, "TRUE", "FALSE")))
cat(sprintf("  trajectory_chunk_size = %d\n", cfg$trajectory_chunk_size))
cat(sprintf("  pilot_n_rep = %d\n", cfg$pilot_n_rep))
cat(sprintf("  pilot_n_total = %d\n", cfg$pilot_n_total))
cat(sprintf("  pilot_cases = %s\n", paste(cfg$pilot_cases, collapse = ",")))
cat(sprintf("  pilot_alpha = %.3f\n", cfg$pilot_alpha))
cat(sprintf("  plot_1q_start = %d\n", cfg$plot_1q_start))
cat(sprintf("  plot_2q_start = %d\n", cfg$plot_2q_start))
cat(sprintf("  plot_t_end = %d\n", cfg$plot_t_end))
cat(sprintf("  baseline_policy_mode = %s\n", cfg$baseline_policy_mode))
cat(sprintf("  resume = %s\n", ifelse(cfg$resume, "TRUE", "FALSE")))
cat(sprintf("  output_dir = %s\n\n", output_dir))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
plots_dir <- file.path(output_dir, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
trajectories_dir <- file.path(output_dir, "trajectories")
if (!dir.exists(trajectories_dir)) dir.create(trajectories_dir, recursive = TRUE)
checkpoint_file <- file.path(output_dir, "checkpoint_latest.rds")
complete_flag <- file.path(output_dir, "RUN_COMPLETE.flag")
pilot_summary_file <- file.path(output_dir, "solver_pilot_summary.csv")
pilot_decision_file <- file.path(output_dir, "solver_pilot_decision.csv")
pilot_steps_file <- file.path(output_dir, "solver_pilot_step_details.csv")

## ============================================================================
## Solver strategy resolution (pilot / dual / fallback)
## ============================================================================

case_grid <- get_case_grid()

if (cfg$solver_selection_mode == "always_dual") {
  cfg$solver_selection_mode_resolved <- "always_dual"
  cfg$production_solver <- "DUAL"
} else if (cfg$solver_selection_mode == "cvx_fallback_pgd") {
  cfg$solver_selection_mode_resolved <- "cvx_fallback_pgd"
  cfg$production_solver <- "CVX"
} else if (cfg$solver_selection_mode == "pilot_choose") {
  if (cfg$resume && file.exists(pilot_decision_file) && file.exists(pilot_summary_file)) {
    pilot_decision <- tryCatch(read.csv(pilot_decision_file, stringsAsFactors = FALSE), error = function(e) data.frame())
    if (nrow(pilot_decision) > 0 && "chosen_solver" %in% names(pilot_decision)) {
      chosen <- as.character(pilot_decision$chosen_solver[1])
      if (!chosen %in% c("CVX", "PGD")) chosen <- "CVX"
      cfg$production_solver <- chosen
      cfg$solver_selection_mode_resolved <- if (chosen == "PGD") "single_pgd" else "cvx_fallback_pgd"
      cat(sprintf("Loaded existing pilot decision: %s\n", chosen))
    } else {
      pilot_res <- run_solver_pilot(case_grid, alphas, cfg)
      write.csv(pilot_res$summary_df, pilot_summary_file, row.names = FALSE)
      write.csv(pilot_res$step_df, pilot_steps_file, row.names = FALSE)
      write.csv(data.frame(chosen_solver = pilot_res$chosen_solver,
                           decision_reason = pilot_res$reason,
                           stringsAsFactors = FALSE), pilot_decision_file, row.names = FALSE)
      cfg$production_solver <- pilot_res$chosen_solver
      cfg$solver_selection_mode_resolved <- if (cfg$production_solver == "PGD") "single_pgd" else "cvx_fallback_pgd"
      cat(sprintf("Pilot completed: chosen solver = %s (%s)\n", pilot_res$chosen_solver, pilot_res$reason))
    }
  } else {
    pilot_res <- run_solver_pilot(case_grid, alphas, cfg)
    write.csv(pilot_res$summary_df, pilot_summary_file, row.names = FALSE)
    write.csv(pilot_res$step_df, pilot_steps_file, row.names = FALSE)
    write.csv(data.frame(chosen_solver = pilot_res$chosen_solver,
                         decision_reason = pilot_res$reason,
                         stringsAsFactors = FALSE), pilot_decision_file, row.names = FALSE)
    cfg$production_solver <- pilot_res$chosen_solver
    cfg$solver_selection_mode_resolved <- if (cfg$production_solver == "PGD") "single_pgd" else "cvx_fallback_pgd"
    cat(sprintf("Pilot completed: chosen solver = %s (%s)\n", pilot_res$chosen_solver, pilot_res$reason))
  }
} else {
  cfg$solver_selection_mode_resolved <- "cvx_fallback_pgd"
  cfg$production_solver <- "CVX"
}

cat(sprintf("Resolved solver execution mode: %s\n", cfg$solver_selection_mode_resolved))
cat(sprintf("Resolved production solver: %s\n\n", cfg$production_solver))

## ============================================================================
## Main Run
## ============================================================================

policies <- c("uniform", "exact", "GI1", "oracle_GI1")

all_results <- list()
oracle_rows <- list()
curve_rows <- list()
eq23_rows <- list()
comment_rows <- list()
plot_rows <- list()
trajectory_rows <- list()

if (cfg$resume && file.exists(complete_flag)) {
  cat("Completion flag exists; run already finished for this output_dir.\n")
  q(save = "no", status = 0)
}

if (cfg$resume && file.exists(checkpoint_file)) {
  cat(sprintf("Loading checkpoint: %s\n", checkpoint_file))
  ck <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
  if (!is.null(ck) && is.list(ck)) {
    all_results <- ck$all_results %||% list()
    oracle_rows <- ck$oracle_rows %||% list()
    curve_rows <- ck$curve_rows %||% list()
    eq23_rows <- ck$eq23_rows %||% list()
    comment_rows <- ck$comment_rows %||% list()
    plot_rows <- ck$plot_rows %||% list()
    trajectory_rows <- ck$trajectory_rows %||% list()
    if (!is.list(all_results)) all_results <- list()
    if (!is.list(oracle_rows)) oracle_rows <- list()
    if (!is.list(curve_rows)) curve_rows <- list()
    if (!is.list(eq23_rows)) eq23_rows <- list()
    if (!is.list(comment_rows)) comment_rows <- list()
    if (!is.list(plot_rows)) plot_rows <- list()
    if (!is.list(trajectory_rows)) trajectory_rows <- list()
    cat(sprintf("Checkpoint loaded. Completed cases: %d\n", length(all_results)))
  } else {
    cat("Checkpoint could not be parsed; starting from scratch.\n")
  }
}

run_start <- Sys.time()

for (ai in seq_along(alphas)) {
  alpha <- alphas[ai]
  alpha_name <- alpha_names[ai]
  a_tag <- alpha_tag(alpha)

  cat(sprintf("\n============================\n"))
  cat(sprintf("Running %s (alpha=%.1f)\n", alpha_name, alpha))
  cat(sprintf("============================\n"))

  alpha_case_results <- list()
  for (exp in case_grid) {
    key_existing <- sprintf("alpha_%s_%s", a_tag, exp$case_id)
    if (!is.null(all_results[[key_existing]])) {
      alpha_case_results[[length(alpha_case_results) + 1L]] <- all_results[[key_existing]]
    }
  }

  for (exp in case_grid) {
    key_existing <- sprintf("alpha_%s_%s", a_tag, exp$case_id)
    if (!is.null(all_results[[key_existing]])) {
      cat(sprintf("\n[alpha=%.1f] Case %02d: %s (already completed; skip)\n",
                  alpha, exp$case_num, exp$case_id))
      next
    }

    cat(sprintf("\n[alpha=%.1f] Case %02d: %s\n", alpha, exp$case_num, exp$case_id))
    t_case <- Sys.time()

    case_obj <- build_case_object(exp, alpha, cfg)

    results_policy <- list()
    for (policy in policies) {
      cat(sprintf("  - policy %-11s ... ", policy))
      flush.console()
      t_pol <- Sys.time()

      traj_file <- NULL
      if (isTRUE(cfg$save_full_trajectory)) {
        traj_file <- file.path(
          trajectories_dir,
          sprintf("alpha_%s_case%02d_%s_policy_%s_trajectory.rds", a_tag, case_obj$case_num, case_obj$case_id, policy)
        )
      }

      mc <- run_policy_monte_carlo(case_obj, policy, cfg, trajectory_file = traj_file)
      results_policy[[policy]] <- mc

      if (!is.null(mc$trajectory_artifact)) {
        trajectory_rows[[length(trajectory_rows) + 1L]] <- data.frame(
          alpha = alpha,
          alpha_name = alpha_name,
          case_num = case_obj$case_num,
          case_id = case_obj$case_id,
          policy = policy,
          path = traj_file,
          n_settings = case_obj$lib$k,
          n_init = case_obj$n_init,
          n_rep = cfg$n_rep,
          n_total = cfg$n_total,
          n_trials = cfg$n_trials,
          n_parts = length(mc$trajectory_artifact$part_files),
          stringsAsFactors = FALSE
        )
      }

      dt <- as.numeric(difftime(Sys.time(), t_pol, units = "secs"))
      cat(sprintf("done (%.1fs)\n", dt))
    }

    case_result <- list(
      alpha = alpha,
      alpha_name = alpha_name,
      case_num = case_obj$case_num,
      case_id = case_obj$case_id,
      system = case_obj$system,
      library = case_obj$library,
      loss = case_obj$loss,
      library_full = case_obj$library_full,
      oracle_limit = case_obj$oracle_limit,
      oracle_mse_limit = case_obj$oracle_mse_limit,
      uniform_limit = case_obj$uniform_limit,
      ratio_uniform_over_oracle = case_obj$ratio_uniform_over_oracle,
      oracle_status = case_obj$oracle_status,
      N = case_obj$N,
      config = cfg,
      results = results_policy
    )

    cmp <- compute_case_comment(case_result)
    case_result$exact_better <- cmp$exact_better
    case_result$gi1_better <- cmp$gi1_better
    case_result$oracle_gi1_better <- cmp$oracle_gi1_better
    case_result$comment <- cmp$comment

    key <- sprintf("alpha_%s_%s", a_tag, case_result$case_id)
    all_results[[key]] <- case_result
    alpha_case_results[[length(alpha_case_results) + 1L]] <- case_result

    oracle_rows[[length(oracle_rows) + 1L]] <- data.frame(
      alpha = alpha,
      alpha_name = alpha_name,
      case_num = case_result$case_num,
      case_id = case_result$case_id,
      system = case_result$system,
      library = case_result$library,
      loss = case_result$loss,
      oracle_limit = case_result$oracle_limit,
      oracle_mse_limit = case_result$oracle_mse_limit,
      uniform_limit = case_result$uniform_limit,
      ratio_uniform_over_oracle = case_result$ratio_uniform_over_oracle,
      oracle_status = case_result$oracle_status,
      stringsAsFactors = FALSE
    )

    n_seq <- seq_len(cfg$n_total)
    n_eff <- n_seq * cfg$n_trials
    for (policy in policies) {
      rp <- case_result$results[[policy]]
      curve_rows[[length(curve_rows) + 1L]] <- data.frame(
        alpha = alpha,
        alpha_name = alpha_name,
        case_num = case_result$case_num,
        case_id = case_result$case_id,
        system = case_result$system,
        library = case_result$library,
        loss = case_result$loss,
        policy = policy,
        n_step = n_seq,
        n_eff = n_eff,
        mean_scaled_proxy = rp$proxy$mean,
        se_scaled_proxy = rp$proxy$se,
        mean_scaled_mse = rp$mse$mean,
        mean_scaled_fro_sq = rp$fro_sq$mean,
        mean_scaled_bures_sq = rp$bures_sq$mean,
        stringsAsFactors = FALSE
      )
    }

    final_idx <- cfg$n_total
    if (case_result$loss %in% c("frobenius", "bures")) {
      for (policy in policies) {
        rp <- case_result$results[[policy]]
        if (case_result$loss == "frobenius") {
          lhs_final <- rp$mse$mean[final_idx]
          lhs_tail <- mean(rp$mse$mean[pmax(1, final_idx - 19):final_idx], na.rm = TRUE)
          rhs <- case_result$oracle_mse_limit
          metric_name <- "scaled_frobenius_mse"
        } else {
          lhs_final <- rp$bures_sq$mean[final_idx]
          lhs_tail <- mean(rp$bures_sq$mean[pmax(1, final_idx - 19):final_idx], na.rm = TRUE)
          rhs <- case_result$oracle_limit
          metric_name <- "scaled_bures_squared"
        }
        eq23_rows[[length(eq23_rows) + 1L]] <- data.frame(
          alpha = alpha,
          alpha_name = alpha_name,
          case_num = case_result$case_num,
          case_id = case_result$case_id,
          system = case_result$system,
          library = case_result$library,
          loss = case_result$loss,
          policy = policy,
          metric = metric_name,
          lhs_final = lhs_final,
          lhs_tail_mean = lhs_tail,
          rhs_oracle = rhs,
          rel_gap_final = ifelse(is.finite(rhs) && rhs != 0, (lhs_final - rhs) / rhs, NA_real_),
          rel_gap_tail = ifelse(is.finite(rhs) && rhs != 0, (lhs_tail - rhs) / rhs, NA_real_),
          stringsAsFactors = FALSE
        )
      }
    }

    comment_rows[[length(comment_rows) + 1L]] <- data.frame(
      alpha = alpha,
      alpha_name = alpha_name,
      case_num = case_result$case_num,
      case_id = case_result$case_id,
      system = case_result$system,
      library = case_result$library,
      loss = case_result$loss,
      exact_better = case_result$exact_better,
      gi1_better = case_result$gi1_better,
      oracle_gi1_better = case_result$oracle_gi1_better,
      uniform_final = cmp$uniform_final,
      exact_final = cmp$exact_final,
      gi1_final = cmp$gi1_final,
      oracle_gi1_final = cmp$oracle_gi1_final,
      comment = case_result$comment,
      stringsAsFactors = FALSE
    )

    case_plot <- file.path(plots_dir, sprintf("alpha_%s_case%02d_%s_scaled_proxy.png",
                                               a_tag, case_result$case_num, case_result$case_id))
    plot_case_scaled_proxy(case_result, case_plot, cfg)
    plot_rows[[length(plot_rows) + 1L]] <- data.frame(
      alpha = alpha,
      case_num = case_result$case_num,
      case_id = case_result$case_id,
      plot_type = "case_scaled_proxy",
      path = case_plot,
      stringsAsFactors = FALSE
    )

    saveRDS(
      list(
        all_results = all_results,
        oracle_rows = oracle_rows,
        curve_rows = curve_rows,
        eq23_rows = eq23_rows,
        comment_rows = comment_rows,
        plot_rows = plot_rows,
        trajectory_rows = trajectory_rows,
        cfg = cfg,
        alpha = alpha,
        case_num = case_result$case_num
      ),
      file.path(output_dir, "checkpoint_latest.rds")
    )

    dt_case <- as.numeric(difftime(Sys.time(), t_case, units = "secs"))
    cat(sprintf("  Case done in %.1fs | %s\n", dt_case, case_result$comment))
  }

  alpha_case_results <- alpha_case_results[order(vapply(alpha_case_results, function(x) as.numeric(x$case_num), numeric(1)))]

  alpha_rel_plot <- file.path(plots_dir, sprintf("alpha_%s_all15_relative_oracle_proxy.png", a_tag))
  plot_alpha_relative_panel(alpha_case_results, alpha_rel_plot, cfg)
  plot_rows[[length(plot_rows) + 1L]] <- data.frame(
    alpha = alpha,
    case_num = NA_integer_,
    case_id = "ALL_15",
    plot_type = "all15_relative_oracle_proxy",
    path = alpha_rel_plot,
    stringsAsFactors = FALSE
  )

  alpha_fro_plot <- file.path(plots_dir, sprintf("alpha_%s_frobenius_mse_vs_oracle.png", a_tag))
  plot_alpha_metric_panel(
    alpha_case_results,
    loss_target = "frobenius",
    metric_slot = "mse",
    oracle_getter = function(cr) cr$oracle_mse_limit,
    ylab = expression(n[eff] %.% "Frobenius MSE"),
    title_text = sprintf("Frobenius MSE vs Oracle (alpha=%.1f)", alpha),
    out_file = alpha_fro_plot,
    cfg = cfg
  )
  plot_rows[[length(plot_rows) + 1L]] <- data.frame(
    alpha = alpha,
    case_num = NA_integer_,
    case_id = "FROBENIUS_5",
    plot_type = "frobenius_mse_vs_oracle",
    path = alpha_fro_plot,
    stringsAsFactors = FALSE
  )

  alpha_bures_plot <- file.path(plots_dir, sprintf("alpha_%s_bures_sq_vs_oracle.png", a_tag))
  plot_alpha_metric_panel(
    alpha_case_results,
    loss_target = "bures",
    metric_slot = "bures_sq",
    oracle_getter = function(cr) cr$oracle_limit,
    ylab = expression(n[eff] %.% "Bures"^2),
    title_text = sprintf("Bures Squared Loss vs Oracle (alpha=%.1f)", alpha),
    out_file = alpha_bures_plot,
    cfg = cfg
  )
  plot_rows[[length(plot_rows) + 1L]] <- data.frame(
    alpha = alpha,
    case_num = NA_integer_,
    case_id = "BURES_5",
    plot_type = "bures_sq_vs_oracle",
    path = alpha_bures_plot,
    stringsAsFactors = FALSE
  )
}

run_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))

oracle_df <- bind_rows_or_empty(oracle_rows)
curve_df <- bind_rows_or_empty(curve_rows)
eq23_df <- bind_rows_or_empty(eq23_rows)
comment_df <- bind_rows_or_empty(comment_rows)
plot_df <- bind_rows_or_empty(plot_rows)
trajectory_df <- bind_rows_or_empty(trajectory_rows)

write.csv(oracle_df, file.path(output_dir, "oracle_limits_caseI_caseII.csv"), row.names = FALSE)
write.csv(curve_df, file.path(output_dir, "curve_summary_caseI_caseII.csv"), row.names = FALSE)
write.csv(eq23_df, file.path(output_dir, "eq23_verification_frobenius_bures.csv"), row.names = FALSE)
write.csv(comment_df, file.path(output_dir, "policy_comparison_comments.csv"), row.names = FALSE)
write.csv(plot_df, file.path(output_dir, "plot_manifest.csv"), row.names = FALSE)
write.csv(trajectory_df, file.path(output_dir, "trajectory_manifest.csv"), row.names = FALSE)

saveRDS(cfg, file.path(output_dir, "run_config.rds"))
saveRDS(all_results, file.path(output_dir, "full_results_caseI_caseII.rds"))
writeLines(as.character(Sys.time()), complete_flag)

cat("\n==============================================================\n")
cat("Run complete.\n")
cat(sprintf("Elapsed: %.1f seconds\n", run_seconds))
cat(sprintf("Saved: %s\n", normalizePath(output_dir)))
cat("Key outputs:\n")
cat(sprintf("  - %s\n", file.path(output_dir, "oracle_limits_caseI_caseII.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "curve_summary_caseI_caseII.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "eq23_verification_frobenius_bures.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "policy_comparison_comments.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "trajectory_manifest.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "plot_manifest.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "plots")))
cat(sprintf("  - %s\n", file.path(output_dir, "trajectories")))
cat("==============================================================\n")
