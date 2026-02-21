## ============================================================================
## 10_run_case11_bures_rep1000_parallel.R
## Focused run for Case 11 (2q_A_bures) with parallel replicas.
## Policies: permutation-baseline(key='uniform') / exact / GI1 / oracle_GI1
## ============================================================================

script_dir <- tryCatch({
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  this_file <- sys.frame(1)$ofile
  if (!is.null(this_file) && this_file != "") return(dirname(normalizePath(this_file)))
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
}

source("01_utilities.R")
source("02_state_basis.R")
source("03_measurement_library_1q.R")
source("04_measurement_library_2q.R")
source("05_mle_cvxr_solver.R")
source("06_fisher_and_metrics.R")
source(file.path(getwd(), "00_simulation_parameters.R"))

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

solve_optimal_design <- function(I_list, G, solver = "SCS") {
  if (!exists("CVXR_AVAILABLE", inherits = TRUE) || !isTRUE(get("CVXR_AVAILABLE", inherits = TRUE))) {
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
      scores[i] <- tryCatch(sum(diag(G %*% safe_inverse(J_cand, ridge = 0))), error = function(e) Inf)
    }
    a_next <- as.integer(names(scores)[which.min(scores)])
  } else {
    J_inv <- tryCatch(safe_inverse(J_reg, ridge = 0), error = function(e) NULL)
    if (is.null(J_inv)) stop("GI1 selection failed: singular accumulated Fisher matrix and ridge=0.")
    for (i in seq_along(I_list)) {
      Ia <- I_list[[i]]
      scores[i] <- as.numeric(n_trials) * sum(diag(G %*% J_inv %*% Ia %*% J_inv))
    }
    a_next <- as.integer(names(scores)[which.max(scores)])
  }

  list(a_next = a_next, scores = scores)
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

mle_fit_success <- function(fit) fit$status %in% c("optimal", "optimal_inaccurate", "trivial")

extract_fit_diag <- function(fit, cvx_solver) {
  if (cvx_solver %in% c("PGD", "AUTO", "AUTO_PGD")) cvx_solver <- "SCS"
  cvx_status <- solver_map_status(fit, cvx_solver, fallback = if (identical(fit$solver_used, cvx_solver)) fit$status else NA_character_)
  pgd_status <- solver_map_status(fit, "PGD", fallback = if (identical(fit$solver_used, "PGD")) fit$status else NA_character_)
  cvx_loss <- solver_map_val(fit, cvx_solver, "train_loss", fallback = NA_real_)
  pgd_loss <- solver_map_val(fit, "PGD", "train_loss", fallback = NA_real_)
  if (!is.finite(cvx_loss) && identical(fit$solver_used, cvx_solver)) cvx_loss <- as.numeric(fit$value)
  if (!is.finite(pgd_loss) && identical(fit$solver_used, "PGD")) pgd_loss <- as.numeric(fit$value)
  selected_loss <- if (identical(fit$solver_used, "PGD")) pgd_loss else cvx_loss
  if (!is.finite(selected_loss)) selected_loss <- as.numeric(fit$value)
  list(
    solver_used = as.character(fit$solver_used),
    cvx_status = as.character(cvx_status),
    pgd_status = as.character(pgd_status),
    cvx_train_loss = as.numeric(cvx_loss),
    pgd_train_loss = as.numeric(pgd_loss),
    selected_train_loss = as.numeric(selected_loss)
  )
}

fit_mle_step <- function(case_obj, Nab, theta_hat, cfg) {
  cvx_backend <- cfg$cvx_solver_preference
  if (cvx_backend %in% c("PGD", "AUTO", "AUTO_PGD")) cvx_backend <- "SCS"

  if (cfg$solver_selection_mode == "always_dual") {
    fit <- fit_theta_dual_compare(
      N = case_obj$N,
      sigmas = case_obj$sigmas,
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
      selected_train_loss = as.numeric(fit$selected_train_loss %||% NA_real_)
    )
    return(list(fit = fit, diag = diag))
  }

  if (cfg$solver_selection_mode == "single_pgd") {
    fit <- fit_theta_cvxr(
      N = case_obj$N,
      sigmas = case_obj$sigmas,
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
    return(list(fit = fit, diag = extract_fit_diag(fit, cvx_backend)))
  }

  fit_cvx <- fit_theta_cvxr(
    N = case_obj$N,
    sigmas = case_obj$sigmas,
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
  )

  if (mle_fit_success(fit_cvx)) {
    return(list(fit = fit_cvx, diag = extract_fit_diag(fit_cvx, cvx_backend)))
  }

  fit_pgd <- fit_theta_cvxr(
    N = case_obj$N,
    sigmas = case_obj$sigmas,
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

  diag <- extract_fit_diag(fit_pgd, cvx_backend)
  diag$cvx_status <- as.character(fit_cvx$status)
  diag$cvx_train_loss <- as.numeric(fit_cvx$value)
  list(fit = fit_pgd, diag = diag)
}

build_case11 <- function(alpha, cfg) {
  lib <- build_library_2q_A()
  basis <- build_pauli_product_basis_2q()
  sigmas <- basis$sigmas
  N <- lib$N

  rho_pure <- fixed_pure_state_2q_phi_phase(phi = 0.2)
  rho_true <- mix_with_maximally_mixed(rho_pure, alpha)
  theta_true <- theta_from_rho(rho_true, sigmas)

  metric_fun <- function(theta) metric_bures(theta, sigmas, N, tol = cfg$bures_tol)
  G_true <- metric_fun(theta_true)

  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)
  I_list <- fisher_info_all_settings(theta_true, sc$S_ab, sc$c_ab, lib$ab_df)
  k <- length(I_list)

  I_uniform <- Reduce("+", I_list) / k
  uniform_limit <- proxy_risk(G_true, I_uniform, ridge = cfg$ridge_selection)

  pi_opt <- solve_optimal_design_mirror(I_list, G_true, ridge = cfg$ridge_selection)$pi
  I_opt <- Reduce("+", lapply(seq_len(k), function(i) pi_opt[i] * I_list[[i]]))
  oracle_limit <- proxy_risk(G_true, I_opt, ridge = cfg$ridge_selection)

  list(
    alpha = alpha,
    lib = lib,
    sigmas = sigmas,
    N = N,
    rho_true = rho_true,
    theta_true = theta_true,
    metric_fun = metric_fun,
    G_true = G_true,
    S_ab = sc$S_ab,
    c_ab = sc$c_ab,
    prob_list_true = born_probs_list(rho_true, lib$Q_list),
    init_settings = 1:15,
    n_init = 15L,
    oracle_limit = oracle_limit,
    uniform_limit = uniform_limit
  )
}

run_one_replicate <- function(case_obj, policy, cfg, rep_seed) {
  set.seed(rep_seed)

  lib <- case_obj$lib
  N <- case_obj$N
  sigmas <- case_obj$sigmas
  d <- length(sigmas)

  Nab <- numeric(lib$M)
  theta_hat <- rep(0, d)
  rho_hat_cur <- diag(N) / N
  scaled_proxy <- numeric(cfg$n_total)
  scaled_bures_sq <- numeric(cfg$n_total)

  if (isTRUE(cfg$save_full_trajectory)) {
    action_idx <- integer(cfg$n_total)
    action_label <- character(cfg$n_total)
    theta_path <- matrix(NA_real_, nrow = cfg$n_total, ncol = d)
    solver_used <- character(cfg$n_total)
    cvx_status <- character(cfg$n_total)
    pgd_status <- character(cfg$n_total)
    cvx_train_loss <- rep(NA_real_, cfg$n_total)
    pgd_train_loss <- rep(NA_real_, cfg$n_total)
    selected_train_loss <- rep(NA_real_, cfg$n_total)
  }

  for (t in seq_len(cfg$n_total)) {
    if (t <= case_obj$n_init) {
      a_t <- case_obj$init_settings[t]
    } else if (policy == "uniform") {
      if (identical(cfg$baseline_policy_mode, "permutation_cyclic")) {
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
          n_trials = cfg$n_trials,
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
          n_trials = cfg$n_trials,
          ridge = cfg$ridge_selection,
          eps = cfg$fisher_eps
        )
      }
      a_t <- sel$a_next
    }

    probs_a <- case_obj$prob_list_true[[a_t]]
    batch_counts <- as.integer(rmultinom(1, size = cfg$n_trials, prob = probs_a)[, 1])
    Nab <- increment_count_batch(Nab, a_t, batch_counts, lib$ab_row)

    fit_out <- fit_mle_step(case_obj, Nab, theta_hat, cfg)
    fit <- fit_out$fit

    if (mle_fit_success(fit)) {
      theta_hat <- fit$theta_hat
      if (!is.null(fit$rho_hat) && is.matrix(fit$rho_hat)) {
        rho_hat_cur <- hermitianize(fit$rho_hat)
      } else {
        rho_hat_cur <- hermitianize(rho_of_theta(theta_hat, sigmas, N))
      }
    }
    rho_hat_cur <- enforce_rho_floor(rho_hat_cur, eta = 0, sigmas = sigmas)$rho_hat
    bures_sq <- bures_distance_density(rho_hat_cur, case_obj$rho_true, tol = cfg$bures_tol)^2

    counts_now <- counts_by_setting(Nab, lib$ab_df)
    tf_true <- total_fisher_info(case_obj$theta_true, counts_now,
                                 case_obj$S_ab, case_obj$c_ab, lib$ab_df,
                                 eps = cfg$fisher_eps)
    raw_proxy <- proxy_risk(case_obj$G_true, tf_true$J_total, ridge = cfg$ridge_selection)
    n_eff <- t * cfg$n_trials
    scaled_proxy[t] <- n_eff * raw_proxy
    scaled_bures_sq[t] <- n_eff * bures_sq

    if (isTRUE(cfg$save_full_trajectory)) {
      action_idx[t] <- a_t
      action_label[t] <- case_obj$lib$setting_labels[a_t]
      theta_path[t, ] <- theta_hat
      solver_used[t] <- as.character(fit_out$diag$solver_used %||% NA_character_)
      cvx_status[t] <- as.character(fit_out$diag$cvx_status %||% NA_character_)
      pgd_status[t] <- as.character(fit_out$diag$pgd_status %||% NA_character_)
      cvx_train_loss[t] <- as.numeric(fit_out$diag$cvx_train_loss %||% NA_real_)
      pgd_train_loss[t] <- as.numeric(fit_out$diag$pgd_train_loss %||% NA_real_)
      selected_train_loss[t] <- as.numeric(fit_out$diag$selected_train_loss %||% NA_real_)
    }
  }

  out_obj <- list(proxy = scaled_proxy, bures_sq = scaled_bures_sq, rep_seed = rep_seed)
  if (isTRUE(cfg$save_full_trajectory)) {
    out_obj$trajectory <- list(
      action_idx = action_idx,
      action_label = action_label,
      theta_hat = theta_path,
      scaled_proxy = scaled_proxy,
      scaled_bures_sq = scaled_bures_sq,
      oracle_limit_step = rep(case_obj$oracle_limit, cfg$n_total),
      oracle_gi1_proxy_step = if (policy == "oracle_GI1") scaled_proxy else rep(NA_real_, cfg$n_total),
      solver_used = solver_used,
      cvx_status = cvx_status,
      pgd_status = pgd_status,
      cvx_train_loss = cvx_train_loss,
      pgd_train_loss = pgd_train_loss,
      selected_train_loss = selected_train_loss
    )
  }

  out_obj
}

run_policy <- function(case_obj, policy, cfg, trajectories_dir = NULL) {
  cat(sprintf("  policy=%s n_rep=%d ... ", policy, cfg$n_rep))
  flush.console()
  t0 <- Sys.time()

  idx <- seq_len(cfg$n_rep)
  rep_fun <- function(r) {
    seed <- stable_seed(sprintf("case11_alpha%.1f_%s_rep_%d", case_obj$alpha, policy, r), base = cfg$seed_base)
    run_one_replicate(case_obj, policy, cfg, seed)
  }

  if (.Platform$OS.type != "windows" && cfg$n_workers > 1L) {
    reps <- parallel::mclapply(idx, rep_fun, mc.cores = cfg$n_workers)
  } else {
    reps <- lapply(idx, rep_fun)
  }

  mat_proxy <- do.call(rbind, lapply(reps, function(x) x$proxy))
  mat_bures <- do.call(rbind, lapply(reps, function(x) x$bures_sq))
  mean_curve <- colMeans(mat_proxy, na.rm = TRUE)
  se_curve <- apply(mat_proxy, 2, sd, na.rm = TRUE) / sqrt(nrow(mat_proxy))
  mean_bures <- colMeans(mat_bures, na.rm = TRUE)
  se_bures <- apply(mat_bures, 2, sd, na.rm = TRUE) / sqrt(nrow(mat_bures))

  traj_file <- NA_character_
  if (isTRUE(cfg$save_full_trajectory) && !is.null(trajectories_dir)) {
    traj_file <- file.path(trajectories_dir, sprintf("alpha_%s_case11_policy_%s_trajectory.rds", gsub("\\.", "p", format(case_obj$alpha, nsmall = 1)), policy))
    traj_payload <- list(
      format = "case11_full_trajectory_v1",
      alpha = case_obj$alpha,
      policy = policy,
      n_rep = cfg$n_rep,
      n_total = cfg$n_total,
      n_trials = cfg$n_trials,
      replicate_id = idx,
      rep_seed = vapply(reps, function(x) as.integer(x$rep_seed), integer(1)),
      action_idx = do.call(rbind, lapply(reps, function(x) x$trajectory$action_idx)),
      action_label = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$action_label, nrow = 1))),
      theta_hat = lapply(reps, function(x) x$trajectory$theta_hat),
      scaled_proxy = do.call(rbind, lapply(reps, function(x) x$trajectory$scaled_proxy)),
      scaled_bures_sq = do.call(rbind, lapply(reps, function(x) x$trajectory$scaled_bures_sq)),
      oracle_limit_step = do.call(rbind, lapply(reps, function(x) x$trajectory$oracle_limit_step)),
      oracle_gi1_proxy_step = do.call(rbind, lapply(reps, function(x) x$trajectory$oracle_gi1_proxy_step)),
      solver_used = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$solver_used, nrow = 1))),
      cvx_status = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$cvx_status, nrow = 1))),
      pgd_status = do.call(rbind, lapply(reps, function(x) matrix(x$trajectory$pgd_status, nrow = 1))),
      cvx_train_loss = do.call(rbind, lapply(reps, function(x) x$trajectory$cvx_train_loss)),
      pgd_train_loss = do.call(rbind, lapply(reps, function(x) x$trajectory$pgd_train_loss)),
      selected_train_loss = do.call(rbind, lapply(reps, function(x) x$trajectory$selected_train_loss))
    )
    saveRDS(traj_payload, traj_file)
  }

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("done (%.1fs)\n", dt))

  list(
    proxy_matrix = mat_proxy,
    mean = mean_curve,
    se = se_curve,
    bures_sq_matrix = mat_bures,
    mean_bures_sq = mean_bures,
    se_bures_sq = se_bures,
    trajectory_file = traj_file
  )
}

plot_case11 <- function(alpha_res, out_file) {
  t_start <- max(1L, as.integer(alpha_res$config$plot_t_start))
  t_end <- as.integer(alpha_res$config$n_total)
  idx <- t_start:t_end
  n_seq <- idx

  u <- alpha_res$results$uniform$mean[idx]
  e <- alpha_res$results$exact$mean[idx]
  g <- alpha_res$results$GI1$mean[idx]
  og <- alpha_res$results$oracle_GI1$mean[idx]
  oracle <- alpha_res$oracle_limit

  y_all <- c(u, e, g, og, oracle)
  y_all <- y_all[is.finite(y_all)]
  y_lim <- range(y_all)
  span <- y_lim[2] - y_lim[1]
  pad <- if (span > 0) 0.05 * span else max(1e-6, abs(y_lim[1]) * 0.05 + 1e-6)
  y_lim <- c(max(0, y_lim[1] - pad), y_lim[2] + pad)

  png(out_file, width = 1100, height = 760, res = 120)
  par(mar = c(5, 5, 4, 2))

  plot(n_seq, u, type = "l", col = "gray40", lwd = 2.5, lty = 1,
       xlab = "Adaptive step (t)",
       ylab = expression(n[eff] %.% "proxy risk"),
       main = sprintf("Case 11 (2q_A_bures), alpha=%.1f, n_rep=%d",
                      alpha_res$alpha, alpha_res$config$n_rep),
       ylim = y_lim)
  lines(n_seq, e, col = "blue", lwd = 2.5, lty = 1)
  lines(n_seq, g, col = "red", lwd = 2.5, lty = 2)
  lines(n_seq, og, col = "darkgreen", lwd = 2.5, lty = 4)
  abline(h = oracle, col = "black", lwd = 2, lty = 3)
  grid(lty = "dotted", col = "gray85")

  legend("topright",
         legend = c("Permutation", "Exact", "GI1", "Oracle GI1", "Oracle limit"),
         col = c("gray40", "blue", "red", "darkgreen", "black"),
         lty = c(1, 1, 2, 4, 3),
         lwd = c(2.5, 2.5, 2.5, 2.5, 2),
         bg = "white")

  mtext(sprintf("Shown t=%d..%d | Final shown values: U=%.4f, E=%.4f, G=%.4f, OG=%.4f, Oracle=%.4f",
                t_start, t_end,
                tail(u, 1), tail(e, 1), tail(g, 1), tail(og, 1), oracle),
        side = 1, line = 3.5, cex = 0.9)

  dev.off()
}

plot_case11_relative_oracle <- function(alpha_res, out_file) {
  t_start <- max(1L, as.integer(alpha_res$config$plot_t_start))
  t_end <- as.integer(alpha_res$config$n_total)
  idx <- t_start:t_end
  t_seq <- idx

  oracle <- alpha_res$oracle_limit
  u <- alpha_res$results$uniform$mean[idx] / oracle
  e <- alpha_res$results$exact$mean[idx] / oracle
  g <- alpha_res$results$GI1$mean[idx] / oracle
  og <- alpha_res$results$oracle_GI1$mean[idx] / oracle

  y_all <- c(u, e, g, og)
  y_all <- y_all[is.finite(y_all)]
  y_max <- if (length(y_all) > 0) max(y_all) else 1.2
  pad <- if (y_max > 1) 0.06 * (y_max - 1) else 0.08
  y_lim <- c(1, y_max + max(pad, 0.05))

  png(out_file, width = 1100, height = 760, res = 120)
  par(mar = c(5, 5, 4, 2))

  plot(t_seq, u, type = "l", col = "gray40", lwd = 2.5, lty = 1,
       xlab = "Adaptive step (t)",
       ylab = expression((n[eff] %.% "proxy risk") / "Oracle limit"),
       main = sprintf("Case 11 Relative Oracle (2q_A_bures), alpha=%.1f, n_rep=%d",
                      alpha_res$alpha, alpha_res$config$n_rep),
       ylim = y_lim)
  lines(t_seq, e, col = "blue", lwd = 2.5, lty = 1)
  lines(t_seq, g, col = "red", lwd = 2.5, lty = 2)
  lines(t_seq, og, col = "darkgreen", lwd = 2.5, lty = 4)
  abline(h = 1, col = "black", lwd = 2, lty = 3)
  grid(lty = "dotted", col = "gray85")

  legend("topright",
         legend = c("Permutation", "Exact", "GI1", "Oracle GI1", "Oracle ratio = 1"),
         col = c("gray40", "blue", "red", "darkgreen", "black"),
         lty = c(1, 1, 2, 4, 3),
         lwd = c(2.5, 2.5, 2.5, 2.5, 2),
         bg = "white")

  mtext(sprintf("Shown t=%d..%d | Final shown ratios: U=%.4f, E=%.4f, G=%.4f, OG=%.4f",
                t_start, t_end, tail(u, 1), tail(e, 1), tail(g, 1), tail(og, 1)),
        side = 1, line = 3.5, cex = 0.9)
  dev.off()
}

plot_case11_bures_vs_oracle <- function(alpha_res, out_file) {
  t_start <- max(1L, as.integer(alpha_res$config$plot_t_start))
  t_end <- as.integer(alpha_res$config$n_total)
  idx <- t_start:t_end
  t_seq <- idx

  u <- alpha_res$results$uniform$mean_bures_sq[idx]
  e <- alpha_res$results$exact$mean_bures_sq[idx]
  g <- alpha_res$results$GI1$mean_bures_sq[idx]
  og <- alpha_res$results$oracle_GI1$mean_bures_sq[idx]
  oracle <- alpha_res$oracle_limit

  y_all <- c(u, e, g, og, oracle)
  y_all <- y_all[is.finite(y_all)]
  y_lim <- range(y_all)
  span <- y_lim[2] - y_lim[1]
  pad <- if (span > 0) 0.05 * span else max(1e-6, abs(y_lim[1]) * 0.05 + 1e-6)
  y_lim <- c(max(0, y_lim[1] - pad), y_lim[2] + pad)

  png(out_file, width = 1100, height = 760, res = 120)
  par(mar = c(5, 5, 4, 2))

  plot(t_seq, u, type = "l", col = "gray40", lwd = 2.5, lty = 1,
       xlab = "Adaptive step (t)",
       ylab = expression(n[eff] %.% "Bures"^2),
       main = sprintf("Case 11 Bures^2 vs Oracle (2q_A_bures), alpha=%.1f, n_rep=%d",
                      alpha_res$alpha, alpha_res$config$n_rep),
       ylim = y_lim)
  lines(t_seq, e, col = "blue", lwd = 2.5, lty = 1)
  lines(t_seq, g, col = "red", lwd = 2.5, lty = 2)
  lines(t_seq, og, col = "darkgreen", lwd = 2.5, lty = 4)
  abline(h = oracle, col = "black", lwd = 2, lty = 3)
  grid(lty = "dotted", col = "gray85")

  legend("topright",
         legend = c("Permutation", "Exact", "GI1", "Oracle GI1", "Oracle limit"),
         col = c("gray40", "blue", "red", "darkgreen", "black"),
         lty = c(1, 1, 2, 4, 3),
         lwd = c(2.5, 2.5, 2.5, 2.5, 2),
         bg = "white")

  mtext(sprintf("Shown t=%d..%d | Final shown values: U=%.4f, E=%.4f, G=%.4f, OG=%.4f, Oracle=%.4f",
                t_start, t_end, tail(u, 1), tail(e, 1), tail(g, 1), tail(og, 1), oracle),
        side = 1, line = 3.5, cex = 0.9)
  dev.off()
}

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
config_profile <- parse_cli_value("--config_profile=", args = args, default = "case11")
config_file <- parse_cli_value("--config_file=", args = args, default = NULL)
cfg <- get_case11_config(profile = config_profile, config_file = config_file, args = args)
cfg$plot_t_start <- as_int_or(parse_cli_value("--plot_t_start=", args = args, default = as.character(cfg$plot_t_start %||% 20L)), cfg$plot_t_start %||% 20L)

# Legacy alias support.
uniform_mode_cli <- parse_cli_value("--uniform_policy_mode=", args = args, default = NA_character_)
if (!is.na(uniform_mode_cli) && nzchar(uniform_mode_cli)) {
  cfg$baseline_policy_mode <- resolve_baseline_policy_mode(uniform_mode_cli)
}
cfg$uniform_policy_mode <- if (identical(cfg$baseline_policy_mode, "permutation_cyclic")) "cyclic" else "random"

if (cfg$solver_selection_mode == "pilot_choose") cfg$solver_selection_mode <- "cvx_fallback_pgd"
if (!(cfg$solver_selection_mode %in% c("cvx_fallback_pgd", "always_dual", "single_pgd"))) cfg$solver_selection_mode <- "cvx_fallback_pgd"
if (cfg$n_workers < 1L) cfg$n_workers <- 1L
if (cfg$plot_t_start < 1L) cfg$plot_t_start <- 1L

options(qt_pgd_max_iter = as.integer(cfg$pgd_max_iter))

output_dir <- parse_cli_value("--output_dir=", args = args, default = cfg$output_dir %||% "results_case11A_bures_final")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
plots_dir <- file.path(output_dir, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
trajectories_dir <- file.path(output_dir, "trajectories")
if (!dir.exists(trajectories_dir)) dir.create(trajectories_dir, recursive = TRUE)

cat("Case11 configuration:\n")
cat(sprintf("  config_profile=%s\n", config_profile))
cat(sprintf("  alpha=%.3f n_rep=%d n_total=%d n_trials=%d\n", cfg$alpha, cfg$n_rep, cfg$n_total, cfg$n_trials))
cat(sprintf("  solver_selection_mode=%s cvx_solver_preference=%s\n", cfg$solver_selection_mode, cfg$cvx_solver_preference))
cat(sprintf("  baseline_policy_mode=%s\n", cfg$baseline_policy_mode))
cat(sprintf("  eta_mle=%.2e ridge_selection=%.2e fisher_eps=%.2e nll_eps_log=%.2e bures_tol=%.2e\n",
            cfg$eta_mle, cfg$ridge_selection, cfg$fisher_eps, cfg$nll_eps_log, cfg$bures_tol))
cat(sprintf("  n_workers=%d pgd_max_iter=%d plot_t_start=%d plot_t_end=%d\n",
            cfg$n_workers, cfg$pgd_max_iter, cfg$plot_t_start, cfg$plot_t_end))
cat(sprintf("  save_full_trajectory=%s output_dir=%s\n\n", ifelse(cfg$save_full_trajectory, "TRUE", "FALSE"), output_dir))

case_obj <- build_case11(cfg$alpha, cfg)
policies <- c("uniform", "exact", "GI1", "oracle_GI1")

res <- list()
for (p in policies) {
  res[[p]] <- run_policy(case_obj, p, cfg, trajectories_dir = trajectories_dir)
}

alpha_res <- list(
  alpha = cfg$alpha,
  oracle_limit = case_obj$oracle_limit,
  uniform_limit = case_obj$uniform_limit,
  config = cfg,
  results = res
)

out_plot <- file.path(plots_dir, sprintf("alpha_%s_case11_proxy_rep%d.png",
                                         gsub("\\.", "p", format(cfg$alpha, nsmall = 1)), cfg$n_rep))
plot_case11(alpha_res, out_plot)
out_relative <- file.path(plots_dir, sprintf("alpha_%s_case11_relative_oracle_proxy_rep%d.png",
                                             gsub("\\.", "p", format(cfg$alpha, nsmall = 1)), cfg$n_rep))
plot_case11_relative_oracle(alpha_res, out_relative)
out_bures <- file.path(plots_dir, sprintf("alpha_%s_case11_bures_sq_vs_oracle_rep%d.png",
                                          gsub("\\.", "p", format(cfg$alpha, nsmall = 1)), cfg$n_rep))
plot_case11_bures_vs_oracle(alpha_res, out_bures)

summary_df <- do.call(rbind, lapply(policies, function(p) {
  data.frame(
    alpha = cfg$alpha,
    policy = p,
    t = seq_len(cfg$n_total),
    n_eff = seq_len(cfg$n_total) * cfg$n_trials,
    mean_scaled_proxy = res[[p]]$mean,
    se_scaled_proxy = res[[p]]$se,
    mean_scaled_bures_sq = res[[p]]$mean_bures_sq,
    se_scaled_bures_sq = res[[p]]$se_bures_sq,
    stringsAsFactors = FALSE
  )
}))

traj_manifest <- data.frame(
  alpha = cfg$alpha,
  policy = policies,
  path = vapply(policies, function(p) as.character(res[[p]]$trajectory_file %||% NA_character_), character(1)),
  stringsAsFactors = FALSE
)

write.csv(summary_df, file.path(output_dir, "case11_rep1000_proxy_curves.csv"), row.names = FALSE)
write.csv(traj_manifest, file.path(output_dir, "case11_trajectory_manifest.csv"), row.names = FALSE)
saveRDS(alpha_res, file.path(output_dir, "case11_rep1000_full_results.rds"))
saveRDS(cfg, file.path(output_dir, "case11_rep1000_config.rds"))

final_rows <- data.frame(
  alpha = cfg$alpha,
  n_rep = cfg$n_rep,
  n_total = cfg$n_total,
  n_trials = cfg$n_trials,
  uniform_final = tail(res$uniform$mean, 1),
  exact_final = tail(res$exact$mean, 1),
  gi1_final = tail(res$GI1$mean, 1),
  oracle_gi1_final = tail(res$oracle_GI1$mean, 1),
  oracle_limit = case_obj$oracle_limit,
  exact_minus_oracle = tail(res$exact$mean, 1) - case_obj$oracle_limit,
  gi1_minus_oracle = tail(res$GI1$mean, 1) - case_obj$oracle_limit,
  oracle_gi1_minus_oracle = tail(res$oracle_GI1$mean, 1) - case_obj$oracle_limit,
  uniform_minus_oracle = tail(res$uniform$mean, 1) - case_obj$oracle_limit,
  stringsAsFactors = FALSE
)
write.csv(final_rows, file.path(output_dir, "case11_rep1000_final_summary.csv"), row.names = FALSE)

cat("Saved outputs:\n")
cat(sprintf("  - %s\n", file.path(output_dir, "case11_rep1000_proxy_curves.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "case11_trajectory_manifest.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "case11_rep1000_final_summary.csv")))
cat(sprintf("  - %s\n", out_plot))
cat(sprintf("  - %s\n", out_relative))
cat(sprintf("  - %s\n", out_bures))
