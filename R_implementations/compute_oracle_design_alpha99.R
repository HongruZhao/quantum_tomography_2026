suppressPackageStartupMessages({
  library(CVXR)
  library(Matrix)
})

setwd("/Users/arthur/myfile/Research simulation/Quantum_Xiaxuan/R_implementations")

source("01_utilities.R")
source("02_state_basis.R")
source("03_measurement_library_1q.R")
source("04_measurement_library_2q.R")
source("05_mle_cvxr_solver.R")
source("06_fisher_and_metrics.R")

alpha <- 0.99

rho_pure_1q_default <- fixed_pure_state_1q_plus()
rho_pure_1q_L1 <- fixed_pure_state_1q_yphase(phi = 0.2)
rho_mix_1q <- fixed_pure_state_1q_0()
rho_pure_2q_phi <- fixed_pure_state_2q_phi_phase(phi = 0.2)

experiment_grid <- list(
  list(system = "1q", lib = "L1", loss = "frobenius",  case_num = 1,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L1", loss = "bures",      case_num = 2,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L1", loss = "observable", case_num = 3,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L2", loss = "frobenius",  case_num = 4,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L2", loss = "bures",      case_num = 5,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L2", loss = "observable", case_num = 6,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L3", loss = "frobenius",  case_num = 7,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L3", loss = "bures",      case_num = 8,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "1q", lib = "L3", loss = "observable", case_num = 9,  pure_state = "cos(0.2)|0> + i sin(0.2)|1>"),
  list(system = "2q", lib = "A",  loss = "frobenius",  case_num = 10, pure_state = "cos(0.2)|00> + i sin(0.2)|11>"),
  list(system = "2q", lib = "A",  loss = "bures",      case_num = 11, pure_state = "cos(0.2)|00> + i sin(0.2)|11>"),
  list(system = "2q", lib = "A",  loss = "observable", case_num = 12, pure_state = "cos(0.2)|00> + i sin(0.2)|11>"),
  list(system = "2q", lib = "B",  loss = "frobenius",  case_num = 13, pure_state = "cos(0.2)|00> + i sin(0.2)|11>"),
  list(system = "2q", lib = "B",  loss = "bures",      case_num = 14, pure_state = "cos(0.2)|00> + i sin(0.2)|11>"),
  list(system = "2q", lib = "B",  loss = "observable", case_num = 15, pure_state = "cos(0.2)|00> + i sin(0.2)|11>")
)

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

solve_optimal_design <- function(I_list, G, solver = "SCS") {
  k <- length(I_list)
  d <- nrow(G)

  pi <- Variable(k)
  W <- Variable(d, d, symmetric = TRUE)

  I_expr <- 0
  for (i in seq_len(k)) {
    I_expr <- I_expr + pi[i] * I_list[[i]]
  }
  I_expr <- (I_expr + t(I_expr)) / 2

  I_block <- bmat(list(
    list(I_expr, diag(d)),
    list(diag(d), W)
  ))

  obj <- Minimize(sum_entries(G * W))
  constraints <- list(sum(pi) == 1, pi >= 0, I_block %>>% 0)

  prob <- Problem(obj, constraints)
  result <- solve(
    prob,
    solver = solver,
    verbose = FALSE,
    max_iters = 50000,
    eps = 1e-6
  )

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

rows <- list()

rho_true_1q_default <- mix_with_state(rho_pure_1q_default, rho_mix_1q, alpha)
rho_true_1q_L1 <- mix_with_state(rho_pure_1q_L1, rho_mix_1q, alpha)
rho_true_2q_phi <- mix_with_maximally_mixed(rho_pure_2q_phi, alpha)

for (exp in experiment_grid) {
  lib <- get_library(exp$lib)
  basis <- get_basis(exp$system)
  sigmas <- basis$sigmas
  N <- lib$N

  if (exp$system == "1q") {
    rho_true <- rho_true_1q_L1
  } else {
    rho_true <- rho_true_2q_phi
  }
  theta_true <- theta_from_rho(rho_true, sigmas)

  metric_fun <- get_metric_function(exp$loss, sigmas, N, lib)
  G_true <- metric_fun(theta_true)

  sc <- build_Sab_cab(sigmas, lib$Q_list, N, lib$ab_df)
  I_list <- fisher_info_all_settings(theta_true, sc$S_ab, sc$c_ab, lib$ab_df)

  k <- length(I_list)
  I_uniform <- Reduce("+", I_list) / k
  uniform_val <- proxy_risk(G_true, I_uniform)

  opt <- solve_optimal_design(I_list, G_true)
  if (opt$status != "optimal" && opt$status != "optimal_inaccurate") {
    warning(sprintf("Case %d (alpha=%.2f) status: %s", exp$case_num, alpha, opt$status))
  }
  pi_opt <- as.numeric(opt$pi)
  I_opt <- Reduce("+", lapply(seq_len(k), function(i) pi_opt[i] * I_list[[i]]))
  opt_val <- proxy_risk(G_true, I_opt)

  if (!is.finite(opt_val) || opt_val <= 0 || opt_val > uniform_val * (1 + 1e-8)) {
    opt_md <- solve_optimal_design_mirror(I_list, G_true)
    pi_opt <- opt_md$pi
    I_opt <- Reduce("+", lapply(seq_len(k), function(i) pi_opt[i] * I_list[[i]]))
    opt_val <- proxy_risk(G_true, I_opt)
  }

  rows[[length(rows) + 1]] <- data.frame(
    alpha = alpha,
    case = exp$case_num,
    system = exp$system,
    library = exp$lib,
    loss = exp$loss,
    pure_state = exp$pure_state,
    k = k,
    oracle_min = opt_val,
    uniform = uniform_val,
    ratio_uniform_over_opt = uniform_val / opt_val
  )
}

out <- do.call(rbind, rows)
write.csv(out, "results/oracle_design_comparison_alpha_099.csv", row.names = FALSE)
print(out)
cat("\nSaved: results/oracle_design_comparison_alpha_099.csv\n")
