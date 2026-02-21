## ============================================================================
## 05_mle_cvxr_solver.R
## Stabilized MLE solver using CVXR with eigenvalue floor constraint
## Based on: 05_mle_cvxr_solver.md
## ============================================================================

source("04_measurement_library_2q.R")

## ==========================================================================
## AFFINE PROBABILITY MODEL: p_{a,b}(θ) = c_{a,b} + (1/2) s(a,b)^T θ
## ==========================================================================

#' Project a vector onto the simplex {x >= 0, sum(x) = s}
#' @param v Numeric vector
#' @param s Target sum (> 0)
#' @return Projected vector
project_simplex <- function(v, s = 1) {
  v <- as.numeric(v)
  n <- length(v)
  if (n == 0) return(v)
  if (s <= 0) return(rep(0, n))

  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u - (cssv - s) / seq_along(u) > 0))
  theta <- (cssv[rho] - s) / rho
  pmax(v - theta, 0)
}

#' Enforce rho ⪰ eta I with trace 1 by spectral projection
#' @param rho Complex Hermitian density estimate
#' @param eta Eigenvalue floor
#' @param sigmas Basis matrices for theta extraction
#' @return List with rho_hat, theta_hat, min_eig
enforce_rho_floor <- function(rho, eta, sigmas) {
  N <- nrow(rho)
  H <- hermitianize(rho)
  ee <- eigen(H, symmetric = FALSE)
  lam <- Re(ee$values)
  U <- ee$vectors

  # Project eigenvalues onto {x_i >= eta, sum x_i = 1}
  slack <- 1 - N * eta
  if (slack <= 0) {
    lam_proj <- rep(1 / N, N)
  } else {
    z <- project_simplex(lam - eta, s = slack)
    lam_proj <- z + eta
  }

  rho_proj <- U %*% diag(lam_proj, nrow = N) %*% Conj(t(U))
  rho_proj <- hermitianize(rho_proj)
  # Normalize trace defensively.
  rho_proj <- rho_proj / Re(traceC(rho_proj))

  min_eig <- min(Re(eigen(rho_proj, symmetric = FALSE, only.values = TRUE)$values))
  theta_proj <- theta_from_rho(rho_proj, sigmas)

  list(rho_hat = rho_proj, theta_hat = theta_proj, min_eig = min_eig)
}

#' Strictly enforce rho >= eta I, with final identity mixing if needed
#' @param rho Density matrix candidate
#' @param eta Eigenvalue floor
#' @param sigmas Basis list for theta conversion
#' @param tol Numerical tolerance
#' @return List(rho_hat, theta_hat, min_eig)
enforce_rho_floor_strict <- function(rho, eta, sigmas, tol = 0) {
  proj <- enforce_rho_floor(rho, eta, sigmas)
  rho_hat <- hermitianize(proj$rho_hat)
  N <- nrow(rho_hat)
  lam_min <- min(Re(eigen(rho_hat, symmetric = FALSE, only.values = TRUE)$values))

  if (!is.finite(lam_min)) {
    rho_hat <- diag(N) / N
    lam_min <- 1 / N
  }

  if (lam_min < eta - tol) {
    denom <- 1 / N - lam_min
    w <- if (denom <= 0) 1 else (eta - lam_min) / denom
    w <- min(1, max(0, w))
    rho_hat <- hermitianize((1 - w) * rho_hat + w * diag(N) / N)
    lam_min <- min(Re(eigen(rho_hat, symmetric = FALSE, only.values = TRUE)$values))
  }

  list(
    rho_hat = rho_hat,
    theta_hat = theta_from_rho(rho_hat, sigmas),
    min_eig = lam_min
  )
}

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
## MULTI-SOLVER STATS (MOSEK / ECOS / SCS / PGD)
## ==========================================================================

.mle_solver_names <- c("MOSEK", "ECOS", "SCS", "PGD")
.mle_solver_stats <- new.env(parent = emptyenv())

init_mle_solver_stats <- function() {
  assign("attempted", setNames(integer(length(.mle_solver_names)), .mle_solver_names), envir = .mle_solver_stats)
  assign("optimal", setNames(integer(length(.mle_solver_names)), .mle_solver_names), envir = .mle_solver_stats)
  assign("selected", setNames(integer(length(.mle_solver_names)), .mle_solver_names), envir = .mle_solver_stats)
  assign("selected_optimal", setNames(integer(length(.mle_solver_names)), .mle_solver_names), envir = .mle_solver_stats)
}

if (!exists("attempted", envir = .mle_solver_stats, inherits = FALSE)) {
  init_mle_solver_stats()
}

reset_mle_solver_stats <- function() {
  init_mle_solver_stats()
  invisible(NULL)
}

record_mle_solver_stat <- function(field, solver_name, increment = 1L) {
  if (!(solver_name %in% .mle_solver_names)) return(invisible(NULL))
  vec <- get(field, envir = .mle_solver_stats, inherits = FALSE)
  vec[solver_name] <- vec[solver_name] + as.integer(increment)
  assign(field, vec, envir = .mle_solver_stats)
  invisible(NULL)
}

get_mle_solver_stats <- function() {
  data.frame(
    solver = .mle_solver_names,
    attempted = as.integer(get("attempted", envir = .mle_solver_stats, inherits = FALSE)[.mle_solver_names]),
    optimal = as.integer(get("optimal", envir = .mle_solver_stats, inherits = FALSE)[.mle_solver_names]),
    selected = as.integer(get("selected", envir = .mle_solver_stats, inherits = FALSE)[.mle_solver_names]),
    selected_optimal = as.integer(get("selected_optimal", envir = .mle_solver_stats, inherits = FALSE)[.mle_solver_names]),
    stringsAsFactors = FALSE
  )
}

## ==========================================================================
## PROJECTED-GRADIENT MLE (NON-CVXR FALLBACK)
## ==========================================================================

fit_theta_pgd <- function(N, sigmas, S_ab, c_ab, Nab,
                          eta = 0,
                          eps_log = 0,
                          max_iter = 1200L,
                          step0 = 1.0,
                          tol = 1e-8,
                          theta_start = NULL,
                          truncate_density = TRUE) {
  d <- ncol(S_ab)
  theta <- if (!is.null(theta_start) && length(theta_start) == d &&
               all(is.finite(theta_start))) {
    as.numeric(theta_start)
  } else {
    rep(0, d)
  }
  f <- nll_theta(theta, S_ab, c_ab, Nab, eps = eps_log)

  if (!is.finite(f)) {
    theta <- theta_from_rho(diag(N) / N, sigmas)
    f <- nll_theta(theta, S_ab, c_ab, Nab, eps = eps_log)
  }

  for (iter in seq_len(max_iter)) {
    p <- probs_from_theta(theta, S_ab, c_ab)
    if (any(!is.finite(p)) || any(p <= 0)) break
    grad <- as.numeric(-0.5 * crossprod(S_ab, Nab / p))
    g2 <- sum(grad^2)
    if (!is.finite(g2) || g2 < 1e-20) break

    step <- step0
    accepted <- FALSE
    for (bt in 1:40) {
      theta_try <- theta - step * grad
      if (truncate_density) {
        rho_try <- rho_of_theta(theta_try, sigmas, N)
        proj <- enforce_rho_floor_strict(rho_try, eta, sigmas)
        theta_next <- proj$theta_hat
      } else {
        theta_next <- theta_try
      }
      p_proj <- probs_from_theta(theta_next, S_ab, c_ab)
      if (any(!is.finite(p_proj)) || any(p_proj <= 0)) {
        step <- step * 0.5
        next
      }
      f_new <- nll_theta(theta_next, S_ab, c_ab, Nab, eps = eps_log)
      if (is.finite(f_new) && f_new <= f - 1e-4 * step * g2) {
        theta <- theta_next
        if (abs(f - f_new) / max(1, abs(f)) < tol) {
          f <- f_new
          accepted <- TRUE
          break
        }
        f <- f_new
        accepted <- TRUE
        break
      }
      step <- step * 0.5
    }
    if (!accepted) break
  }

  rho_hat <- rho_of_theta(theta, sigmas, N)
  min_eig_raw <- min(Re(eigen(hermitianize(rho_hat), symmetric = FALSE, only.values = TRUE)$values))
  floor_enforced <- FALSE
  min_eig_hat <- min_eig_raw
  if (truncate_density && (!is.finite(min_eig_raw) || min_eig_raw < eta)) {
    proj <- enforce_rho_floor_strict(rho_hat, eta, sigmas)
    rho_hat <- proj$rho_hat
    theta <- proj$theta_hat
    min_eig_hat <- proj$min_eig
    floor_enforced <- TRUE
  }

  list(
    theta_hat = theta,
    rho_hat = rho_hat,
    status = "optimal",
    value = f,
    train_loss = f,
    floor_enforced = floor_enforced,
    min_eig_raw = min_eig_raw,
    min_eig_hat = min_eig_hat,
    solver_used = "PGD"
  )
}

## ==========================================================================
## STABILIZED MLE WITH EIGENVALUE FLOOR (CVXR)
## ==========================================================================

#' Fit theta via MLE with explicit eigenvalue-floor control
#' Solves: min_θ -Σ_{a,b} N_{a,b} log(p_{a,b}(θ))
#'         s.t. ρ(θ) ⪰ η I
#' @param N Hilbert space dimension
#' @param sigmas List of basis matrices
#' @param S_ab Coefficient matrix (M x d)
#' @param c_ab Constant vector (length M)
#' @param Nab Count vector (length M)
#' @param eta Eigenvalue floor (default 0; requires eta < 1/N)
#' @param solver CVXR/optimizer selector ("AUTO", "AUTO_PGD", "MOSEK", "ECOS", "SCS", "PGD")
#' @param eps_log Optional additive offset inside log-likelihood (default 0 = exact likelihood)
#' @param verbose Print solver output
#' @return List with theta_hat, rho_hat, status, value
fit_theta_cvxr <- function(N, sigmas, S_ab, c_ab, Nab,
                           eta = 0,
                           solver = c("AUTO", "AUTO_PGD", "MOSEK", "ECOS", "SCS", "PGD"),
                           eps_log = 0,
                           verbose = FALSE,
                           theta_start = NULL,
                           truncate_density = TRUE,
                           allow_eta_zero = TRUE) {

  solver <- match.arg(solver)
  cvxr_available <- if (exists("CVXR_AVAILABLE", inherits = TRUE)) {
    isTRUE(get("CVXR_AVAILABLE", inherits = TRUE))
  } else {
    requireNamespace("CVXR", quietly = TRUE)
  }

  if (!cvxr_available && solver != "PGD") {
    warning("CVXR not available; switching to solver = 'PGD'.")
    solver <- "PGD"
  }

  pgd_max_iter <- suppressWarnings(as.integer(getOption("qt_pgd_max_iter", 1200L)))
  if (is.na(pgd_max_iter) || pgd_max_iter < 10L) pgd_max_iter <- 1200L

  d <- length(sigmas)
  eta_min <- if (allow_eta_zero) 0 else 1e-3

  # Explicit eigenvalue floor control.
  if (!is.finite(eta) || eta < eta_min) {
    warning(sprintf("eta must be >= %.1e; using eta = %.1e", eta_min, eta_min))
    eta <- eta_min
  }
  if (eta >= 1/N) {
    stop(sprintf("eta = %.4g must satisfy eta < 1/N = %.4g", eta, 1/N))
  }

  if (solver == "PGD") {
    record_mle_solver_stat("attempted", "PGD")
    fit_pgd <- fit_theta_pgd(
      N = N, sigmas = sigmas, S_ab = S_ab, c_ab = c_ab, Nab = Nab,
      eta = eta, eps_log = eps_log, max_iter = pgd_max_iter,
      theta_start = theta_start,
      truncate_density = truncate_density
    )
    record_mle_solver_stat("optimal", "PGD")
    record_mle_solver_stat("selected", "PGD")
    record_mle_solver_stat("selected_optimal", "PGD")
    return(list(
      theta_hat = fit_pgd$theta_hat,
      rho_hat = fit_pgd$rho_hat,
      status = fit_pgd$status,
      value = fit_pgd$value,
      floor_enforced = fit_pgd$floor_enforced,
      min_eig_raw = fit_pgd$min_eig_raw,
      min_eig_hat = fit_pgd$min_eig_hat,
      solver_used = "PGD",
      solver_mode = solver,
      solver_status_map = data.frame(
        solver = .mle_solver_names,
        status = c("not_requested", "not_requested", "not_requested", "optimal"),
        cvx_value = c(NA_real_, NA_real_, NA_real_, fit_pgd$value),
        train_loss = c(NA_real_, NA_real_, NA_real_, fit_pgd$train_loss),
        floor_enforced = c(FALSE, FALSE, FALSE, fit_pgd$floor_enforced),
        min_eig_raw = c(NA_real_, NA_real_, NA_real_, fit_pgd$min_eig_raw),
        min_eig_hat = c(NA_real_, NA_real_, NA_real_, fit_pgd$min_eig_hat),
        stringsAsFactors = FALSE
      )
    ))
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
  # No data: return maximally mixed estimate.
  if (sum(Nab) <= 0) {
    theta_hat <- rep(0, d)
    rho_hat <- rho_of_theta(theta_hat, sigmas, N)
    return(list(
      theta_hat = theta_hat,
      rho_hat = rho_hat,
      status = "trivial",
      value = 0,
      floor_enforced = FALSE,
      min_eig_raw = min(Re(eigen(hermitianize(rho_hat), symmetric = FALSE, only.values = TRUE)$values)),
      min_eig_hat = min(Re(eigen(hermitianize(rho_hat), symmetric = FALSE, only.values = TRUE)$values)),
      solver_used = NA_character_,
      solver_mode = solver,
      solver_status_map = data.frame(
        solver = .mle_solver_names,
        status = "not_requested",
        cvx_value = NA_real_,
        train_loss = NA_real_,
        floor_enforced = FALSE,
        min_eig_raw = NA_real_,
        min_eig_hat = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  # Exact objective by default; optional epsilon shift can be enabled explicitly.
  obj <- if (eps_log > 0) {
    -sum_entries(Nab * log(p_expr + eps_log))
  } else {
    -sum_entries(Nab * log(p_expr))
  }

  # Constraints
  rho_floor <- Variable(2 * N, 2 * N, PSD = TRUE)
  constraints <- list(
    rho_var == A_affine,  # rho PSD via variable
    p_expr >= 0,          # Probabilities non-negative
    rho_floor == A_affine - eta * SI  # ρ(θ) ⪰ ηI
  )

  # Define and solve problem
  prob <- Problem(Minimize(obj), constraints)

  is_solver_available <- function(solver_name) {
    switch(solver_name,
      "MOSEK" = requireNamespace("Rmosek", quietly = TRUE),
      "ECOS" = requireNamespace("ECOSolveR", quietly = TRUE),
      "SCS" = TRUE,
      "PGD" = TRUE,
      FALSE
    )
  }

  solver_candidates <- if (solver == "AUTO_PGD") {
    c("MOSEK", "ECOS", "SCS", "PGD")
  } else if (solver == "AUTO") {
    c("MOSEK", "ECOS", "SCS")
  } else {
    solver
  }
  solver_records <- list()
  best <- NULL

  for (solver_try in solver_candidates) {
    record_mle_solver_stat("attempted", solver_try)

    if (!is_solver_available(solver_try)) {
      solver_records[[solver_try]] <- list(
        status = "unavailable",
        cvx_value = NA_real_,
        train_loss = NA_real_,
        floor_enforced = FALSE,
        min_eig_raw = NA_real_,
        min_eig_hat = NA_real_
      )
      next
    }

    if (solver_try == "PGD") {
      fit_pgd <- fit_theta_pgd(
        N = N, sigmas = sigmas, S_ab = S_ab, c_ab = c_ab, Nab = Nab,
        eta = eta, eps_log = eps_log, max_iter = pgd_max_iter,
        theta_start = theta_start,
        truncate_density = truncate_density
      )
      status_try <- fit_pgd$status
      value_try <- fit_pgd$value
      train_loss_try <- fit_pgd$train_loss
      floor_enforced_try <- fit_pgd$floor_enforced
      min_eig_raw_try <- fit_pgd$min_eig_raw
      min_eig_hat_try <- fit_pgd$min_eig_hat
      theta_try <- fit_pgd$theta_hat
      rho_try <- fit_pgd$rho_hat

      if (status_try %in% c("optimal", "optimal_inaccurate")) {
        record_mle_solver_stat("optimal", solver_try)
      }

      solver_records[[solver_try]] <- list(
        status = status_try,
        cvx_value = value_try,
        train_loss = train_loss_try,
        floor_enforced = floor_enforced_try,
        min_eig_raw = min_eig_raw_try,
        min_eig_hat = min_eig_hat_try
      )

      is_candidate <- status_try %in% c("optimal", "optimal_inaccurate") && is.finite(train_loss_try)
      if (is_candidate) {
        candidate <- list(
          solver = solver_try,
          status = status_try,
          value = value_try,
          train_loss = train_loss_try,
          theta_hat = theta_try,
          rho_hat = rho_try,
          floor_enforced = floor_enforced_try,
          min_eig_raw = min_eig_raw_try,
          min_eig_hat = min_eig_hat_try
        )

        choose_candidate <- FALSE
        if (is.null(best)) {
          choose_candidate <- TRUE
        } else if (!is.finite(best$train_loss) && is.finite(candidate$train_loss)) {
          choose_candidate <- TRUE
        } else if (is.finite(best$train_loss) && is.finite(candidate$train_loss) &&
                   candidate$train_loss < best$train_loss - 1e-9) {
          choose_candidate <- TRUE
        } else if (is.finite(best$train_loss) && is.finite(candidate$train_loss) &&
                   abs(candidate$train_loss - best$train_loss) <= 1e-9 &&
                   candidate$status == "optimal" && best$status != "optimal") {
          choose_candidate <- TRUE
        }

        if (choose_candidate) best <- candidate
      }
      next
    }

    res <- tryCatch({
      solve(prob, solver = solver_try, verbose = verbose, max_iters = 5000)
    }, error = function(e) {
      list(status = "error", value = NA_real_, error_message = e$message)
    })

    status_try <- if (!is.null(res$status)) as.character(res$status) else "error"
    value_try <- if (!is.null(res$value) && is.finite(res$value)) as.numeric(res$value) else NA_real_
    theta_try <- rep(NA_real_, d)
    rho_try <- matrix(NA_complex_, N, N)
    floor_enforced_try <- FALSE
    min_eig_raw_try <- NA_real_
    min_eig_hat_try <- NA_real_
    train_loss_try <- NA_real_

    if (status_try %in% c("optimal", "optimal_inaccurate")) {
      record_mle_solver_stat("optimal", solver_try)
      theta_try <- tryCatch(drop(res$getValue(theta_var)), error = function(e) rep(NA_real_, d))

      if (length(theta_try) == d && all(is.finite(theta_try))) {
        rho_try <- rho_of_theta(theta_try, sigmas, N)
        min_eig_raw_try <- min(Re(eigen(hermitianize(rho_try), symmetric = FALSE, only.values = TRUE)$values))
        if (truncate_density && (!is.finite(min_eig_raw_try) || min_eig_raw_try < eta)) {
          proj <- enforce_rho_floor_strict(rho_try, eta, sigmas)
          rho_try <- proj$rho_hat
          theta_try <- proj$theta_hat
          min_eig_hat_try <- proj$min_eig
          floor_enforced_try <- TRUE
        } else {
          min_eig_hat_try <- min_eig_raw_try
        }
        train_loss_try <- nll_theta(theta_try, S_ab, c_ab, Nab, eps = eps_log)
        if (!is.finite(train_loss_try)) train_loss_try <- value_try
      } else {
        status_try <- "failed_extract"
      }
    }

    solver_records[[solver_try]] <- list(
      status = status_try,
      cvx_value = value_try,
      train_loss = train_loss_try,
      floor_enforced = floor_enforced_try,
      min_eig_raw = min_eig_raw_try,
      min_eig_hat = min_eig_hat_try
    )

    is_candidate <- status_try %in% c("optimal", "optimal_inaccurate") && is.finite(train_loss_try)
    if (is_candidate) {
      candidate <- list(
        solver = solver_try,
        status = status_try,
        value = value_try,
        train_loss = train_loss_try,
        theta_hat = theta_try,
        rho_hat = rho_try,
        floor_enforced = floor_enforced_try,
        min_eig_raw = min_eig_raw_try,
        min_eig_hat = min_eig_hat_try
      )

      choose_candidate <- FALSE
      if (is.null(best)) {
        choose_candidate <- TRUE
      } else if (!is.finite(best$train_loss) && is.finite(candidate$train_loss)) {
        choose_candidate <- TRUE
      } else if (is.finite(best$train_loss) && is.finite(candidate$train_loss) &&
                 candidate$train_loss < best$train_loss - 1e-9) {
        choose_candidate <- TRUE
      } else if (is.finite(best$train_loss) && is.finite(candidate$train_loss) &&
                 abs(candidate$train_loss - best$train_loss) <= 1e-9 &&
                 candidate$status == "optimal" && best$status != "optimal") {
        choose_candidate <- TRUE
      }

      if (choose_candidate) best <- candidate
    }
  }

  empty_record <- list(
    status = "not_requested",
    cvx_value = NA_real_,
    train_loss = NA_real_,
    floor_enforced = FALSE,
    min_eig_raw = NA_real_,
    min_eig_hat = NA_real_
  )
  solver_status_map <- do.call(rbind, lapply(.mle_solver_names, function(s) {
    rec <- solver_records[[s]]
    if (is.null(rec)) rec <- empty_record
    data.frame(
      solver = s,
      status = as.character(rec$status),
      cvx_value = as.numeric(rec$cvx_value),
      train_loss = as.numeric(rec$train_loss),
      floor_enforced = as.logical(rec$floor_enforced),
      min_eig_raw = as.numeric(rec$min_eig_raw),
      min_eig_hat = as.numeric(rec$min_eig_hat),
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(best)) {
    record_mle_solver_stat("selected", best$solver)
    if (best$status %in% c("optimal", "optimal_inaccurate")) {
      record_mle_solver_stat("selected_optimal", best$solver)
    }
    return(list(
      theta_hat = best$theta_hat,
      rho_hat = best$rho_hat,
      status = best$status,
      value = best$value,
      floor_enforced = best$floor_enforced,
      min_eig_raw = best$min_eig_raw,
      min_eig_hat = best$min_eig_hat,
      solver_used = best$solver,
      solver_mode = solver,
      solver_status_map = solver_status_map
    ))
  }

  list(
    theta_hat = rep(NA_real_, d),
    rho_hat = matrix(NA_complex_, N, N),
    status = "all_solvers_failed",
    value = NA_real_,
    floor_enforced = FALSE,
    min_eig_raw = NA_real_,
    min_eig_hat = NA_real_,
    solver_used = NA_character_,
    solver_mode = solver,
    solver_status_map = solver_status_map
  )
}

solver_metric_from_map <- function(fit, solver_name, field, fallback = NA_real_) {
  if (is.null(fit$solver_status_map) || !is.data.frame(fit$solver_status_map)) {
    return(fallback)
  }
  sm <- fit$solver_status_map
  idx <- which(sm$solver == solver_name)
  if (length(idx) == 0L || !(field %in% names(sm))) return(fallback)
  val <- suppressWarnings(as.numeric(sm[[field]][idx[1]]))
  if (length(val) == 0L || is.na(val)) fallback else val
}

solver_status_from_map <- function(fit, solver_name, fallback = NA_character_) {
  if (is.null(fit$solver_status_map) || !is.data.frame(fit$solver_status_map)) {
    return(fallback)
  }
  sm <- fit$solver_status_map
  idx <- which(sm$solver == solver_name)
  if (length(idx) == 0L || !("status" %in% names(sm))) return(fallback)
  as.character(sm$status[idx[1]])
}

is_mle_success_status <- function(status) {
  identical(status, "optimal") || identical(status, "optimal_inaccurate") || identical(status, "trivial")
}

empty_solver_fit <- function(d, N, status = "error") {
  list(
    theta_hat = rep(NA_real_, d),
    rho_hat = matrix(NA_complex_, N, N),
    status = status,
    value = NA_real_,
    floor_enforced = FALSE,
    min_eig_raw = NA_real_,
    min_eig_hat = NA_real_,
    solver_used = NA_character_,
    solver_mode = NA_character_,
    solver_status_map = NULL
  )
}

#' Fit MLE by evaluating CVX and PGD on the same data, then selecting lower train loss.
#' If CVX fails, PGD is used.
#' @param cvx_solver CVX backend (typically "SCS", "ECOS", or "MOSEK")
#' @return Same structure as fit_theta_cvxr plus dual-solver diagnostics
fit_theta_dual_compare <- function(N, sigmas, S_ab, c_ab, Nab,
                                   eta = 0,
                                   cvx_solver = "SCS",
                                   eps_log = 0,
                                   verbose = FALSE,
                                   theta_start = NULL,
                                   truncate_density = TRUE,
                                   allow_eta_zero = TRUE) {
  d <- length(sigmas)

  t0 <- proc.time()[["elapsed"]]
  fit_cvx <- tryCatch(
    fit_theta_cvxr(
      N = N, sigmas = sigmas, S_ab = S_ab, c_ab = c_ab, Nab = Nab,
      eta = eta, solver = cvx_solver, eps_log = eps_log, verbose = verbose,
      theta_start = theta_start, truncate_density = truncate_density,
      allow_eta_zero = allow_eta_zero
    ),
    error = function(e) {
      out <- empty_solver_fit(d, N, status = "error")
      out$error_message <- conditionMessage(e)
      out
    }
  )
  cvx_runtime_sec <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  fit_pgd <- tryCatch(
    fit_theta_cvxr(
      N = N, sigmas = sigmas, S_ab = S_ab, c_ab = c_ab, Nab = Nab,
      eta = eta, solver = "PGD", eps_log = eps_log, verbose = verbose,
      theta_start = theta_start, truncate_density = truncate_density,
      allow_eta_zero = allow_eta_zero
    ),
    error = function(e) {
      out <- empty_solver_fit(d, N, status = "error")
      out$error_message <- conditionMessage(e)
      out
    }
  )
  pgd_runtime_sec <- proc.time()[["elapsed"]] - t1

  cvx_status <- if (!is.null(fit_cvx$status)) as.character(fit_cvx$status) else NA_character_
  pgd_status <- if (!is.null(fit_pgd$status)) as.character(fit_pgd$status) else NA_character_
  cvx_train_loss <- solver_metric_from_map(fit_cvx, cvx_solver, "train_loss", fallback = as.numeric(fit_cvx$value))
  pgd_train_loss <- solver_metric_from_map(fit_pgd, "PGD", "train_loss", fallback = as.numeric(fit_pgd$value))

  cvx_ok <- is_mle_success_status(cvx_status) && is.finite(cvx_train_loss)
  pgd_ok <- is_mle_success_status(pgd_status) && is.finite(pgd_train_loss)

  selected <- NULL
  selected_tag <- NA_character_
  if (cvx_ok && pgd_ok) {
    if (pgd_train_loss < cvx_train_loss - 1e-9) {
      selected <- fit_pgd
      selected_tag <- "PGD"
    } else {
      selected <- fit_cvx
      selected_tag <- "CVX"
    }
  } else if (cvx_ok) {
    selected <- fit_cvx
    selected_tag <- "CVX"
  } else if (pgd_ok) {
    selected <- fit_pgd
    selected_tag <- "PGD"
  } else {
    selected <- fit_pgd
    selected_tag <- "PGD"
  }

  selected_loss <- if (selected_tag == "CVX") cvx_train_loss else pgd_train_loss

  selected$cvx_solver <- cvx_solver
  selected$cvx_status <- cvx_status
  selected$pgd_status <- pgd_status
  selected$cvx_train_loss <- cvx_train_loss
  selected$pgd_train_loss <- pgd_train_loss
  selected$selected_train_loss <- selected_loss
  selected$cvx_runtime_sec <- as.numeric(cvx_runtime_sec)
  selected$pgd_runtime_sec <- as.numeric(pgd_runtime_sec)
  selected$solver_dual_selected <- selected_tag
  selected$cvx_solver_status <- solver_status_from_map(fit_cvx, cvx_solver, fallback = cvx_status)
  selected$pgd_solver_status <- solver_status_from_map(fit_pgd, "PGD", fallback = pgd_status)
  selected$cvx_solver_used <- if (!is.null(fit_cvx$solver_used)) as.character(fit_cvx$solver_used) else NA_character_
  selected$pgd_solver_used <- if (!is.null(fit_pgd$solver_used)) as.character(fit_pgd$solver_used) else "PGD"
  selected
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
                              eta = 1e-3, solver = "SCS", verbose = FALSE) {

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
#' @param eps Optional probability floor in NLL evaluation (default 0 = disabled)
#' @return NLL value
nll_theta <- function(theta, S_ab, c_ab, Nab, eps = 0) {
  p <- probs_from_theta(theta, S_ab, c_ab)
  if (any(!is.finite(p))) return(Inf)
  if (eps > 0) {
    p <- pmax(p, eps)
  } else if (any(p <= 0)) {
    return(Inf)
  }
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
# fit <- fit_mle_from_data(L1, sigmas, a_seq, b_seq, eta = 1e-3, solver = "SCS")
# cat("Solver status:", fit$status, "\n")
# cat("Estimated theta:", round(fit$theta_hat, 4), "\n")
# cat("Theta error:", round(sqrt(sum((fit$theta_hat - theta_true)^2)), 4), "\n")

cat("05_mle_cvxr_solver.R loaded successfully\n")
