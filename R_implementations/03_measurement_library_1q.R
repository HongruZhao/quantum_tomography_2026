## ============================================================================
## 03_measurement_library_1q.R
## One-qubit measurement libraries (Libraries 1-3)
## Based on: 03_measurement_libraries_one_qubit.md
## ============================================================================

source("02_state_basis.R")

## ==========================================================================
## SHARED HELPER: Projectors for a Bloch direction
## ==========================================================================

#' Create PVM projectors for measurement along a Bloch direction
#' Q_+ = (1/2)(I + u·σ), Q_- = (1/2)(I - u·σ)
#' @param u Direction vector (length 3, will be normalized)
#' @param pauli Output from build_pauli_basis_1q() (optional)
#' @return List with Qp, Qm (projectors), u (normalized), B (observable)
pvm_projectors_for_direction_1q <- function(u, pauli = NULL) {
  if (is.null(pauli)) pauli <- build_pauli_basis_1q()

  u <- as.numeric(u)
  stopifnot(length(u) == 3)

  # Normalize to unit vector
  u <- u / sqrt(sum(u^2))

  X <- pauli$sigmas$X
  Y <- pauli$sigmas$Y
  Z <- pauli$sigmas$Z
  I2 <- pauli$I

  # Observable B = u·σ
  B <- u[1] * X + u[2] * Y + u[3] * Z

  # Projectors
  Qp <- 0.5 * (I2 + B)  # eigenvalue +1

  Qm <- 0.5 * (I2 - B)  # eigenvalue -1

  list(
    Qp = hermitianize(Qp),
    Qm = hermitianize(Qm),
    u = u,
    B = hermitianize(B)
  )
}

## ==========================================================================
## HELPER: Build ab indexing structure
## ==========================================================================

#' Build standardized measurement library indexing structure
#' Creates ab_df and ab_row for mapping (a,b) to flat index
#' @param Q_list List of POVMs, Q_list[[a]] is list of effects for setting a
#' @param setting_labels Character vector of setting names (optional)
#' @param N Hilbert space dimension (optional)
#' @param extra Additional data to store (optional)
#' @return Standardized library object
build_ab_indexing <- function(Q_list, setting_labels = NULL, N = NULL, extra = NULL) {
  k <- length(Q_list)
  r_vec <- vapply(Q_list, length, integer(1))
  total_rows <- sum(r_vec)

  # Build data frame with (a, b, row) mapping
  a_col <- rep(1:k, times = r_vec)
  b_col <- unlist(lapply(r_vec, function(r) 1:r))
  ab_df <- data.frame(a = a_col, b = b_col, row = 1:total_rows)

  # Build ab_row: ab_row[[a]][b] gives flat row index
  ab_row <- vector("list", k)
  idx <- 1
  for (a in 1:k) {
    ra <- r_vec[a]
    ab_row[[a]] <- idx:(idx + ra - 1)
    idx <- idx + ra
  }

  out <- list(
    k = k,
    r_vec = r_vec,
    Q_list = Q_list,
    ab_df = ab_df,
    ab_row = ab_row,
    M = total_rows  # total number of (a,b) cells
  )

  if (!is.null(setting_labels)) out$setting_labels <- setting_labels
  if (!is.null(N)) out$N <- N
  if (!is.null(extra)) out$extra <- extra

  out
}

## ==========================================================================
## LIBRARY 1: Pauli PVMs (k = 3)
## ==========================================================================

#' Build Library 1: Pauli PVMs (X, Y, Z measurements)
#' k1 = 3 settings, each with r = 2 outcomes
#' @return Standardized measurement library object
build_library_1q_L1 <- function() {
  pauli <- build_pauli_basis_1q()

  # Unit vectors for X, Y, Z directions
  U <- list(
    X = c(1, 0, 0),  # Measure X
    Y = c(0, 1, 0),  # Measure Y
    Z = c(0, 0, 1)   # Measure Z
  )

  P_list <- lapply(U, function(u) pvm_projectors_for_direction_1q(u, pauli))
  Q_list <- lapply(P_list, function(P) {
    list(P$Qp, P$Qm)  # outcomes: b=1 (+1), b=2 (-1)
  })

  lib <- build_ab_indexing(Q_list, setting_labels = names(U), N = 2)
  lib$library_name <- "L1_Pauli"
  lib$description <- "Pauli PVMs: measure X, Y, Z"
  lib$directions <- U
  lib$Obs <- lapply(P_list, function(P) P$B)
  lib$pauli <- pauli

  lib
}

## ==========================================================================
## LIBRARY 2: Nine-axis PVMs (k = 9)
## ==========================================================================

#' Build Library 2: Nine-axis projective PVMs
#' Directions: ex, ey, ez, ex±ey, ex±ez, ey±ez (normalized)
#' k2 = 9 settings, each with r = 2 outcomes
#' @return Standardized measurement library object
build_library_1q_L2 <- function() {
  pauli <- build_pauli_basis_1q()

  # Nine direction vectors (not all normalized yet)
  U <- list(
    ex = c(1, 0, 0),
    ey = c(0, 1, 0),
    ez = c(0, 0, 1),
    ex_py = c(1, 1, 0),   # (X+Y)/sqrt(2)
    ex_my = c(1, -1, 0),  # (X-Y)/sqrt(2)
    ex_pz = c(1, 0, 1),   # (X+Z)/sqrt(2)
    ex_mz = c(1, 0, -1),  # (X-Z)/sqrt(2)
    ey_pz = c(0, 1, 1),   # (Y+Z)/sqrt(2)
    ey_mz = c(0, 1, -1)   # (Y-Z)/sqrt(2)
  )

  # pvm_projectors_for_direction_1q normalizes automatically
  P_list <- lapply(U, function(u) pvm_projectors_for_direction_1q(u, pauli))
  Q_list <- lapply(P_list, function(P) list(P$Qp, P$Qm))

  lib <- build_ab_indexing(Q_list, setting_labels = names(U), N = 2)
  lib$library_name <- "L2_NineAxis"
  lib$description <- "Nine-axis PVMs: X, Y, Z, (X±Y)/√2, (X±Z)/√2, (Y±Z)/√2"
  lib$directions <- U
  lib$Obs <- lapply(P_list, function(P) P$B)
  lib$pauli <- pauli

  lib
}

## ==========================================================================
## LIBRARY 3: Random 4-axis PVMs (k = 4)
## ==========================================================================

#' Build Library 3: Random projective PVMs
#' Four random unit vectors on the Bloch sphere
#' k3 = 4 settings, each with r = 2 outcomes
#' @param seed Random seed for reproducibility
#' @return Standardized measurement library object
build_library_1q_L3 <- function(seed = 1) {
  # Use local RNG scope so L3 construction doesn't perturb global randomness.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  pauli <- build_pauli_basis_1q()

  # Draw 4 random unit vectors
  U <- lapply(1:4, function(j) sample_random_unit_vector(3))
  names(U) <- sprintf("u%02d", 1:4)

  P_list <- lapply(U, function(u) pvm_projectors_for_direction_1q(u, pauli))
  Q_list <- lapply(P_list, function(P) list(P$Qp, P$Qm))

  lib <- build_ab_indexing(Q_list, setting_labels = names(U), N = 2)
  lib$library_name <- "L3_Random4"
  lib$description <- sprintf("Random 4-axis PVMs (seed=%d)", seed)
  lib$directions <- U
  lib$seed <- seed
  lib$Obs <- lapply(P_list, function(P) P$B)
  lib$pauli <- pauli

  lib
}

## ==========================================================================
## VALIDATION FUNCTIONS
## ==========================================================================

#' Validate a measurement library (POVM completeness and PSD)
#' @param lib Measurement library object
#' @param tol Numerical tolerance
#' @return List with validation results
validate_library_1q <- function(lib, tol = 1e-10) {
  results <- list(
    valid = TRUE,
    messages = character(),
    setting_details = list()
  )

  N <- lib$N

  for (a in 1:lib$k) {
    Q_a <- lib$Q_list[[a]]
    ra <- length(Q_a)

    # Check POVM completeness: Σ_b Q_{a,b} = I
    Q_sum <- Reduce(`+`, Q_a)
    completeness_err <- max(Mod(Q_sum - diag(N)))

    # Check each effect is PSD
    psd_ok <- sapply(Q_a, function(Q) {
      eigs <- eigen(hermitianize(Q), only.values = TRUE)$values
      min(Re(eigs)) >= -tol
    })

    setting_valid <- (completeness_err < tol) && all(psd_ok)

    results$setting_details[[a]] <- list(
      setting = a,
      label = lib$setting_labels[a],
      completeness_error = completeness_err,
      effects_psd = all(psd_ok),
      valid = setting_valid
    )

    if (!setting_valid) {
      results$valid <- FALSE
      results$messages <- c(results$messages,
                            sprintf("Setting %d (%s) failed validation",
                                    a, lib$setting_labels[a]))
    }
  }

  results
}

## ==========================================================================
## BORN PROBABILITY FUNCTIONS
## ==========================================================================

#' Compute Born probabilities for all outcomes at a given setting
#' p_{a,b}(ρ) = tr(ρ Q_{a,b})
#' @param rho Density matrix
#' @param Q_a List of POVM effects for setting a
#' @return Probability vector (normalized)
born_probs_setting <- function(rho, Q_a) {
  p <- vapply(Q_a, function(Q) Re(traceC(rho %*% Q)), numeric(1))
  # Normalize for numerical safety
  p <- pmax(p, 0)  # Clip negative due to numerics
  p / sum(p)
}

#' Compute Born probabilities for all settings
#' @param rho Density matrix
#' @param Q_list List of POVMs
#' @return List of probability vectors
born_probs_list <- function(rho, Q_list) {
  lapply(Q_list, function(Q_a) born_probs_setting(rho, Q_a))
}

## ==========================================================================
## SAMPLING FUNCTIONS
## ==========================================================================

#' Sample outcome from Born rule at a given setting
#' @param rho Density matrix
#' @param Q_a List of POVM effects for setting a
#' @return Sampled outcome index b
sample_outcome <- function(rho, Q_a) {
  probs <- born_probs_setting(rho, Q_a)
  sample.int(length(probs), size = 1, prob = probs)
}

#' Sample sequence of outcomes given settings
#' @param a_seq Vector of setting indices
#' @param rho_true True density matrix
#' @param Q_list List of POVMs
#' @param seed Random seed (optional)
#' @return Vector of outcome indices
sample_outcome_sequence <- function(a_seq, rho_true, Q_list, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  # Pre-compute probabilities for all settings
  prob_list <- born_probs_list(rho_true, Q_list)

  b_seq <- integer(length(a_seq))
  for (m in seq_along(a_seq)) {
    a <- a_seq[m]
    probs <- prob_list[[a]]
    b_seq[m] <- sample.int(length(probs), size = 1, prob = probs)
  }

  b_seq
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Example usage (uncomment to test)
# cat("\n=== Library 1 (Pauli PVMs) ===\n")
# L1 <- build_library_1q_L1()
# cat("k =", L1$k, "settings,", L1$M, "total cells\n")
# cat("Settings:", paste(L1$setting_labels, collapse=", "), "\n")
# val1 <- validate_library_1q(L1)
# cat("Valid:", val1$valid, "\n")
#
# cat("\n=== Library 2 (Nine-axis) ===\n")
# L2 <- build_library_1q_L2()
# cat("k =", L2$k, "settings,", L2$M, "total cells\n")
#
# cat("\n=== Library 3 (Random 4-axis) ===\n")
# L3 <- build_library_1q_L3(seed=42)
# cat("k =", L3$k, "settings,", L3$M, "total cells\n")
# cat("Random directions:\n")
# for (nm in names(L3$directions)) {
#   cat(sprintf("  %s: (%.3f, %.3f, %.3f)\n", nm,
#               L3$directions[[nm]][1], L3$directions[[nm]][2], L3$directions[[nm]][3]))
# }

cat("03_measurement_library_1q.R loaded successfully\n")
