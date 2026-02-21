## ============================================================================
## 04_measurement_library_2q.R
## Two-qubit measurement libraries (Libraries A and B)
## Based on: 04_measurement_libraries_two_qubit.md
## ============================================================================

source("03_measurement_library_1q.R")

## ==========================================================================
## LIBRARY A: 15 Two-qubit Pauli-parity PVMs (binary, 2+2 rank)
## ==========================================================================

#' Build Library A: Binary Pauli-parity measurements
#' Measures all 15 non-identity two-qubit Pauli operators
#' Each P ∈ {I,X,Y,Z}⊗2 \ {I⊗I} has P² = I, spectrum {+1,+1,-1,-1}
#' Q_+ = (I+P)/2, Q_- = (I-P)/2
#' k = 15 settings, each with r = 2 outcomes
#' @param order Character vector specifying ordering of 15 Pauli products
#' @return Standardized measurement library object
build_library_2q_A <- function(
    order = c("IX","IY","IZ","XI","XX","XY","XZ","YI","YX","YY","YZ","ZI","ZX","ZY","ZZ")
) {
  pauli <- build_pauli_basis_1q()
  single <- list(I = pauli$I, X = pauli$sigmas$X, Y = pauli$sigmas$Y, Z = pauli$sigmas$Z)
  I4 <- diag(4) + 0i

  P_list <- list()
  Q_list <- list()

  for (lab in order) {
    a <- substr(lab, 1, 1)
    b <- substr(lab, 2, 2)
    # Two-qubit Pauli operator
    P <- kron(single[[a]], single[[b]])
    P_list[[lab]] <- P

    # Binary PVM from parity measurement
    Qp <- 0.5 * (I4 + P)  # eigenvalue +1
    Qm <- 0.5 * (I4 - P)  # eigenvalue -1
    Q_list[[lab]] <- list(hermitianize(Qp), hermitianize(Qm))
  }

  lib <- build_ab_indexing(Q_list, setting_labels = order, N = 4)
  lib$library_name <- "A_PauliParity"
  lib$description <- "15 two-qubit Pauli-parity PVMs (binary, 2+2 rank)"
  lib$P_list <- P_list
  lib$Obs <- P_list
  lib$pauli_1q <- pauli

  lib
}

## ==========================================================================
## LIBRARY B: 5 MUB bases (rank-1 PVM, 4 outcomes each)
## ==========================================================================

#' Build Library B: 5 Mutually Unbiased Bases (MUBs)
#' For N=4, there exist N+1=5 MUBs
#' Each basis is a 4-outcome rank-1 PVM
#' k = 5 settings, each with r = 4 outcomes
#' @return Standardized measurement library object
build_library_2q_B <- function() {
  pauli <- build_pauli_basis_1q()

  # One-qubit basis states
  ket0 <- matrix(c(1, 0), ncol = 1) + 0i
  ket1 <- matrix(c(0, 1), ncol = 1) + 0i

  # |+⟩, |-⟩ (X eigenstates)
  ketp <- (ket0 + ket1) / sqrt(2)
  ketm <- (ket0 - ket1) / sqrt(2)

  # |+i⟩, |-i⟩ (Y eigenstates)
  ketpi <- (ket0 + 1i * ket1) / sqrt(2)
  ketmi <- (ket0 - 1i * ket1) / sqrt(2)

  # Two-qubit computational basis states
  ket00 <- kron(ket0, ket0)
  ket01 <- kron(ket0, ket1)
  ket10 <- kron(ket1, ket0)
  ket11 <- kron(ket1, ket1)

  # ===== Five MUB Bases =====

  # Basis 1: Z⊗Z (computational basis)
  BZ <- list(ket00, ket01, ket10, ket11)

  # Basis 2: X⊗X (product basis)
  BX <- list(
    kron(ketp, ketp),  # |++⟩
    kron(ketp, ketm),  # |+-⟩
    kron(ketm, ketp),  # |-+⟩
    kron(ketm, ketm)   # |--⟩
  )

  # Basis 3: Y⊗Y (product basis)
  BY <- list(
    kron(ketpi, ketpi),  # |+i,+i⟩
    kron(ketpi, ketmi),  # |+i,-i⟩
    kron(ketmi, ketpi),  # |-i,+i⟩
    kron(ketmi, ketmi)   # |-i,-i⟩
  )

  # Basis 4: Entangled MUB (BE1)
  BE1 <- list(
    0.5 * (ket00 + ket01 + 1i*ket10 - 1i*ket11),
    0.5 * (ket00 + ket01 - 1i*ket10 + 1i*ket11),
    0.5 * (ket00 - ket01 + 1i*ket10 + 1i*ket11),
    0.5 * (ket00 - ket01 - 1i*ket10 - 1i*ket11)
  )

  # Basis 5: Entangled MUB (BE2)
  BE2 <- list(
    0.5 * (ket00 + 1i*ket01 - ket10 + 1i*ket11),
    0.5 * (ket00 - 1i*ket01 + ket10 + 1i*ket11),
    0.5 * (ket00 + 1i*ket01 + ket10 - 1i*ket11),
    0.5 * (ket00 - 1i*ket01 - ket10 - 1i*ket11)
  )

  bases <- list(BZ = BZ, BX = BX, BY = BY, BE1 = BE1, BE2 = BE2)

  # Build Q_list: each setting a is a list of 4 rank-1 projectors
  Q_list <- lapply(bases, function(B) {
    lapply(B, function(psi) hermitianize(proj_rank1(psi)))
  })

  # Build observables with distinct eigenvalues for each basis
  eigs <- c(3, 1, -1, -3)
  O_list <- lapply(bases, function(B) {
    obs <- matrix(0, 4, 4)
    for (i in 1:4) {
      obs <- obs + eigs[i] * hermitianize(proj_rank1(B[[i]]))
    }
    hermitianize(obs)
  })

  lib <- build_ab_indexing(Q_list, setting_labels = names(bases), N = 4)
  lib$library_name <- "B_MUB5"
  lib$description <- "5 Mutually Unbiased Bases (4 outcomes each)"
  lib$bases <- bases
  lib$Obs <- O_list
  lib$Obs_eigs <- eigs
  lib$pauli_1q <- pauli

  lib
}

## ==========================================================================
## VALIDATION FUNCTIONS
## ==========================================================================

#' Validate a two-qubit measurement library
#' @param lib Measurement library object
#' @param tol Numerical tolerance
#' @return List with validation results
validate_library_2q <- function(lib, tol = 1e-10) {
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

#' Validate MUB unbiasedness for Library B
#' Check |⟨ψ|φ⟩|² = 1/4 for states from different bases
#' @param lib Library B object
#' @param tol Numerical tolerance
#' @return List with validation results
validate_mub_unbiasedness <- function(lib, tol = 1e-10) {
  if (is.null(lib$bases)) {
    return(list(valid = FALSE, message = "Library does not contain basis vectors"))
  }

  bases <- lib$bases
  k <- length(bases)
  expected_overlap <- 1/4

  all_overlaps <- list()
  max_deviation <- 0

  for (a1 in 1:(k-1)) {
    for (a2 in (a1+1):k) {
      basis1 <- bases[[a1]]
      basis2 <- bases[[a2]]

      overlaps <- matrix(0, 4, 4)
      for (i in 1:4) {
        for (j in 1:4) {
          overlap <- Mod(Conj(t(basis1[[i]])) %*% basis2[[j]])^2
          overlaps[i, j] <- overlap
          deviation <- abs(overlap - expected_overlap)
          max_deviation <- max(max_deviation, deviation)
        }
      }
      all_overlaps[[paste0(names(bases)[a1], "_", names(bases)[a2])]] <- overlaps
    }
  }

  list(
    valid = max_deviation < tol,
    max_deviation = max_deviation,
    expected = expected_overlap,
    overlaps = all_overlaps
  )
}

#' Validate basis orthonormality for Library B
#' Check that each basis is orthonormal
#' @param lib Library B object
#' @param tol Numerical tolerance
#' @return List with validation results
validate_basis_orthonormality_2q <- function(lib, tol = 1e-10) {
  if (is.null(lib$bases)) {
    return(list(valid = FALSE, message = "Library does not contain basis vectors"))
  }

  bases <- lib$bases
  results <- list()
  all_valid <- TRUE

  for (nm in names(bases)) {
    B <- bases[[nm]]
    # Form matrix with basis vectors as columns
    Psi <- do.call(cbind, B)
    # Check Ψ†Ψ = I
    gram <- Conj(t(Psi)) %*% Psi
    deviation <- max(Mod(gram - diag(4)))

    results[[nm]] <- list(
      orthonormal = deviation < tol,
      max_deviation = deviation
    )

    if (deviation >= tol) all_valid <- FALSE
  }

  list(valid = all_valid, basis_results = results)
}

## ==========================================================================
## EXAMPLES
## ==========================================================================

# Example usage (uncomment to test)
# cat("\n=== Library A (Pauli Parity) ===\n")
# LA <- build_library_2q_A()
# cat("k =", LA$k, "settings,", LA$M, "total cells\n")
# cat("Settings:", paste(LA$setting_labels[1:5], collapse=", "), "...\n")
# valA <- validate_library_2q(LA)
# cat("POVM valid:", valA$valid, "\n")
#
# cat("\n=== Library B (5 MUBs) ===\n")
# LB <- build_library_2q_B()
# cat("k =", LB$k, "settings,", LB$M, "total cells\n")
# cat("Settings:", paste(LB$setting_labels, collapse=", "), "\n")
# valB <- validate_library_2q(LB)
# cat("POVM valid:", valB$valid, "\n")
#
# # Check MUB properties
# mub_val <- validate_mub_unbiasedness(LB)
# cat("MUB unbiasedness valid:", mub_val$valid, "\n")
# cat("Max overlap deviation:", mub_val$max_deviation, "\n")
#
# ortho_val <- validate_basis_orthonormality_2q(LB)
# cat("Basis orthonormality valid:", ortho_val$valid, "\n")

cat("04_measurement_library_2q.R loaded successfully\n")
