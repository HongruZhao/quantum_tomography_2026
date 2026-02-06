# One-qubit measurement libraries (Libraries 1–3)

This file specifies how to build the three one-qubit measurement libraries used in the simulation.

All three libraries operate on \(N=2\) and use 2-outcome projective measurements (binary PVMs).
Each setting corresponds to measuring along a Bloch-sphere direction \(u\in\mathbb R^3\), \(\|u\|=1\),
with projectors:
\[
Q_{u,+}=\tfrac12(I+u\cdot\sigma),\qquad
Q_{u,-}=\tfrac12(I-u\cdot\sigma),
\]
where \(u\cdot\sigma = u_x X + u_y Y + u_z Z\).

## Shared helper: projectors for a Bloch direction

Implement:

```r
pvm_projectors_for_direction_1q <- function(u, pauli) {
  # u: numeric length-3, not necessarily unit
  # pauli: output of build_pauli_basis_1q()
  u <- as.numeric(u)
  u <- u / sqrt(sum(u^2))
  X <- pauli$sigmas$X; Y <- pauli$sigmas$Y; Z <- pauli$sigmas$Z
  I2 <- pauli$I

  B <- u[1]*X + u[2]*Y + u[3]*Z  # Hermitian observable u·σ
  Qp <- 0.5*(I2 + B)
  Qm <- 0.5*(I2 - B)
  list(Qp = Qp, Qm = Qm, u = u, B = B)
}
```

All libraries below return a standardized measurement-library object:

```r
list(
  N = 2,
  k = k,                     # number of settings
  r_vec = rep(2, k),          # outcomes per setting
  setting_labels = ...,       # length k
  Q_list = ...,               # list length k, each element is list of 2 projectors
  ab_df = ...,                # data frame (a,b,row)
  ab_row = ...                # list: ab_row[[a]][b] gives row index in Nab
)
```

You can reuse the `ab_df/ab_row` construction pattern from `build_measurements_from_basis()`.

---

## Library 1: Pauli PVMs (k1 = 3)

Settings correspond to measuring \(X\), \(Y\), \(Z\), i.e. directions:
- \(e_x=(1,0,0)\)
- \(e_y=(0,1,0)\)
- \(e_z=(0,0,1)\)

Specification:

```r
build_library_1q_L1 <- function() {
  pauli <- build_pauli_basis_1q()
  U <- list(ex=c(1,0,0), ey=c(0,1,0), ez=c(0,0,1))
  Q_list <- lapply(U, function(u) {
    P <- pvm_projectors_for_direction_1q(u, pauli)
    list(P$Qp, P$Qm)  # outcomes b=1 (+), b=2 (-)
  })
  setting_labels <- names(U)
  build_ab_indexing(Q_list, setting_labels = setting_labels, N = 2)
}
```

---

## Library 2: Nine-axis PVMs (k2 = 9)

The axis set:
\[
U_9=\{e_x,e_y,e_z,e_x\pm e_y,e_x\pm e_z,e_y\pm e_z\}.
\]
For each direction \(u\in U_9\), use the normalized unit vector \(\hat u = u/\|u\|_2\).

Specification:

```r
build_library_1q_L2 <- function() {
  pauli <- build_pauli_basis_1q()
  U <- list(
    ex=c(1,0,0), ey=c(0,1,0), ez=c(0,0,1),
    ex_py=c(1,1,0), ex_my=c(1,-1,0),
    ex_pz=c(1,0,1), ex_mz=c(1,0,-1),
    ey_pz=c(0,1,1), ey_mz=c(0,1,-1)
  )
  Q_list <- lapply(U, function(u) {
    P <- pvm_projectors_for_direction_1q(u, pauli)
    list(P$Qp, P$Qm)
  })
  setting_labels <- names(U)
  build_ab_indexing(Q_list, setting_labels = setting_labels, N = 2)
}
```

The nine observables are:
\[
X,\;Y,\;Z,\;\frac{X\pm Y}{\sqrt2},\;\frac{X\pm Z}{\sqrt2},\;\frac{Y\pm Z}{\sqrt2}.
\]

---

## Library 3: Random 4-axis PVMs (k3 = 4)

Draw four random unit vectors \(u^{(1)},\dots,u^{(4)}\) uniformly from the sphere and fix them for
the whole simulation. Use the same projector formula.

Specification:

```r
sample_random_unit_vector_3d <- function() {
  z <- rnorm(3)
  z / sqrt(sum(z^2))
}

build_library_1q_L3 <- function(seed = 1) {
  set.seed(seed)
  pauli <- build_pauli_basis_1q()
  U <- lapply(1:4, function(j) sample_random_unit_vector_3d())
  names(U) <- sprintf("u%02d", 1:4)
  Q_list <- lapply(U, function(u) {
    P <- pvm_projectors_for_direction_1q(u, pauli)
    list(P$Qp, P$Qm)
  })
  setting_labels <- names(U)
  build_ab_indexing(Q_list, setting_labels = setting_labels, N = 2, extra = list(U = U))
}
```

### Reproducibility requirement
- `seed` must be stored with the library object so the same random axes can be reused across policies.

---

## Required shared helper: `build_ab_indexing(Q_list, ...)`

Implement a generic helper that constructs:

- `ab_df`: rows `(a,b,row)`
- `ab_row`: list mapping `(a,b)` to flat index
- `r_vec`: number of outcomes per setting

This should be compatible with `build_Sab_cab()` and `counts_from_ab()` in CVXR1121.R.

Suggested implementation:

```r
build_ab_indexing <- function(Q_list, setting_labels = NULL, N = NULL, extra = NULL) {
  k <- length(Q_list)
  r_vec <- vapply(Q_list, length, integer(1))
  rows <- 1:sum(r_vec)
  a_col <- rep(1:k, times = r_vec)
  b_col <- unlist(lapply(r_vec, function(r) 1:r))
  ab_df <- data.frame(a = a_col, b = b_col, row = rows)
  ab_row <- vector("list", k)
  idx <- 1
  for (a in 1:k) {
    ra <- r_vec[a]
    ab_row[[a]] <- idx:(idx+ra-1)
    idx <- idx + ra
  }
  out <- list(k = k, r_vec = r_vec, Q_list = Q_list, ab_df = ab_df, ab_row = ab_row)
  if (!is.null(setting_labels)) out$setting_labels <- setting_labels
  if (!is.null(N)) out$N <- N
  if (!is.null(extra)) out$extra <- extra
  out
}
```

