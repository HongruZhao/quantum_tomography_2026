# Fisher information and loss metric matrices (Frobenius, Bures, observable-expectation)

This file specifies how to compute:

- the classical Fisher information matrices \(I_a(\theta)\),
- the accumulated plug-in information \(J_n\),
- the loss metric matrix \(G(\hat\theta_n)\) for the three losses.

These objects determine both:
- the proxy risk \(R_n(G)=\mathrm{tr}(GJ_n^{-1})\),
- the adaptive selection scores.

---

## 1) Fisher information in Bloch coordinates

Given a measurement setting \(a\) with effects \(\{Q_{a,b}\}_{b=1}^{r_a}\), define:

- \(s_j(a,b)=\mathrm{tr}(\sigma_j Q_{a,b})\)
- \(c_{a,b}=\mathrm{tr}(Q_{a,b})/N\)

Then:
\[
p_{a,b}(\theta)=c_{a,b}+\frac12\sum_{j=1}^d \theta_j s_j(a,b).
\]

The per-shot Fisher information at setting \(a\) is:
\[
[I_a(\theta)]_{jk}=\frac14\sum_{b=1}^{r_a}\frac{s_j(a,b)s_k(a,b)}{p_{a,b}(\theta)}.
\]

### Implementation in CVXR1121.R
Already implemented via:
- `build_Sab_cab()`
- `fisher_info_setting()`

## 2) Total plug-in Fisher information

Let \(N_a=\sum_b N_{a,b}\) be the total shots at setting \(a\).
Define the total plug-in FI accumulated up to time \(n\):
\[
J_n = \sum_{a=1}^k N_a I_a(\hat\theta_n).
\]

Already implemented by `total_fisher_info(theta, counts, ...)`.

---

## 3) Proxy risk used for plots and design

For a given loss with metric matrix \(G(\hat\theta_n)\),
\[
R_n(G)=\mathrm{tr}\big(G(\hat\theta_n)\,J_n^{-1}\big).
\]

This is computed for every \(n\) after initialization (and can be computed with ridge before invertibility if desired).

---

## 4) Frobenius (Hilbert–Schmidt) metric matrix \(G_F\)

### Definition
For squared Frobenius loss \(L_F(\rho,\hat\rho)=\|\rho-\hat\rho\|_F^2\), with Bloch model
\(\rho(\theta)=I/N + \tfrac12\sum_j \theta_j\sigma_j\), we have:
\[
\|\rho(\theta)-\rho(\theta')\|_F^2
=\frac14(\theta-\theta')^\top \Big[\mathrm{tr}(\sigma_j\sigma_k)\Big]_{jk}\,(\theta-\theta').
\]
Therefore a valid Frobenius metric matrix in \(\theta\)-coordinates is:
\[
G_F = \frac14\Big[\mathrm{tr}(\sigma_j\sigma_k)\Big]_{jk}.
\]

If the basis satisfies \(\mathrm{tr}(\sigma_j\sigma_k)=2\delta_{jk}\), then \(G_F=\tfrac12 I_d\).
For the 2-qubit scaled Pauli-product basis in `02_state_basis_one_two_qubit.md`,
this becomes a constant scalar times identity as well.

### Required function
```r
metric_frobenius <- function(sigmas) {
  d <- length(sigmas)
  G <- matrix(0, d, d)
  for (j in 1:d) for (k in 1:d) {
    G[j,k] <- Re(sum(diag(sigmas[[j]] %*% sigmas[[k]]))) / 4
  }
  G
}
```

---

## 5) Bures / SLD metric matrix \(G_B(\theta)\)

### Lyapunov operator and its inverse
For full-rank \(\rho\), define:
\[
\Omega_\rho(X)=\frac12(\rho X + X\rho).
\]
On the traceless Hermitian subspace, \(\Omega_\rho\) is invertible for \(\rho\succ 0\).

In an eigenbasis \(\rho=U\mathrm{diag}(\lambda)U^\dagger\):
\[
[\Omega_\rho^{-1}(Y)]_{ij} = \frac{2}{\lambda_i+\lambda_j} Y_{ij}.
\]

### Bures metric in Bloch coordinates
Because \(\partial_{\theta_j}\rho(\theta)=\sigma_j/2\), the metric matrix is:
\[
[G_B(\theta)]_{jk} = \frac{1}{16}\mathrm{tr}\big(\sigma_j \Omega_{\rho(\theta)}^{-1}(\sigma_k)\big).
\]

### Required functions

#### (a) Lyapunov inverse via eigen-decomposition
```r
lyapunov_inv <- function(rho, Y, tol = 1e-12) {
  # rho: NxN Hermitian, assumed full rank
  # Y:   NxN (typically Hermitian)
  eg <- eigen(rho)
  U <- eg$vectors
  lam <- pmax(Re(eg$values), tol)  # guard against numerical noise
  Ytilde <- Conj(t(U)) %*% Y %*% U
  # elementwise division by (lam_i + lam_j)/2  -> multiply by 2/(lam_i+lam_j)
  denom <- outer(lam, lam, "+")
  Xtilde <- (2/denom) * Ytilde
  X <- U %*% Xtilde %*% Conj(t(U))
  hermitianize(X)
}
```

#### (b) Bures metric matrix
```r
metric_bures <- function(theta, sigmas, N, tol = 1e-12) {
  rho <- rho_of_theta(theta, sigmas, N)
  d <- length(sigmas)
  G <- matrix(0, d, d)
  # precompute Omega^{-1}(sigma_k) for each k
  Omega_inv_sig <- vector("list", d)
  for (k in 1:d) Omega_inv_sig[[k]] <- lyapunov_inv(rho, sigmas[[k]], tol = tol)
  for (j in 1:d) for (k in 1:d) {
    G[j,k] <- Re(sum(diag(sigmas[[j]] %*% Omega_inv_sig[[k]]))) / 16
  }
  (G + t(G))/2
}
```

### Stability requirement
This computation assumes \(\rho(\hat\theta_n)\) is full rank. Therefore:
- the stabilized MLE constraint \(\hat\rho_n\succeq\eta I\) is strongly recommended,
- else add a spectral floor `tol` when forming denominators.

---

## 6) Observable-expectation metric matrix \(G_O\)

### Loss definition
Given fixed Hermitian observables \(O_1,\dots,O_m\),
\[
L_O(\rho,\hat\rho)=\sum_{j=1}^m \big(\mathrm{tr}((\rho-\hat\rho)O_j)\big)^2.
\]

### Metric matrix in Bloch coordinates
Define vectors \(o_j\in\mathbb R^d\) by:
\[
[o_j]_k=\mathrm{tr}(\sigma_k O_j).
\]
Then:
\[
G_O=\frac14\sum_{j=1}^m o_j o_j^\top \succeq 0,
\]
which is constant (does not depend on \(\theta\)).

### Completeness requirement
To make the loss identifying, the traceless parts of \(\{O_j\}\) should span \(\mathrm{Herm}_0(\mathbb C^N)\).

### Required function
```r
metric_observable <- function(sigmas, Obs) {
  d <- length(sigmas)
  m <- length(Obs)
  G <- matrix(0, d, d)
  for (j in 1:m) {
    oj <- numeric(d)
    for (k in 1:d) oj[k] <- Re(sum(diag(sigmas[[k]] %*% Obs[[j]])))
    G <- G + tcrossprod(oj, oj)
  }
  G / 4
}
```

### Default observable families for this project (recommended)

To ensure completeness without overthinking:
- **One qubit:** \(O=\{X,Y,Z\}\).
- **Two qubits:** the 15 Pauli products \(\sigma_\alpha\otimes\sigma_\beta\) excluding \(I\otimes I\)
  (same ordering as Library A).

If you instead want task-driven behavior (different from Frobenius), choose a smaller set such as:
- single-qubit marginals: \(X\otimes I, Y\otimes I, Z\otimes I, I\otimes X, I\otimes Y, I\otimes Z\)
and accept that \(G_O\) will be rank-deficient (then the proxy risk targets only that subspace).

---

## 7) Proxy risk computation helper

Implement:

```r
proxy_risk <- function(G, J, ridge = 0) {
  d <- nrow(J)
  Juse <- J + ridge*diag(d)
  solveJ <- solve(Juse)
  sum(diag(G %*% solveJ))
}
```

During initialization, use a small `ridge > 0` if `J` is singular.

