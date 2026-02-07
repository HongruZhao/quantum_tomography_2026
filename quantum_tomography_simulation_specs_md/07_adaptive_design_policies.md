# Adaptive design policies (uniform, exact, GI1) and initialization for invertible Fisher information

This file specifies the sequential measurement-design rules used in the simulation.

Notation:
- \(k\): number of available measurement settings in the chosen library.
- \(I_a(\theta)\): per-shot Fisher information at setting \(a\).
- \(N_a\): number of shots allocated to setting \(a\) so far.
- \(J_n = \sum_a N_a I_a(\hat\theta_n)\): accumulated plug-in information.
- \(G(\hat\theta_n)\): loss metric matrix for the chosen loss.

We compare three policies within each library:

1. Uniform sampling.
2. Adaptive exact one-step lookahead.
3. Adaptive GI1 (first-order proxy).

---

## 1) Initialization phase: ensure \(J_n\) invertible

### Why initialization is needed
For small \(n\), \(J_n\) can be singular if the trajectory has not explored enough settings to identify all
\(d=N^2-1\) Bloch coordinates. Adaptive rules require \(J_n^{-1}\).

### Deterministic initialization (current implementation)
Initialization is deterministic and library-specific, chosen to span all degrees of freedom so \(J_n\) is
invertible by construction:
- **L1:** Pauli \(X, Y, Z\)
- **L2:** Pauli \(X, Y, Z\) (first 3 of 9)
- **L3:** all 4 random axes (fixed seed)
- **A:** all 15 Pauli-parity PVMs
- **B:** all 5 MUB bases

No eigenvalue threshold is used in the current implementation.

### Ridge inverse during initialization (optional)
If you must compute selection scores before invertibility, replace:
\[
J_n^{-1}\ \leadsto\ (J_n+\gamma I)^{-1}
\]
for a small \(\gamma>0\). Once invertible, set \(\gamma=0\) permanently.

---

## 2) Policy 1: Uniform sampling

At each time \(t\), pick:
\[
a_{t}\sim\text{Unif}\{1,\dots,k\}.
\]

No dependence on \(\hat\theta_t\).

---

## 3) Policy 2: Adaptive exact one-step lookahead

After initialization (so \(J_n\) invertible), define the proxy risk:
\[
R_n(G) = \mathrm{tr}(G(\hat\theta_n)\,J_n^{-1}).
\]

The exact one-step lookahead chooses the next setting by minimizing the risk after adding one more
shot at setting \(a\):
\[
a_{n+1}^{\mathrm{exact}} \in \arg\min_{a\in[k]}\ \mathrm{tr}\Big(G(\hat\theta_n)\,(J_n + I_a(\hat\theta_n))^{-1}\Big).
\]

### Implementation notes
- Compute \(I_a(\hat\theta_n)\) for all settings (already supported).
- For each candidate setting a:
  - form `Jcand = Jn + Ia`,
  - solve `Jcand^{-1}`,
  - compute score `tr(G %*% solve(Jcand))`,
  - pick argmin.

---

## 4) Policy 3: Adaptive GI1 (first-order proxy)

Use the first-order approximation:
\[
(J_n + I_a)^{-1}\approx J_n^{-1} - J_n^{-1}I_aJ_n^{-1},
\]
drop terms independent of \(a\), and choose:
\[
a_{n+1}^{\mathrm{GI1}} \in \arg\max_{a\in[k]}\ \mathrm{tr}\Big(G(\hat\theta_n)\,J_n^{-1}I_a(\hat\theta_n)J_n^{-1}\Big).
\]

### Implementation notes
- Precompute `Jinv = solve(Jn)` once.
- For each candidate `Ia`, compute scalar score:
  - `sum(diag(G %*% Jinv %*% Ia %*% Jinv))`
- Pick argmax.

GI1 is much cheaper than exact because it avoids inverting \(k\) different matrices.

---

## 5) Specialization: how \(G\) is chosen for each loss

This is where we differ from the existing CVXR1121.R A/D-opt code.

### (a) Frobenius loss
Use \(G=G_F\) (constant).
For an orthonormal basis, \(G_F=\tfrac12 I\) and both exact and GI1 reduce to the standard A-optimal
rules:
- exact: minimize \(\mathrm{tr}((J_n+I_a)^{-1})\)
- GI1: maximize \(\mathrm{tr}(J_n^{-1} I_a J_n^{-1})\)

### (b) Bures loss
Use \(G=G_B(\hat\theta_n)\) computed from \(\rho(\hat\theta_n)\) via the Lyapunov inverse.
This is state-dependent and must be recomputed when \(\hat\theta_n\) updates.

### (c) Observable-expectation loss
Use \(G=G_O\) computed once from the chosen observables \(\{O_j\}\).

---

## 6) Required new selector API

Implement a metric-aware selector:

```r
select_next_setting_metric <- function(
  theta_hat, counts,
  S_ab, c_ab, ab_df,
  G,
  method = c("exact","GI1"),
  ridge = 0
) {
  # returns: list(a_next = ..., scores = numeric(k))
}
```

Key requirements:

- `G` is a numeric d×d matrix.
- The function computes:
  - `Ia_list[a]` = fisher_info_setting(theta_hat, a, ...)
  - `Jn` = total_fisher_info(theta_hat, counts, ...)
- If `method == "exact"`:
  - score[a] = tr( G %*% solve(Jn + Ia_list[[a]] + ridge*I) )
  - choose argmin
- If `method == "GI1"`:
  - Jinv = solve(Jn + ridge*I)
  - score[a] = tr( G %*% Jinv %*% Ia_list[[a]] %*% Jinv )
  - choose argmax

---

## 7) Required adaptive loop wrapper

Implement:

```r
adaptive_design_sequence_metric <- function(
  theta_init, rho_sampling,
  Q_list, S_ab, c_ab, ab_df,
  n_total,
  metric_fun,           # function(theta_hat) -> G matrix
  policy = c("uniform","exact","GI1"),
  ridge_init = 1e-6,
  eta_mle = 1e-3,
  solver = "SCS",
  check_every = 1,
  seed = NULL
) {
  # returns: list(a_seq, b_seq, Nab_hist?, theta_hat_hist?, risk_hist)
}
```

Behavior:

1. Run deterministic initialization with the library-specific fixed settings.
2. For t > n0:
   - choose `a_t` by the chosen policy:
     - uniform: random
     - exact / GI1: call `select_next_setting_metric()` with current `G`
   - sample outcome `b_t` from `rho_sampling`
   - update counts
   - refit MLE (`fit_theta_cvxr(..., eta=eta_mle)`)
   - compute risk `R_t(G)` and store.

---

## 8) Output expectations

For plotting, we need at minimum:

- `risk[n]` for n = 1..n_total (or `NA` before initialization if you prefer).
- optionally `theta_hat[n,]` (for debugging and extra plots).
- counts or empirical design weights (optional).

The simulation controller will average `risk[n]` across Monte Carlo replicates and plot.
