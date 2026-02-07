# Existing R code inventory (CVXR1121.R) and required extensions

This file documents the functions that already exist in **CVXR1121.R** and how they map to the
notation in the PDF. It also specifies minimal extensions required to match the simulation study in
**Jan29_simple_simulation (1).pdf**.

## 1) Utility helpers (already implemented)

### `traceC(A)`
- Complex-safe trace: `Re(sum(diag(A)))`.

### `hermitianize(A)`
- Enforces Hermiticity: \((A + A^\dagger)/2\).

### `real_embed(M)`
Implements the standard real embedding \( \mathbb C^{N\times N} \to \mathbb R^{2N\times 2N}\):
\[
M \mapsto \begin{bmatrix}\Re M & -\Im M\\ \Im M & \Re M\end{bmatrix}.
\]
Used to express complex PSD constraints inside CVXR as real PSD constraints.

## 2) Basis construction (already implemented for generic SU(N))

### `build_suN_basis(N)`
Builds Kimura-style traceless Hermitian generators \(\{\sigma_j\}_{j=1}^{N^2-1}\) satisfying
\(\operatorname{tr}(\sigma_j\sigma_k)=2\delta_{jk}\). This is the construction in the PDF (Eqs. 10–14).

**Note:** For this project we will also implement a *Pauli / Pauli-product* basis builder, because the
simulation section in the PDF specifies Pauli bases for N=2 and Pauli-product bases for N=4.

## 3) Measurement decomposition (already implemented for observables → spectral projectors)

### `spectral_projectors(B, tol)`
Given a Hermitian operator \(B\), compute its spectral projectors \(Q_{b}\), grouping eigenvalues within tolerance.

### `build_measurements_from_basis(sigmas)`
Convenience: treats each \(\sigma_a\) as an observable \(B_a\) and builds the corresponding PVM projectors.

**Important:** The simulation libraries in the PDF are *not* “measure each sigma in the SU(N) basis”.
Instead, they are Pauli-axis / Pauli-parity / MUB-basis measurements. We therefore need custom
library builders that directly return `Q_list`, `ab_df`, `ab_row`.

## 4) Bloch parametrization and Born probabilities (already implemented)

### `rho_of_theta(theta, sigmas, N)`
Implements the affine expansion:
\[
\rho(\theta) = \frac{I}{N} + \frac12\sum_{j=1}^d \theta_j \sigma_j.
\]
This matches the PDF Bloch model (Eq. 15). The only caveat is basis scaling: if you choose a basis
with different normalization, it is still valid, but the meaning of \(\theta\) changes accordingly.

### `born_probs_list(rho, Q_list)`
Returns the list \(\{p_{a,\cdot}\}_{a=1}^k\) where
\(p_{a,b}=\operatorname{tr}(\rho Q_{a,b})\) and each \(p_{a,\cdot}\) is normalized to sum to 1.

### `sample_b_sequence(a_seq, prob_by_a)`
Samples outcomes \(b_t\) given selected settings \(a_t\) and probabilities.

## 5) Affine probability model in theta (already implemented)

### `build_Sab_cab(sigmas, Q_list, N, ab_df)`
Computes:
- \(c_{a,b}=\operatorname{tr}(Q_{a,b})/N\)
- \(s_j(a,b)=\operatorname{tr}(\sigma_j Q_{a,b})\)

and stores them in:
- `c_ab` (length M)
- `S_ab` (M × d)

so that
\[
p_{a,b}(\theta)=c_{a,b}+\frac12\,s(a,b)^\top \theta
\]
as in the PDF (Eq. 18).

## 6) CVXR MLE in theta (standardized to CVXR1121-style with fixed floor)

### `fit_theta_cvxr(N, sigmas, S_ab, c_ab, Nab, eta=1e-3, solver, eps_log, verbose)`
Solves:
\[
\min_\theta -\sum_{a,b} N_{a,b}\log\big(c_{a,b}+\tfrac12 S_{a,b}\theta + \varepsilon_{\log}\big)
\quad\text{s.t.}\quad \rho(\theta)\succeq 0
\]
where PSD is imposed via real embedding and a CVXR PSD variable.

The production solver enforces the eigenvalue floor constraint
\[
\rho(\theta)\succeq \eta I
\]
with \(\eta=10^{-3}\) by default (and never below \(10^{-3}\)).

This follows the CVXR1121 formulation with:
- affine probabilities \(p(\theta)=c+\tfrac12S\theta\),
- CVXR PSD variable for \(\rho_{\mathbb R}(\theta)\),
- explicit floor constraint \(\rho_{\mathbb R}(\theta)-\eta I_{\mathbb R}\succeq 0\),
- probability constraint \(p(\theta)\ge 0\).

## 7) Fisher information and design selectors (partially implemented)

### Fisher information
Already implemented:

- `fisher_info_setting(theta, a, S_ab, c_ab, ab_df)`  
  Implements the PDF FI formula (Eq. 19):
  \[
  I_a(\theta)=\frac14 S_a^\top \operatorname{diag}(1/p_{a,\cdot}(\theta)) S_a.
  \]

- `total_fisher_info(theta, counts, S_ab, c_ab, ab_df)`  
  Implements \(J_n=\sum_a N_a I_a(\theta)\) (PDF Eq. 24).

### Adaptive selection
Already implemented, but only for *unweighted* A-optimality and D-optimality:

- `select_next_setting(theta, ..., criterion = c("D","A"), method = c("exact","first-order"))`

**Required extension:** the PDF simulation uses **metric-weighted A-optimality** with
\(R_n(G)=\operatorname{tr}(GJ_n^{-1})\) and the two greedy rules:

- Exact: \(\arg\min_a \operatorname{tr}(G(J_n+I_a)^{-1})\)
- GI1:   \(\arg\max_a \operatorname{tr}(GJ_n^{-1}I_aJ_n^{-1})\)

for three different choices of \(G\): Frobenius, Bures, observable-expectation.

We will implement:

- `select_next_setting_metric(theta, ..., G, method = c("exact","GI1"), ridge=...)`
- `adaptive_design_sequence_metric(...)`

so that `G` can be:
- constant (Frobenius, observable),
- dependent on current theta (Bures).

## 8) Simulation-loop helpers already present (optional reuse)

The following exist and can be reused or adapted:

- `uniform_design_sequence(...)`  
- `adaptive_design_sequence(...)` (currently D/A only, not metric-weighted)
- plotting helpers (`plot_*`)

For this project, it is acceptable to:
- reuse their structure, but
- replace selection criteria with the metric-weighted scores required by the PDF.

## Summary of “must-add” functions

Minimum new functionality to match the PDF simulation section:

1. **Pauli / Pauli-product basis builder** (N=2 and N=4).
2. **Measurement library builders** for Libraries 1–3 and A–B.
3. **Stabilized MLE**: use the CVXR1121-style solver with eigenvalue floor \(\eta=10^{-3}\).
4. **Loss metric matrix functions**:
   - `G_frobenius(sigmas)` (constant),
   - `G_bures(theta, sigmas, N)` (state-dependent),
   - `G_observable(sigmas, Obs)` (constant).
5. **Metric-weighted selection rules** (exact and GI1).
6. **Initialization routine** to ensure \(J_n\) invertible (PDF Eq. 25).
7. **Simulation controller** producing 15 Monte Carlo mean risk plots.
