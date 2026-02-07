# Quantum tomography adaptive-design simulation (specification)

This document set translates **Jan29_simple_simulation (1).pdf** into implementable, R-friendly specifications.
It is intended to be used as input to an automated “markdown → R code” generation workflow, reusing
functions already present in **CVXR1121.R**.

## Scope (what we are implementing)

We simulate sequential quantum-state tomography under a **finite measurement library**, comparing:

- **Uniform sampling** over measurement settings (baseline).
- **Adaptive exact one-step lookahead** (minimizes a metric-weighted A-optimal proxy after adding one shot).
- **Adaptive GI1** (first-order proxy to the exact lookahead).

We do this for:

- **One qubit** (N = 2, d = 3) with **Libraries 1–3** (3, 9, and 4 binary PVM settings).
- **Two qubits** (N = 4, d = 15) with **Library A** (15 binary 2+2 PVMs) and **Library B** (5 rank-1 PVM bases).

Loss / design metric choices:

1. **Frobenius (Hilbert–Schmidt) squared loss**.
2. **Bures squared loss** (via the SLD/Bures Riemannian metric).
3. **Observable-expectation squared loss** (task-driven, user-provided observables).

For each (system, library, loss), we plot the Monte Carlo mean of the proxy risk
\[
R_n(G) := \operatorname{tr}\big(G(\hat\theta_n)\,J_n^{-1}\big)
\]
versus sample size \(n\), with three curves (uniform / exact / GI1).

## High-level architecture

At runtime, each (system, library) defines:

- A Bloch basis \(\{\sigma_j\}_{j=1}^d\) and the affine parametrization \(\rho(\theta)\).
- A measurement library: settings \(a \in [k]\), each with POVM effects \(\{Q_{a,b}\}_{b=1}^{r_a}\).
- A map between (a,b) cells and a flat count vector \(N_{a,b}\) (as in CVXR1121.R).

Per replicate, we simulate:

1. Draw a true full-rank state \(\rho^\star\).
2. Run an initialization phase to ensure the plug-in total Fisher information \(J_n\) is invertible.
3. Run sequential measurements up to \(n_{\max}\), selecting \(a_{n+1}\) by the chosen policy.
4. After each step, compute:
   - stabilized MLE \(\hat\theta_n\) (via CVXR),
   - per-setting FI matrices \(I_a(\hat\theta_n)\),
   - total FI \(J_n = \sum_a N_a I_a(\hat\theta_n)\),
   - loss metric matrix \(G(\hat\theta_n)\),
   - proxy risk \(R_n(G)\).
5. Aggregate \(R_n(G)\) across Monte Carlo replicates.

## File map

- `01_existing_r_functions.md`  
  What is already implemented in CVXR1121.R, how it maps to the PDF notation, and what needs extension.

- `02_state_basis_one_two_qubit.md`  
  Bloch parametrizations and basis construction for N=2 and N=4 (Pauli / Pauli-product).

- `03_measurement_libraries_one_qubit.md`  
  Libraries 1–3 (PVMs on X/Y/Z, 9-axis, random 4-axis).

- `04_measurement_libraries_two_qubit.md`  
  Libraries A and B (15 Pauli-parity measurements; 5 explicit MUB bases).

- `05_mle_cvxr_solver.md`  
  Stabilized MLE with eigenvalue floor \(\rho \succeq \eta I\) in CVXR; how to implement via real embedding.

- `06_fisher_information_and_metrics.md`  
  Fisher information formulas and the three loss metric matrices \(G_F\), \(G_B(\theta)\), \(G_O\).

- `07_adaptive_design_policies.md`  
  Initialization, exact and GI1 scores, ridge handling, and required new selector APIs.

- `08_simulation_controller.md`  
  Simulation parameters, loops, outputs, and plotting conventions.

- `09_validation_and_sanity_checks.md`  
  Optional but strongly recommended checks: POVM validity, FI PSD, MUB unbiasedness, invariances.

## Conventions (important)

### Indices and data structures

We follow the PDF notation but implement using the `ab_df / ab_row / Nab` layout from CVXR1121.R:

- `Q_list[[a]][[b]]` is the matrix \(Q_{a,b}\).
- `ab_df` is a data frame with columns: `a`, `b`, `row` mapping each (a,b) to a row index.
- `Nab[row]` is the count \(N_{a,b}\) for the corresponding (a,b).

### Hermitian complex matrices in R

- Use `0+0i` complex matrices.
- Use `Conj(t(M))` for Hermitian transpose.
- Always enforce Hermiticity with `hermitianize()` (already in CVXR1121.R).

### Solver stability

The implementation uses a stabilized estimator constraint:
\[
\hat\rho \succeq \eta I,\quad 0<\eta<1/N,
\]
with standard choice \(\eta = 10^{-3}\), to avoid near-zero probabilities and Bures-metric blow-ups.
This is enforced in the CVXR formulation.
