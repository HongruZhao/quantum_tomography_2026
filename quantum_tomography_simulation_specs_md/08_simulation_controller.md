# Simulation controller: experiments, parameters, and outputs

This file specifies how to orchestrate the simulation study and produce the 15 plots described in the PDF.

We have:

- Systems: 1 qubit (N=2), 2 qubits (N=4)
- Libraries: 3 one-qubit + 2 two-qubit = 5 libraries
- Losses: Frobenius, Bures, Observable-expectation = 3
- Policies: Uniform, Exact, GI1 = 3

Total experiment panels: (3 + 2) × 3 = 15 panels.
Each panel has 3 curves (uniform/exact/GI1).

---

## 1) Simulation parameters (global)

Define a single list that controls everything:

```r
sim_cfg <- list(
  # Monte Carlo
  n_rep = 200,          # number of replicates (choose 50–500 depending on compute)
  seed_base = 1,        # base RNG seed
  case_seeds = NULL,    # optional named list mapping case_id -> fixed seed

  # sample sizes
  n_total = 500,        # total shots per replicate (trajectory length)
  check_every = 1,      # recompute MLE/metrics every this many steps (can be >1)

  # stabilization / conditioning
  eta_mle = 1e-4,       # eigenvalue floor for MLE (must satisfy eta < 1/N)
  ridge_init = 1e-6,    # ridge used before invertible (optional)

  # CVXR options
  solver = "SCS",
  eps_log = 1e-12,
  verbose_solver = FALSE
)
```

### Suggested defaults and scaling
- For N=2, eta_mle can be ~1e-4 to 1e-3.
- For N=4, eta_mle must be < 0.25; 1e-4 is safe.
- `n_total=500` is a reasonable starting point; increase to 1e3 if compute allows.
- `n_rep=200` yields smooth mean curves; reduce if the CVXR loop is slow.

---

## 2) True-state generation (full rank)

Each replicate requires a true state \(\rho^\star\succ0\).
We need a generator that:
- outputs `rho_true` and `theta_true` (in the chosen basis),
- respects eigenvalue floor optionally (for consistency with the stabilized MLE assumptions).

### Recommended generator: random density matrix + mixing with identity
Algorithm:
1. Draw a random complex Ginibre matrix \(G\in\mathbb C^{N\times N}\) with i.i.d. entries.
2. Form \(W = GG^\dagger\) (Wishart PSD).
3. Normalize: \(\rho_0 = W/\mathrm{tr}(W)\).
4. Mix with identity: \(\rho^\star = (1-\epsilon)\rho_0 + \epsilon I/N\) with small \(\epsilon\) (e.g. 0.05).
   This ensures full rank and prevents extremely small eigenvalues.

R skeleton:

```r
random_density_fullrank <- function(N, eps_mix = 0.05) {
  G <- matrix(rnorm(N*N) + 1i*rnorm(N*N), N, N)
  W <- G %*% Conj(t(G))
  rho0 <- W / Re(sum(diag(W)))
  rho <- (1-eps_mix)*rho0 + eps_mix*diag(N)/N
  hermitianize(rho)
}
```

Then compute theta via `theta_from_rho(rho, sigmas)`.

---

## 3) Experiment grid (systems × libraries × losses)

Define a registry of libraries with constructors:

### One-qubit
- L1: `build_library_1q_L1()`
- L2: `build_library_1q_L2()`
- L3: `build_library_1q_L3(seed = fixed)`

### Two-qubit
- A:  `build_library_2q_A()`
- B:  `build_library_2q_B()`

Each library object contains:
- `N`, `k`, `r_vec`, `Q_list`, `ab_df`, `ab_row`, `setting_labels`.

Also define basis per system:
- 1q: `build_pauli_basis_1q()$sigmas`
- 2q: `build_pauli_product_basis_2q()$sigmas`

---

## 4) Loss metric setup per experiment

For each (system, library, loss), define `metric_fun(theta_hat)`:

### Frobenius
- `G` constant: `Gf <- metric_frobenius(sigmas)`
- `metric_fun <- function(theta) Gf`

### Bures
- `metric_fun <- function(theta) metric_bures(theta, sigmas, N)`

### Observable-expectation
- Choose observable family `Obs` for the system (see 06_fisher_information_and_metrics.md).
- `Go <- metric_observable(sigmas, Obs)`
- `metric_fun <- function(theta) Go`

---

## 5) Per-replicate simulation function (single trajectory)

Define a function that runs one replicate and returns a risk curve:

```r
run_one_replicate <- function(
  N, sigmas,
  library_obj,
  loss_name,
  metric_fun,
  policy,                 # "uniform"|"exact"|"GI1"
  sim_cfg
) {
  # 1) sample true state
  rho_true <- random_density_fullrank(N)
  theta_true <- theta_from_rho(rho_true, sigmas)

  # 2) precompute S_ab and c_ab for this library
  sc <- build_Sab_cab(sigmas, library_obj$Q_list, N, library_obj$ab_df)
  S_ab <- sc$S_ab; c_ab <- sc$c_ab

  # 3) run adaptive design sequence with sampling from rho_true
  out <- adaptive_design_sequence_metric(
    theta_init = rep(0, length(sigmas)),  # or any feasible start
    rho_sampling = rho_true,
    Q_list = library_obj$Q_list, S_ab = S_ab, c_ab = c_ab, ab_df = library_obj$ab_df,
    n_total = sim_cfg$n_total,
    metric_fun = metric_fun,
    policy = policy,
    ridge_init = sim_cfg$ridge_init,
    eta_mle = sim_cfg$eta_mle,
    solver = sim_cfg$solver,
    check_every = sim_cfg$check_every,
    seed = NULL
  )
  out$risk  # numeric length n_total (or include NA for early n)
}
```

---

## 6) Monte Carlo aggregation

For each panel (system, library, loss):
1. Run `n_rep` replicates for each of the 3 policies.
2. Compute Monte Carlo mean risk curve:
   \[
   \bar R_n = \frac1{n_{\mathrm{rep}}}\sum_{r=1}^{n_{\mathrm{rep}}} R_{n}^{(r)}.
   \]
3. Store results in a tidy data frame with columns:

- `system` ("1q"|"2q")
- `library` ("L1"|"L2"|"L3"|"A"|"B")
- `loss` ("frobenius"|"bures"|"observable")
- `policy` ("uniform"|"exact"|"GI1")
- `n` (1..n_total)
- `risk_mean`
- optionally `risk_sd`, `risk_se`

Example:

```r
results <- data.frame()
for (policy in c("uniform","exact","GI1")) {
  Rmat <- replicate(sim_cfg$n_rep, run_one_replicate(..., policy=policy, ...))
  mean_curve <- rowMeans(Rmat)
  results <- rbind(results,
    data.frame(system=..., library=..., loss=..., policy=policy,
               n=1:sim_cfg$n_total, risk_mean=mean_curve))
}
```

Use `set.seed(sim_cfg$seed_base + panel_id)` for reproducibility.

---

## 7) Plotting (15 panels)

Each panel is a single (system, library, loss) combination.
Plot x-axis = n, y-axis = Monte Carlo mean risk \(\bar R_n(G)\), with three curves.

Recommended conventions:
- same y-scale within a panel,
- different panels can have different y-scales,
- include legend for policy,
- add subtitle with (system, library, loss).

If using ggplot2:

```r
ggplot(panel_df, aes(x=n, y=risk_mean, linetype=policy)) +
  geom_line() + theme_bw() +
  labs(title=..., x="n", y="mean proxy risk")
```

---

## 8) Performance notes (important)

### CVXR inside sequential loop can be expensive
Full re-optimization at every n can be slow.

Mitigations:
- Use `check_every > 1` and only refit MLE and update metrics every few steps.
- Between refits, keep \(\hat\theta\) fixed and update counts, Jn using the fixed Ia(theta_hat).
  (This is an approximation; acceptable for speed.)

### Potential alternative
Fit MLE at a grid of sample sizes (e.g., every 10 shots), compute risk only at those points, and interpolate.

For the first implementation, keep the literal approach and only then optimize.

---

## 9) Required outputs

At minimum, the simulation should write:

- `results.csv` (tidy format)  
- `plots/` directory with 15 image files (PDF/PNG), named e.g.:
  - `1q_L1_frobenius.pdf`
  - `2q_B_bures.pdf`
  etc.

Also save the simulation config `sim_cfg` for reproducibility.
