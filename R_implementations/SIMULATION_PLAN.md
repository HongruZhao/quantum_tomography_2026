# Quantum Tomography Simulation Plan

## Overview
This plan describes how to run the full adaptive quantum tomography simulation (15 cases × 3 policies) using the current R implementation. It reflects the deterministic initialization, oracle-risk plotting, and the shared plot window specified in `quantum_tomography_tutorial.html` and the current R code.

---

## 1. Simulation Parameters (current defaults in `10_run_simulation.R`)

### 1.1 Monte Carlo Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `n_rep` | 10 | Number of Monte Carlo replicates per case |
| `seed_base` | 2026 | Base random seed for reproducibility |
| `case_seeds` | NULL | Optional named list mapping `case_id` → fixed seed |
| `system_seeds` | list(`1q`=2026, `2q`=2027) | Optional named list mapping system → fixed seed for true state |
| `true_state_alpha` | 0.99 | Mix weight for pure state vs maximally mixed |

### 1.2 Sample Size Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `n_total` | 100 | Total samples per trajectory |
| `check_every` | 1 | Recompute MLE every sample (most accurate, slowest) |

### 1.3 Numerical Stabilization Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `eta_mle` | 1e-4 | Eigenvalue floor for MLE constraint (ρ ≥ ηI) |
| `ridge_init` | 1e-8 | Ridge used only during initialization (selection + oracle risk) |

### 1.4 Solver Options
| Parameter | Value | Description |
|-----------|-------|-------------|
| `solver` | "SCS" | CVXR solver |
| `verbose_solver` | FALSE | Suppress solver output |

### 1.5 Plot Window
| Parameter | Value | Description |
|-----------|-------|-------------|
| `plot_start` | 20 | Start of shared plot window |
| `plot_end` | 100 | End of shared plot window |

### 1.6 Random Seed Policy (order-independent)
- Each case uses a stable `case_id` string: `"{system}_{library}_{loss}"` (e.g., `1q_L1_frobenius`).
- **True state (non-random)**:
  - 1‑qubit: fixed pure state \(|+\rangle\langle +|\) (Bloch vector \((1,0,0)\)).
  - 2‑qubit: fixed pure Bell state \(|\Phi^+\rangle\langle \Phi^+|\).
  - `rho_true = alpha * rho_pure + (1 - alpha) * I/N`, using `true_state_alpha`.
- **Replicate seeds** (per policy):
  - `rep_base_seed = case_seeds[[case_id]]` if provided, otherwise `seed_base`.
  - `stable_seed(case_id + "_" + policy + "_rep_" + r, rep_base_seed)`.
- L3 library axes are generated with a fixed seed (`seed=42`) and **do not** perturb global RNG state.

---

## 2. The 15 Experiment Cases

### 2.1 One-Qubit Systems (N=2, d=3)

| Case | Library | Loss | k (settings) | Description |
|------|---------|------|--------------|-------------|
| 1 | L1 | Frobenius | 3 | Pauli PVMs (X,Y,Z) + Hilbert–Schmidt loss |
| 2 | L1 | Bures | 3 | Pauli PVMs + Bures loss |
| 3 | L1 | Observable | 3 | Pauli PVMs + Observable-expectation loss |
| 4 | L2 | Frobenius | 9 | Nine-axis PVMs + Hilbert–Schmidt loss |
| 5 | L2 | Bures | 9 | Nine-axis PVMs + Bures loss |
| 6 | L2 | Observable | 9 | Nine-axis PVMs + Observable-expectation loss |
| 7 | L3 | Frobenius | 4 | Random 4-axis PVMs + Hilbert–Schmidt loss |
| 8 | L3 | Bures | 4 | Random 4-axis PVMs + Bures loss |
| 9 | L3 | Observable | 4 | Random 4-axis PVMs + Observable-expectation loss |

### 2.2 Two-Qubit Systems (N=4, d=15)

| Case | Library | Loss | k (settings) | r (outcomes) | Description |
|------|---------|------|--------------|--------------|-------------|
| 10 | A | Frobenius | 15 | 2 | 15 Pauli-parity PVMs + Hilbert–Schmidt loss |
| 11 | A | Bures | 15 | 2 | 15 Pauli-parity PVMs + Bures loss |
| 12 | A | Observable | 15 | 2 | 15 Pauli-parity PVMs + Observable-expectation loss |
| 13 | B | Frobenius | 5 | 4 | 5 MUB bases + Hilbert–Schmidt loss |
| 14 | B | Bures | 5 | 4 | 5 MUB bases + Bures loss |
| 15 | B | Observable | 5 | 4 | 5 MUB bases + Observable-expectation loss |

---

## 3. Policies Compared

Each case compares three policies:

1. **Uniform Sampling**: random selection among settings.
2. **Exact One-Step Lookahead**: minimize \(\mathrm{tr}(G(\hat\theta_n)(J_n + I_a)^{-1})\).
3. **GI1 (First-Order Proxy)**: maximize \(\mathrm{tr}(G(\hat\theta_n) J_n^{-1} I_a J_n^{-1})\).

---

## 4. Initialization (Deterministic, Fixed)

Initialization is **deterministic** to make \(J_n\) invertible before adaptive selection:

- **L1**: X, Y, Z (3 samples)
- **L2**: X, Y, Z (first 3 of 9) (3 samples)
- **L3**: all 4 random axes (fixed seed) (4 samples)
- **A**: all 15 Pauli-parity PVMs (15 samples)
- **B**: all 5 MUB bases (5 samples)

---

## 5. Oracle Risk vs Adaptive Selection

- **Adaptive selection** uses MLE plug-in quantities \(J_n(\hat\theta_n)\) and \(G(\hat\theta_n)\).
- **Plotting** uses **oracle risk** at the true parameter \(\theta^*\):
  \[R_n^{\text{oracle}}(G) = \mathrm{tr}(G(\theta^*) J_n(\theta^*)^{-1})\]
- Plots show **scaled** oracle risk: \(n \cdot R_n^{\text{oracle}}(G)\).
- **True state is fixed per system (1q vs 2q)** across all replicates, libraries, losses, and policies, so Monte Carlo variation comes only from measurement outcomes and MLE noise.

---

## 6. Outputs

All outputs are written under `R_implementations/results/`:

```
R_implementations/
├── results/
│   ├── plots/
│   │   ├── case01_1q_L1_frobenius.png
│   │   ├── case02_1q_L1_bures.png
│   │   ├── ...
│   │   ├── case15_2q_B_observable.png
│   │   └── all_15_cases_combined.png
│   ├── results.csv
│   ├── simulation_config.rds
│   └── full_results.rds
```

---

## 7. How to Run

### Option A: Full 15‑Case Simulation
```r
setwd("/Users/arthur/myfile/Research simulation/Quantum_Xiaxuan/R_implementations")
source("10_run_simulation.R")
```
This runs all 15 cases and generates all plots.

### Option B: Single Case (example)
```r
setwd("/Users/arthur/myfile/Research simulation/Quantum_Xiaxuan/R_implementations")
source("10_run_simulation.R")

exp <- list(system = "2q", lib = "A", loss = "bures", case_num = 11)
res <- run_experiment(exp, simulation_config)
plot_experiment(res, output_dir = "results/plots", plot_start = 20, plot_end = 100)
```

---

## 8. Expected Plot Content

Each plot shows:
- **X-axis**: number of samples \(n\), plotted only for \(n=20..100\)
- **Y-axis**: scaled oracle risk \(n \cdot R_n^{\text{oracle}}(G)\)
- **Three curves**: Uniform (gray), Exact (blue), GI1 (red dashed)

---

## 9. Runtime Notes

- `check_every = 1` is slow but most accurate.
- For faster runs, increase `check_every` (e.g., 5 or 10) or reduce `n_rep`.

---

**Status**: READY TO RUN
